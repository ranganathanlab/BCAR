# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True
cimport cython
from libc.math cimport log10
from cython.view cimport array as cvarray

def build_consensus(list sequences, int max_qscore=40, double min_miscall_rate=0.01):
    cdef int L, i, j, b, max_calls, total, qscore, idx, call
    cdef int most_frequent_base = -1
    cdef double miscall_rate, posterior_prob, likelihood_star, sum_likelihoods
    cdef int[:, :] base_counts
    cdef int[:] gap_counts
    cdef double[::1] likelihoods
    cdef int[::1] consensus_seq
    cdef int[::1] qscores
    cdef dict base_equivalence = {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 1, 6: 2, 7: 3, 8: 4, 9: 9}
    
    # Determine the maximum sequence length
    L = max(len(seq) for seq in sequences)
    
    # Initialize base and gap counts
    base_counts = cvarray(shape=(L, 5), itemsize=sizeof(int), format="i", mode="c")
    gap_counts = cvarray(shape=(L,), itemsize=sizeof(int), format="i", mode="c")
    base_counts[:, :] = 0
    gap_counts[:] = 0
    
    # Equalize sequence lengths by padding gaps
    for seq in sequences:
        while len(seq) < L:
            seq.append(0)  # Pad with gaps
    
    # Count bases and gaps at each position
    for i in range(L):
        for seq in sequences:
            call = base_equivalence[seq[i]]
            if call == 0:
                gap_counts[i] += 1
            else:
                base_counts[i, call] += 1
    
    # Estimate total miscall rate across the sequence
    total_miscalls = 0
    total_calls = 0
    for i in range(L):
        total_calls += sum(base_counts[i, :])
        total_miscalls += sum(base_counts[i, :]) - max(base_counts[i, :])
    miscall_rate = max(total_miscalls / total_calls, min_miscall_rate)
    
    # Allocate consensus sequence and Q-scores
    consensus_seq = cvarray(shape=(L,), itemsize=sizeof(int), format="i", mode="c")
    qscores = cvarray(shape=(L,), itemsize=sizeof(int), format="i", mode="c")
    consensus_seq[:] = 0
    qscores[:] = 0

    # Iterate over each position to compute consensus and Q-scores
    idx = 0
    for i in range(L):
        max_calls = max(base_counts[i, :])
        
        # Skip positions where gaps are the majority
        if gap_counts[i] > max_calls:
            continue
        
        # Check for ambiguity
        if list(base_counts[i, :]).count(max_calls) > 1:
            consensus_seq[idx] = 9  # "N"
            qscores[idx] = 0
            idx += 1
            continue
        
        # Determine the most frequent base
        for b in range(1, 5):
            if base_counts[i, b] == max_calls:
                most_frequent_base = b
        
        # Compute Bayesian posterior probability
        total = sum(base_counts[i, :])
        likelihoods = cvarray(shape=(4,), itemsize=sizeof(double), format="d", mode="c")
        for b in range(1, 5):
            likelihoods[b - 1] = (1 - miscall_rate) ** base_counts[i, b] * \
                                 (miscall_rate / 3) ** (total - base_counts[i, b])
        sum_likelihoods = sum(likelihoods)
        posterior_prob = max(likelihoods) / sum_likelihoods if sum_likelihoods > 0 else 0
        
        # Compute Q-score
        if posterior_prob < 1:
            qscore = int(-10 * log10(1 - posterior_prob))
        else:
            qscore = max_qscore
        qscore = min(qscore, max_qscore)
        
        # Append results
        consensus_seq[idx] = most_frequent_base
        qscores[idx] = qscore
        idx += 1
    
    # Resize arrays to final length
    consensus_seq = consensus_seq[:idx]
    qscores = qscores[:idx]

    return list(consensus_seq), list(qscores)

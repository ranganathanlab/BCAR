# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True
cimport cython
from libc.stdlib cimport malloc, free
from libc.string cimport memset

@cython.cfunc
cdef int max3(int a, int b, int c):
    return max(a, max(b, c))

@cython.cfunc
cdef int get_score(int b1, int b2, int match, int mismatch, int gap_penalty):
    if b1 == 0 or b2 == 0:  # Gap
        return -gap_penalty
    elif b1 == b2:  # Match
        return match
    else:  # Mismatch
        return -mismatch

@cython.boundscheck(False)
@cython.wraparound(False)
def needleman_wunsch(
    list query, list reference,
    int match=2, int mismatch=1, int gap_penalty=2
):
    cdef int q_len = len(query)
    cdef int r_len = len(reference)
    cdef int *score_matrix = <int*>malloc((q_len + 1) * (r_len + 1) * sizeof(int))
    cdef int *trace_matrix = <int*>malloc((q_len + 1) * (r_len + 1) * sizeof(int))
    cdef int i, j, index, score_diag, score_up, score_left
    
    if score_matrix is NULL or trace_matrix is NULL:
        raise MemoryError("Unable to allocate memory for matrices.")
    
    memset(score_matrix, 0, (q_len + 1) * (r_len + 1) * sizeof(int))
    memset(trace_matrix, 0, (q_len + 1) * (r_len + 1) * sizeof(int))
    
    # Initialization
    for i in range(q_len + 1):
        score_matrix[i * (r_len + 1)] = -i * gap_penalty
    for j in range(r_len + 1):
        score_matrix[j] = -j * gap_penalty
    
    # Fill score and traceback matrices
    for i in range(1, q_len + 1):
        for j in range(1, r_len + 1):
            index = i * (r_len + 1) + j
            score_diag = score_matrix[(i - 1) * (r_len + 1) + (j - 1)] + get_score(query[i - 1], reference[j - 1], match, mismatch, gap_penalty)
            score_up = score_matrix[(i - 1) * (r_len + 1) + j] - gap_penalty
            score_left = score_matrix[i * (r_len + 1) + (j - 1)] - gap_penalty
            score_matrix[index] = max3(score_diag, score_up, score_left)
            if score_matrix[index] == score_diag:
                trace_matrix[index] = 1  # Diagonal
            elif score_matrix[index] == score_up:
                trace_matrix[index] = 2  # Up
            else:
                trace_matrix[index] = 3  # Left
    
    # Traceback
    cdef list aligned_query = []
    i, j = q_len, r_len
    while i > 0 or j > 0:
        if i > 0 and j > 0 and trace_matrix[i * (r_len + 1) + j] == 1:  # Diagonal
            aligned_query.append(query[i - 1])
            i -= 1
            j -= 1
        elif i > 0 and (j == 0 or trace_matrix[i * (r_len + 1) + j] == 2):  # Up
            aligned_query.append(0)  # Gap in reference
            i -= 1
        else:  # Left
            aligned_query.append(reference[j - 1] + 4)  # Insertion
            j -= 1
    
    aligned_query.reverse()
    free(score_matrix)
    free(trace_matrix)
    return aligned_query

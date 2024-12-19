# BCAR: Barcode Collapse by Aligning Reads

BCAR is a Python tool designed to process raw high-throughput sequencing data. Its primary purpose is to generate consensus reads for each barcode present in the data, eliminating idiosyncratic sequencing errors. This ensures accurate representation of the sequences associated with each barcode.

## Features
- Uses a fast implementation of the Needleman-Wunsch alignment algorithm to handle indels (insertions and deletions) between reads.
- Supports scenarios where barcodes are located at fixed positions within reads.
- Accepts any number of input .fastq files

## Requirements
- Current version only supports paired reads
- The barcode must occur at a **fixed location** in the **forward** read.
- If the barcode is on the reverse read, simply pass the reverse read to the --fwd option and vice versa. 
- If the barcode position is not fixed but is adjacent to a constant sequence, you can preprocess your data using tools like [CutAdapt](https://cutadapt.readthedocs.io/) to trim your reads, ensuring that the barcodes are at a fixed location.

## Installation

### 1. Clone the Repository
```bash
git clone https://github.com/yourusername/BCAR.git
cd BCAR
```

### 2. Build and Install
BCAR relies on a Cython implementation of the Needleman-Wunsch alignment algorithm. Follow these steps to build and install:

```bash
python setup.py build
python setup.py install
```

Ensure that you have `Cython` installed in your Python environment:
```bash
pip install cython
```

## Usage
Once installed, BCAR can be invoked as a command-line script. The typical usage is:

```bash
bcar [-h] --fwd FWD_FASTQS [FWD_FASTQS ...] --rev REV_FASTQS [REV_FASTQS ...] [--out1 OUTPUT_FASTQ_FWD] [--out2 OUTPUT_FASTQ_REV] [--BC-start BC_START] [--BC-len BC_LEN] [--min-qscore MIN_QSCORE] [--min-agree MIN_AGREEMENT] [--min-count MIN_COUNT]
```

### Options:
- `--fwd FWD_FASTQS [FWD_FASTQS ...]`: List of `.fastq` or `.fastq.gz` files containing forward reads.
- `--rev REV_FASTQS [REV_FASTQS ...]`: List of `.fastq` or `.fastq.gz` files containing reverse reads.
- `--out1 OUTPUT_FASTQ_FWD`: Output `.fastq` file for the consensus forward reads.
- `--out2 OUTPUT_FASTQ_REV`: Output `.fastq` file for the consensus reverse reads.
- `--BC-start BC_START`: Position in the sequence where the barcode begins.
- `--BC-len BC_LEN`: Length of the barcode.
- `--min-qscore MIN_QSCORE`: Minimum quality score for bases in the barcode or mapped region; reads below this threshold are discarded.
- `--min-agree MIN_AGREEMENT`: Fraction of bases at each position in the consensus that must agree, otherwise the barcode is discarded.
- `--min-count MIN_COUNT`: Minimum number of times a barcode must appear to be included in the output.

### Example Workflow
1. **Pre-trim Reads (if necessary)**:
   Use CutAdapt or another preprocessing tool to ensure that barcodes occur at a fixed location.
   ```bash
   cutadapt -g ^<constant_sequence> -o trimmed_reads.fastq input_reads.fastq
   ```

2. **Run BCAR**:
   ```bash
   bcar --fwd forward_reads_1.fastq.gz forward_reads_2.fastq.gz \
        --rev reverse_reads_1.fastq.gz reverse_reads_2.fastq.gz \
        --out1 consensus_forward.fastq --out2 consensus_reverse.fastq \
        --BC-start 10 --BC-len 16 --min-qscore 20 --min-agree 0.8 --min-count 5
   ```

## How BCAR Works
BCAR leverages a fast Needleman-Wunsch alignment implementation to align reads associated with the same barcode. This accounts for sequencing errors such as substitutions, insertions, and deletions, ensuring accurate consensus reads.

## Contributing
Contributions are welcome! Feel free to open issues or submit pull requests on the GitHub repository.

## License
This project is licensed under the MIT License. See the `LICENSE` file for details.

---

For further information, contact [Your Name](mailto:your.email@example.com).



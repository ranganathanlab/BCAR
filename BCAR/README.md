# BCAR: Barcode Collapse by Aligning Reads

BCAR is a Python tool designed to process raw high-throughput sequencing data. Its primary purpose is to generate consensus reads for each barcode present in the data, eliminating idiosyncratic sequencing errors. This ensures accurate representation of the sequences associated with each barcode.

## Features
- Processes paired-end sequencing data.
- Uses a fast implementation of the Needleman-Wunsch alignment algorithm to handle indels (insertions and deletions) between reads.
- Supports scenarios where barcodes are located at fixed positions within reads.

## Requirements
- The barcode must occur at a **fixed location** in the read.
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
bcar --input <input_file> --output <output_file> --barcode-position <position>
```

### Options:
- `--input`: Path to the input paired-end read file (in FASTQ format).
- `--output`: Path to the output file where consensus reads will be written.
- `--barcode-position`: The fixed position of the barcode within the reads.

### Example Workflow
1. **Pre-trim Reads (if necessary)**:
   Use CutAdapt or another preprocessing tool to ensure that barcodes occur at a fixed location.
   ```bash
   cutadapt -g ^<constant_sequence> -o trimmed_reads.fastq input_reads.fastq
   ```

2. **Run BCAR**:
   ```bash
   bcar --input trimmed_reads.fastq --output consensus_reads.fastq --barcode-position 10
   ```

## How BCAR Works
BCAR leverages a fast Needleman-Wunsch alignment implementation to align reads associated with the same barcode. This accounts for sequencing errors such as substitutions, insertions, and deletions, ensuring accurate consensus reads.

## Contributing
Contributions are welcome! Feel free to open issues or submit pull requests on the GitHub repository.

## License
This project is licensed under the MIT License. See the `LICENSE` file for details.

---

For further information, contact [Your Name](mailto:your.email@example.com).



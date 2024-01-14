# TrimAdapt Tool
Welcome to the **README Documentation and User Manual** for the TrimAdapt Tool. This comprehensive guide is designed to assist you in understanding and effectively utilizing the tool for your needs.

To get started quickly, we recommend you to first navigate to the 'INTRODUCTION' and 'INSTALL' section. Here, you will find detailed instructions on how to properly install the TrimAdapt Tool on your system. The instructions are designed to be easy to follow, ensuring a smooth installation process.

Once the installation is complete, proceed to the 'USE' section. This section provides a concise tutorial on how to use the TrimAdapt Tool. It covers the basic operations and functionalities of the tool, enabling you to start using it effectively right away.

Remember, this documentation is here to help you. Whether you're a first-time user or an experienced professional, this user manual is designed to make your experience with the TrimAdapt Tool as seamless and productive as possible. 

# INTRODUCTION
The TrimAdapt tool, redesigned and re-implemented to be user-friendly, is a powerful tool for the preprocessing of Next-Generation Sequencing (NGS) data. It efficiently trims adapter sequences from both single-end and paired-end FASTQ files, a common requirement in NGS data preprocessing.

TrimAdapt, similar to Trimmomatic, utilizes a command-line interface for trimming adapters from single-end and paired-end FASTQ files in sequencing data. It is implemented using Python, ensuring error-free operation and enhanced accuracy of adapter identification.

The precision of adapter detection has been improved by refining the detection algorithms. This enhancement increases the tool’s efficiency and reliability, ensuring accurate identification and removal of adapters in complex data sets.

TrimAdapt offers a wide range of customizable trimming parameters, allowing fine-tuning of the trimming process according to specific needs. This provides users with greater control and flexibility.

The tool also implements statistical reporting of the trimming results, providing a comprehensive understanding of the trimming process and its outcomes. Detailed statistics about the number and types of adapters removed, the length of sequences trimmed, and other relevant information are readily available.

The results of the trimming process are visualized through bar and pie charts, providing a clear and intuitive representation of the data. This aids in result interpretation and decision-making.

Furthermore, a comprehensive evaluation has confirmed that the performance of the re-implemented TrimAdapt tool is comparable to that of Trimmomatic, ensuring that the tool is not only user-friendly but also efficient and reliable.

## Features
- **Supports both single-end and paired-end FASTQ files**: The TrimAdapt tool can process both single-end and paired-end data.
- **Customizable adapter trimming**: Users can specify any number of adapter sequences to be trimmed.
- **Flexible alignment scoring**: Match, mismatch, gap opening, and gap extension scores can be modified according to the user's needs.
- **Auto-detect gzipped files**: If the input FASTQ files are gzipped, the script will automatically handle them.

## Dependencies
- Python 3.x
- Biopython
- gzip (for gzipped FASTQ files)
- matplotlib (for optional plotting functionality)

> # Note
>
> ## Python Standard Library Packages
>
> The following packages are part of the standard Python library and do not require separate installation:
>
> - `os`
> - `sys`
> - `gzip`
> - `argparse`
>
> These packages provide a wide range of functionality, from interacting with the operating system (`os`, `sys`), to handling compressed files (`gzip`), to parsing command-line arguments (`argparse`). They are readily available in any standard Python environment.

Ensure all dependencies are installed in your Python environment. If not, they can typically be installed via pip:

```bash
pip install biopython matplotlib
```

## Installation
1. Clone or download the script from the repository.
2. Ensure all dependencies are installed.
3. The script is ready to be executed from the command line.
4. TrimAdapt can be download  as a Docker image that comes with all the dependencies pre-installed. This makes it easy for the user to run the tool without any hassle. To download the Docker image, use the following command (Detailed instructions on how to pull and use the Docker image will be provided in the below after CLI sections.)

## Usage
The script is executed from the command line with the following syntax:

```bash
python TrimAdapt.py [arguments]
```

### Required Arguments
- `--data_type`: Specify 'single' for single-end data or 'paired' for paired-end data.
- `--input_file`: Path to the input FASTQ file (for single-end data).
- `--output_file`: Path for the output trimmed FASTQ file (for single-end data).
- `--adapters`: Comma-separated adapter sequences or a path to a FASTA file containing adapter sequences.

### Optional Arguments
- `--match_score`: Score for matching bases (default: 2).
- `--mismatch_score`: Penalty for mismatching bases (default: -2).
- `--open_gap_score`: Penalty for opening a gap (default: -2).
- `--extend_gap_score`: Penalty for extending a gap (default: -0.1).
- For paired-end data, additional `--input_file1`, `--input_file2`, `--output_file1`, and `--output_file2` arguments are required.

### Running the Script
Here is an example of running the script for single-end data:

```bash
python TrimAdapt.py --data_type single --input_file sample.fastq --output_file trimmed_sample.fastq --adapters "ACTG,GTCA" --match_score 1 --mismatch_score -1 --open_gap_score -1 --extend_gap_score -0.5
```

For paired-end data, specify both input and output files for each end:

```bash
python TrimAdapt.py --data_type paired --input_file1 sample_R1.fastq --input_file2 sample_R2.fastq --output_file1 trimmed_R1.fastq --output_file2 trimmed_R2.fastq --adapters "ACTG,GTCA"
```

## Output
The script generates trimmed FASTQ files as specified in the `--output_file` argument(s). If plotting is enabled and matplotlib is installed, it also generates a bar plot and pie chart summarizing the trimming results.

## Pulling the Docker Image
Before running the container, you will need to pull it from the Docker registry. Use the following command to pull the TrimAdapt Docker image:

```sh
docker pull qaswara98/ubuntu:TrimAdapt
```

## Usage (Docker image)

### Data Preparation
Ensure that your FASTQ files and any adapter sequence files are in a known directory on your host system. You will need to mount this directory to the Docker container to make the files accessible to TrimAdapt.

### Single-End Data Example
To run the TrimAdapt tool for single-end data, use the following command structure:

```sh
docker run -v /path/to/your/data:/data qaswara98/ubuntu:TrimAdapt python3 TrimAdapt.py \
  --data_type single \
  --input_file /data/adapter_content.fastq \
  --output_file /data/trimmed_sample.fastq \
  --adapters "ADAPTER_SEQUENCE" \
  --match_score 2 \
  --mismatch_score -1 \
  --open_gap_score -2 \
  --extend_gap_score -1
```

Replace `/path/to/your/data` with the actual path to the directory containing your FASTQ file and adapter sequences. Replace `"ADAPTER_SEQUENCE"` with your actual adapter sequence.

### Paired-End Data Example
For paired-end data, use the following command structure:

```sh
docker run -v /path/to/your/data:/data qaswara98/ubuntu:TrimAdapt python3 TrimAdapt.py \
  --data_type paired \
  --input_file1 /data/sample1_forward.fastq \
  --input_file2 /data/sample1_reverse.fastq \
  --output_file1 /data/trimmed_forward.fastq \
  --output_file2 /data/trimmed_reverse.fastq \
  --adapters "ADAPTER_SEQUENCE" \
  --match_score 2 \
  --mismatch_score -1 \
  --open_gap_score -2 \
  --extend_gap_score -1
```

Again, replace the placeholders with paths and values specific to your data.

## Notes

- Ensure that you use absolute paths when specifying file locations in the Docker command.
- The adapter sequence and alignment scores (`--match_score`, `--mismatch_score`, `--open_gap_score`, `--extend_gap_score`) should be adjusted according to your specific requirements.
- The `-v /path/to/your/data:/data` part of the command mounts your data directory (`/path/to/your/data`) to the `/data` directory in the container. Ensure that paths inside the container start with `/data/` to correctly point to these files.

# Best Practices
- Always test the script with a small subset of your data first to ensure it performs as expected.
- Backup important data before processing.
- Regularly update the script and dependencies to incorporate improvements and bug fixes.

# Troubleshooting
- Ensure all paths (to input files, output files, adapter files) are correct and accessible by the script.
- Ensure all dependencies are correctly installed.
- If you encounter errors, check the console output for hints about what might be wrong, adjust parameters or inputs accordingly, and try again.

# License
The script is released under the GPL-3.0 license, as typically included in the repository.
# FINAL NOTES
If you decide to test this program, I would appreciate hearing about your experience. Please share any feedback or issues you encounter. Additionally, don’t hesitate to adapt the code to better suit your needs. Your input is valuable and can help improve this program.

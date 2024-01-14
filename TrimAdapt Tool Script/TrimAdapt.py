"""
This TrimAdapt tool script is designed to trim adapter sequences from FASTQ files, a common task in bioinformatics
pre-processing of NGS (Next Generation Sequencing) data. It supports both single-end and paired-end
FASTQ files and allows for customizable adapter sequences and alignment scoring.
"""

import os
import sys
import gzip
import matplotlib.pyplot as plt
from Bio import SeqIO, Seq, Align
import argparse

# Define default alignment scores
DEFAULT_MATCH_SCORE = 2
DEFAULT_MISMATCH_SCORE = -2
DEFAULT_OPEN_GAP_SCORE = -2
DEFAULT_EXTEND_GAP_SCORE = -0.1

# Define thresholds for clipping palindromic and simple sequences
PALINDROME_CLIP_THRESHOLD = 30
SIMPLE_CLIP_THRESHOLD = 10

def manual_trim(sequence, adapter):
    """
    Trims a sequence at the first occurrence of an adapter sequence.

    Parameters:
        sequence (str): The nucleotide sequence (read) from a FASTQ file.
        adapter (str): The adapter sequence to be trimmed off.

    Returns:
        tuple: (trimmed_sequence, found)
            trimmed_sequence (str) - The sequence after trimming the adapter.
            found (bool) - True if the adapter was found and trimmed, else False.
    """
    if adapter in sequence:
        return sequence.split(adapter)[0], True
    return sequence, False

def trim_adapters(read, adapters, match_score, mismatch_score, open_gap_score, extend_gap_score):
    """
    Trims adapter sequences from a read using pairwise alignment.
    Iterates through provided adapter sequences and trims the first identified adapter.

    Parameters:
        read (str): The original sequence read from the FASTQ file.
        adapters (list): A list of adapter sequences to be trimmed.
        match_score, mismatch_score, open_gap_score, extend_gap_score (float): Alignment scoring parameters.

    Returns:
        tuple: (trimmed_read, adapter_found)
            trimmed_read (str) - The read after trimming adapters.
            adapter_found (bool) - True if any adapter was found and trimmed, else False.
    """
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = match_score
    aligner.mismatch_score = mismatch_score
    aligner.open_gap_score = open_gap_score
    aligner.extend_gap_score = extend_gap_score

    trimmed_sequence = read
    adapter_found = False

    for adapter in adapters:
        if len(trimmed_sequence) == 0:
            break

        alignments = aligner.align(trimmed_sequence, adapter)

        if not alignments:
            trimmed_sequence, manual_found = manual_trim(trimmed_sequence, adapter)
            adapter_found = adapter_found or manual_found
            continue

        best_alignment = alignments[0]
        score = best_alignment.score
        if score >= PALINDROME_CLIP_THRESHOLD or score >= SIMPLE_CLIP_THRESHOLD:
            start_of_adapter_in_read = best_alignment.aligned[0][0][0]
            end_of_adapter_in_read = best_alignment.aligned[0][0][1]
            trimmed_sequence = trimmed_sequence[:start_of_adapter_in_read] + trimmed_sequence[end_of_adapter_in_read:]
            adapter_found = True
            break

    return trimmed_sequence, adapter_found

def process_fastq_single(input_file, output_file, adapters, match_score=DEFAULT_MATCH_SCORE, mismatch_score=DEFAULT_MISMATCH_SCORE, open_gap_score=DEFAULT_OPEN_GAP_SCORE, extend_gap_score=DEFAULT_EXTEND_GAP_SCORE):
    """
    Processes a single FASTQ file, trimming adapters from each read and writing the trimmed reads to an output file.

    Parameters:
        input_file (str): Path to the input FASTQ file.
        output_file (str): Path for the output FASTQ file with trimmed reads.
        adapters (list): Adapter sequences to trim from reads.
        match_score, mismatch_score, open_gap_score, extend_gap_score (float): Alignment scoring parameters.
    """
    adapter_found = 0
    total_reads = 0

    with (gzip.open(input_file, 'rt') if input_file.endswith('.gz') else open(input_file, 'r')) as fin:
        reads = []
        for record in SeqIO.parse(fin, "fastq"):
            total_reads += 1
            trimmed_sequence, found = trim_adapters(str(record.seq), adapters, match_score, mismatch_score, open_gap_score, extend_gap_score)

            if found:
                adapter_found += 1
                quality_scores = record.letter_annotations["phred_quality"][:len(trimmed_sequence)]
                record.letter_annotations.clear()
                record.seq = Seq.Seq(trimmed_sequence)
                record.letter_annotations["phred_quality"] = quality_scores

            reads.append(record)

        with (gzip.open(output_file, 'wt') if output_file.endswith('.gz') else open(output_file, 'w')) as fout:
            SeqIO.write(reads, fout, "fastq")

    adapter_not_found = total_reads - adapter_found
    untrimmed_reads = total_reads - adapter_found
    print(f"Total reads processed: {total_reads}")
    print(f"Reads with adapter sequences found and trimmed: {adapter_found}")
    print(f"Reads with no adapter sequences found (Untrimmed Reads): {untrimmed_reads}")

    # Define the data for the bar plot
    labels = ['Processed', 'Trimmed', 'Untrimmed Reads']
    values_bar = [total_reads, adapter_found, untrimmed_reads]
    colors = ['skyblue', 'orange', 'green']

    # Create bar plot
    plt.figure(figsize=(8, 4))
    plt.bar(labels, values_bar, color=colors)
    plt.title('Read Processing Summary (Single-End)')
    plt.ylabel('Number of Reads')
    plt.savefig('bar_plot_single_end.png')

    # Define the data for the pie chart 
    labels_pie = ['Trimmed', 'Untrimmed Reads']
    values_pie = [adapter_found, untrimmed_reads]

    # Create pie chart
    plt.figure(figsize=(8, 8))
    plt.pie(values_pie, labels=labels_pie, colors=colors[1:], autopct='%1.1f%%', startangle=140)
    plt.title('Read Processing Summary (Single-End)')
    plt.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle
    plt.savefig('pie_chart_single_end.png')

def process_fastq_paired(input_file1, input_file2, output_file1, output_file2, adapters, match_score=DEFAULT_MATCH_SCORE, mismatch_score=DEFAULT_MISMATCH_SCORE, open_gap_score=DEFAULT_OPEN_GAP_SCORE, extend_gap_score=DEFAULT_EXTEND_GAP_SCORE):
    """
    Processes paired-end FASTQ files, trimming adapters from each read pair and writing the trimmed reads to output files.

    Parameters:
        input_file1, input_file2 (str): Paths to the input FASTQ files for the paired-end data.
        output_file1, output_file2 (str): Paths for the output FASTQ files with trimmed reads.
        adapters (list): Adapter sequences to trim from reads.
        match_score, mismatch_score, open_gap_score, extend_gap_score (float): Alignment scoring parameters.
    """
    forward_adapter_found = 0
    reverse_adapter_found = 0
    total_pairs = 0

    with (gzip.open(input_file1, 'rt') if input_file1.endswith('.gz') else open(input_file1, 'r')) as fin1,\
        (gzip.open(input_file2, 'rt') if input_file2.endswith('.gz') else open(input_file2, 'r')) as fin2:

        trimmed_forward_reads = []
        trimmed_reverse_reads = []

        for forward_read, reverse_read in zip(SeqIO.parse(fin1, "fastq"), SeqIO.parse(fin2, "fastq")):
            total_pairs += 1
            trimmed_forward_sequence, found_forward = trim_adapters(str(forward_read.seq), adapters, match_score, mismatch_score, open_gap_score, extend_gap_score)
            trimmed_reverse_sequence, found_reverse = trim_adapters(str(reverse_read.seq), adapters, match_score, mismatch_score, open_gap_score, extend_gap_score)

            if found_forward or found_reverse:
                forward_adapter_found += 1
                forward_quality_scores = forward_read.letter_annotations["phred_quality"][:len(trimmed_forward_sequence)]
                forward_read.letter_annotations.clear()
                forward_read.seq = Seq.Seq(trimmed_forward_sequence)
                forward_read.letter_annotations["phred_quality"] = forward_quality_scores

                reverse_adapter_found += 1
                reverse_quality_scores = reverse_read.letter_annotations["phred_quality"][:len(trimmed_reverse_sequence)]
                reverse_read.letter_annotations.clear()
                reverse_read.seq = Seq.Seq(trimmed_reverse_sequence)
                reverse_read.letter_annotations["phred_quality"] = reverse_quality_scores

            trimmed_forward_reads.append(forward_read)
            trimmed_reverse_reads.append(reverse_read)

        with (gzip.open(output_file1, 'wt') if output_file1.endswith('.gz') else open(output_file1, 'w')) as fout1,\
            (gzip.open(output_file2, 'wt') if output_file2.endswith('.gz') else open(output_file2, 'w')) as fout2:
            SeqIO.write(trimmed_forward_reads, fout1, "fastq")
            SeqIO.write(trimmed_reverse_reads, fout2, "fastq")

    forward_adapter_not_found = total_pairs - forward_adapter_found
    reverse_adapter_not_found = total_pairs - reverse_adapter_found
    untrimmed_forward_reads = total_pairs - forward_adapter_found
    untrimmed_reverse_reads = total_pairs - reverse_adapter_found
    print(f"Total read pairs processed: {total_pairs}")
    print(f"Read pairs with adapter sequences found and trimmed (forward): {forward_adapter_found}")
    print(f"Read pairs with adapter sequences found and trimmed (reverse): {reverse_adapter_found}")
    print(f"Forward reads with no adapter sequences found (Untrimmed Reads): {untrimmed_forward_reads}")
    print(f"Reverse reads with no adapter sequences found (Untrimmed Reads): {untrimmed_reverse_reads}")
    
    
    # Define the data for the bar plot
    paired_labels_bar = ['Processed Pairs', 'Trimmed (Forward)', 'Trimmed (Reverse)', 'Untrimmed (Forward)', 'Untrimmed (Reverse)']
    paired_values_bar = [total_pairs, forward_adapter_found, reverse_adapter_found, untrimmed_forward_reads, untrimmed_reverse_reads]
    colors_paired_bar = ['skyblue', 'orange', 'green', 'red', 'purple']

    # Create bar plot for paired-end data
    plt.figure(figsize=(10, 5))
    plt.bar(paired_labels_bar, paired_values_bar, color=colors_paired_bar)
    plt.title('Read Pair Processing Summary (Paired-End)')
    plt.ylabel('Number of Read Pairs')
    plt.savefig('bar_plot_paired_end.png')

    # Define the data for the pie chart
    trimmed_pairs = forward_adapter_found + reverse_adapter_found
    untrimmed_pairs = total_pairs - trimmed_pairs

    paired_labels_pie = ['Trimmed Pairs', 'Untrimmed Pairs']
    paired_values_pie = [trimmed_pairs, untrimmed_pairs]
    colors_paired_pie = ['orange', 'green']
    
    # Create pie chart for paired-end data
    plt.figure(figsize=(8, 8))
    plt.pie(paired_values_pie, labels=paired_labels_pie, colors=colors_paired_pie, autopct='%1.1f%%', startangle=140)
    plt.title('Read Pair Processing Summary (Paired-End)')
    plt.axis('equal')  # Equal aspect ratio ensures pie is drawn as a circle
    plt.savefig('pie_chart_paired_end.png')

def validate_file(file):
    """Check if the file exists and is readable."""
    if not os.path.isfile(file):
        print(f"Error: File '{file}' not found.")
        sys.exit(1)

def main(args):
    """
    Main function to process the command line arguments and call the appropriate functions based on data type.
    Validates the provided files and then processes them based on specified data type (single or paired).

    Parameters:
        args (Namespace): Parsed command line arguments.
    """
    adapters = []
    if args.adapter_file:
        validate_file(args.adapter_file)
        for record in SeqIO.parse(args.adapter_file, "fasta"):
            adapters.append(str(record.seq))
    elif args.adapters:
        adapters = [adapter.strip() for adapter in args.adapters.split(',')]

    if args.data_type == 'single':
        validate_file(args.input_file)
        process_fastq_single(args.input_file, args.output_file, adapters, args.match_score, args.mismatch_score, args.open_gap_score, args.extend_gap_score)
    elif args.data_type == 'paired':
        validate_file(args.input_file1)
        validate_file(args.input_file2)
        process_fastq_paired(args.input_file1, args.input_file2, args.output_file1, args.output_file2, adapters, args.match_score, args.mismatch_score, args.open_gap_score, args.extend_gap_score)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Trim adapters from sequencing data.')

    # Required arguments
    parser.add_argument('--data_type', choices=['single', 'paired'], required=True, help='Type of data (single or paired).')
    adapter_group = parser.add_mutually_exclusive_group(required=True)
    adapter_group.add_argument('--adapters', help='Comma-separated adapter sequences.')
    adapter_group.add_argument('--adapter_file', help='FASTA file containing adapter sequences.')
    parser.add_argument('--input_file', help='Input FASTQ file for single-end data.')
    parser.add_argument('--output_file', help='Output FASTQ file for single-end data.')
    parser.add_argument('--input_file1', help='First input FASTQ file for paired-end data.')
    parser.add_argument('--input_file2', help='Second input FASTQ file for paired-end data.')
    parser.add_argument('--output_file1', help='First output FASTQ file for paired-end data.')
    parser.add_argument('--output_file2', help='Second output FASTQ file for paired-end data.')

    # Optional arguments for alignment scores
    parser.add_argument('--match_score', type=float, default=DEFAULT_MATCH_SCORE, help='Score for matching bases. Default is 2.')
    parser.add_argument('--mismatch_score', type=float, default=DEFAULT_MISMATCH_SCORE, help='Score for mismatching bases. Default is -2.')
    parser.add_argument('--open_gap_score', type=float, default=DEFAULT_OPEN_GAP_SCORE, help='Score for opening a gap. Default is -2.')
    parser.add_argument('--extend_gap_score', type=float, default=DEFAULT_EXTEND_GAP_SCORE, help='Score for extending a gap. Default is -0.1.')
    
    args = parser.parse_args()
    main(args)

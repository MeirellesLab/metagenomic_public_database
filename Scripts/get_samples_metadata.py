## Author: Eduardo HC Galvao - eduardohcgalvao@gmail.com
## Version: 3.0 (2023/03/07)
## Created (yy/mm/dd): 2023/03/01

## USAGE: python3 get_samples_metadata.py INPUT_PATH REFERENCE_FILE.csv

## From a given directory list the samples files and compare to the "REFERENCE_FILE.csv"
## Then run a series of functions to collect metadata: File_size, file number of reads, GC%, average read size, number of bases.

## INPUT_PATH ==> Main directory containing ".fasta" files metagenomes samples.

## REFERENCE_FILE.csv ==> ".csv" file separated by "," with the first column being the samples name ID.
## Assuming the file contains all samples of the given database.

import os, sys, shutil, math, datetime

# Assign input path as a variable.
input_path = sys.argv[1]

# Get reference file as variable.
reference_file = sys.argv[2]

# Get actual time for log file.
def update_time():
    now = datetime.datetime.now()
    current_time = now.strftime("%H:%M:%S")
    return current_time

# Open a ".csv" file containing the metadata for reference and comparison.
def get_reference_dict(reference_file):
    current_time = update_time()
    print(current_time + " - Getting metadata from reference file.\n")
    reference_dict = {}
    with open(reference_file, "r+") as file:
        lines = file.read().splitlines()
        for line in lines:
            sample, metadata_categories = line.split(";", 1)
            reference_dict[sample] = str(metadata_categories)
    return reference_dict

# Function to get the size of a file.
def get_file_size(file):
    current_time = update_time()
    print(current_time + " - Calculating " + str(file).split("_", 1)[0]  + " size.")
    try:
        file_size = os.path.getsize(file)
    except:
        file_size = "error_calculating_file_size"
    return str(file_size)

# Function to get the number of Reads of a file.
def get_file_number_of_reads(input_file):
    current_time = update_time()
    print(current_time + " - Calculating " + str(input_file).split("_", 1)[0] + " number of reads.")
    with open(input_file, "r") as file:
        count = 0
        for line in file:
            for character in line:
                if character == ">":
                    count += 1
        sample_number_of_reads = str(count)
    return sample_number_of_reads

# Function to get the GC content and number of bases.
def get_sample_GC_content(input_file):
    current_time = update_time()
    print(current_time + " - Calculating " + str(input_file).split("_", 1)[0] + " %GC content.")
    sequence = ""
    GC_count = 0
    sample_number_of_bases = 0
    with open(input_file, "r") as file:
        for line in file:
            if not line.startswith(">"):
                sequence += line.strip()
                GC_count += line.count("G") + line.count("C")
                sample_number_of_bases += line.count("G") + line.count("C") + line.count("A") + line.count("T")
    sample_GC_content = GC_count/float(len(sequence))
    return sample_GC_content, sample_number_of_bases

# Function to get the average read size from a fasta file.
def get_average_read_size(input_file):
    current_time = update_time()
    print(current_time + " - Calculating " + str(input_file).split("_", 1)[0]  + " average read size.")
    total_base_pairs = 0
    total_reads = 0
    with open(input_file, "r") as file:
        for line in file:
            if line.startswith(">"):
                total_reads += 1
            else:
                total_base_pairs += len(line.strip())
    average_read_size = int(total_base_pairs) / int(total_reads)
    return average_read_size

# Function to get the standard deviation of the read lenghts from a fasta file.
# The function also returns the maximum and minimum read size.
def get_read_lenght_std(input_file):
    current_time = update_time()
    print(current_time + " - Calculating " + str(input_file).split("_", 1)[0]  + " read size standard deviation.")
    read_lenghts = []
    with open(input_file, 'r') as file:
        sequence = ""
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if sequence:
                    read_lenghts.append(len(sequence))
                    sequence = ""
            else:
                sequence += line
        if sequence:
            read_lenghts.append(len(sequence))
    mean_lenght = sum(read_lenghts) / len(read_lenghts)
    variance = sum((length - mean_lenght) ** 2 for length in read_lenghts) / len(read_lenghts)
    std = math.sqrt(variance)
    maximum_read_size = max(read_lenghts)
    minimum_read_size = min(read_lenghts)
    return std, maximum_read_size, minimum_read_size

# Function to write metadata and if there are remaining samples.
def write_metadata(file, file_size, sample_number_of_reads, sample_GC_content, sample_number_of_bases, average_read_size, read_lenght_std, maximum_read_size, minimum_read_size):
    current_time = update_time()
    print(current_time + " - Writing output of metadata to \"metadata.csv\".\n")
    # Writing output metadata in order: file, file size, number of reads in file
    # File %GC content, number of bases in sample, file average read size.
    # read_lenght_std, maximum_read_size, minimum_read_size.
    with open("metadata.csv", "a+") as metadata_output_file:
        metadata_output_file.write(str(file) + ",")
        metadata_output_file.write(str(file_size) + ",")
        metadata_output_file.write(str(sample_number_of_reads) + ",")
        metadata_output_file.write(str(sample_GC_content) + ",")
        metadata_output_file.write(str(sample_number_of_bases) + ",")
        metadata_output_file.write(str(average_read_size) + ",")
        metadata_output_file.write(str(read_lenght_std) + ",")
        metadata_output_file.write(str(maximum_read_size) + ",")
        metadata_output_file.write(str(minimum_read_size) + "\n")
    return

# Function to annotate samples that isn't in the given input directory.
def annotate_sample_that_isnt_in_the_given_directory(sample):
    with open("not_found_samples.csv", "a+") as not_found_samples:
        current_time = update_time()
        print(current_time + " - Writing output of not found samples in the given directory to \"not_found_samples.csv\".")
        not_found_samples.write(str(sample) + " not_found_in_the_given_directory\n")
    return

def main(input_path, reference_file):
    current_time = update_time()
    print(current_time + " - START!\n")
    reference_dict = get_reference_dict(reference_file)
    samples_in_directory_dict = {}
    for file in os.listdir(input_path):
        if file.endswith(".fasta"):
            file_size = get_file_size(file)
            sample_number_of_reads = get_file_number_of_reads(file)
            sample_GC_content, sample_number_of_bases = get_sample_GC_content(file)
            average_read_size = get_average_read_size(file)
            read_lenght_std, maximum_read_size, minimum_read_size = get_read_lenght_std(file)
            write_metadata(file, file_size, sample_number_of_reads, sample_GC_content, sample_number_of_bases, average_read_size, read_lenght_std, maximum_read_size, minimum_read_size)
            # Add file id to "samples_in_directory_dict" to compare if all samples in reference file were detected in the given directory.
            samples_in_directory_dict[str(file).split("_", 1)[0]] = "present_in_the_given_directory"
    samples_not_in_directory_count = 0
    for item in reference_dict:
        if not item in samples_in_directory_dict:
            samples_not_in_directory_count += 1
            annotate_sample_that_isnt_in_the_given_directory(item)
    current_time = update_time()
    print(current_time + " - Number of samples that aren't in the given directory in comparison to reference file: " + str(samples_not_in_directory_count))
    print("\n" + str(current_time) + " - END!")
    return

if __name__ == "__main__":
    main(input_path, reference_file)

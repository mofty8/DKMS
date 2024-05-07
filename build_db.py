import argparse
from functions import *

def main():
    parser = argparse.ArgumentParser(description="Process all FASTQ files in a directory for specified flanking sequences.")
    parser.add_argument("directory", help="Directory containing FASTQ files")
    parser.add_argument("left_flank", help="Left flanking sequence")
    parser.add_argument("right_flank", help="Right flanking sequence")

    args = parser.parse_args()

    # process all FASTQ files in the specified directory with the given flanking sequences
    process_all_fastq_files(args.directory, args.left_flank, args.right_flank)
    print(f"Completed processing all FASTQ files in {args.directory}")

if __name__ == "__main__":
    main()

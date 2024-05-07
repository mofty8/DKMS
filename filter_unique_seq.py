from functions import *
import argparse
import csv



def process_sequence_filter(record, file_path, left_flank, right_flank, max_distance, unique_seqs, insert_flag):
    """Process and filter sequences based on edit distance."""
    sequence_str = str(record.seq)
    if sequence_str.startswith(left_flank) and sequence_str.endswith(right_flank):
        core_sequence = sequence_str[len(left_flank):-len(right_flank)]
        if is_sequence_unique(core_sequence, max_distance):
            unique_seqs[record.description] = core_sequence
            if insert_flag:
                insert_unique_sequence_into_db(core_sequence, file_path, left_flank, right_flank)
        else:
            existing_sequence = sequences_collection.find_one({"sequence_data": core_sequence})
            if existing_sequence:
                data_to_export = {
                    "sequence_id": existing_sequence.get("sequence_id", ""),
                    "sequence_data": existing_sequence["sequence_data"],
                    "occurrences": existing_sequence["occurrences"],
                    "file_occurrences": existing_sequence.get("file_occurrences", {})
                }
                with open('non_unique_sequences.csv', 'a', newline='') as file:
                    writer = csv.writer(file)
                    # Check if file is empty to write headers
                    if file.tell() == 0:
                        writer.writerow(data_to_export.keys())
                    writer.writerow(data_to_export.values())            
                

def process_paired_fastq_filter(forward_path, reverse_path, left_flank, right_flank, max_distance, unique_seqs, insert_flag):
    forward_reads = SeqIO.parse(forward_path, 'fastq')
    reverse_reads = SeqIO.parse(reverse_path, 'fastq')
    for forward_record, reverse_record in zip(forward_reads, reverse_reads):
        process_sequence_filter(forward_record, forward_path, left_flank, right_flank, max_distance, unique_seqs, insert_flag)
        process_sequence_filter(reverse_record, reverse_path, left_flank, right_flank, max_distance, unique_seqs, insert_flag)

def export_unique_sequences(unique_seqs):
    with open("unique_sequences.fasta", "w") as fasta_file:
        for header, seq in unique_seqs.items():
            fasta_file.write(f">{header}\n{seq}\n")


def main():
    parser = argparse.ArgumentParser(description="Filter FASTQ files and optionally insert unique sequences into the database.")
    parser.add_argument("forward_path", help="Path to the forward FASTQ file")
    parser.add_argument("reverse_path", help="Path to the reverse FASTQ file")
    parser.add_argument("left_flank", help="Left flanking sequence")
    parser.add_argument("right_flank", help="Right flanking sequence")
    parser.add_argument("max_distance", type=int, help="Maximum allowed edit distance for sequence uniqueness")
    parser.add_argument("--insert", action='store_true', help="Flag to insert unique sequences into the database")

    args = parser.parse_args()

    unique_seqs = {}
    process_paired_fastq_filter(args.forward_path, args.reverse_path, args.left_flank, args.right_flank, args.max_distance, unique_seqs, args.insert)
    export_unique_sequences(unique_seqs)
    print(f"Processed {len(unique_seqs)} unique sequences.")





if __name__ == "__main__":
    main()


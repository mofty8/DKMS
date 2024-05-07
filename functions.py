
from Bio import SeqIO
import pymongo
from pymongo import UpdateOne
import hashlib
import os
import logging
from mongo_connection import sequences_collection


def generate_sequence_id(sequence):
    """Generate a unique hash for a sequence."""
    sequence_str = str(sequence)
    return hashlib.sha256(sequence_str.encode()).hexdigest()

def process_sequence(sequence, file_path, left_flank, right_flank):
    """Extract and process sequences flanked by specific sequences from FASTQ files."""
    sequence_str = str(sequence)
    if sequence_str.startswith(left_flank) and sequence_str.endswith(right_flank):
        core_sequence = sequence_str[len(left_flank):-len(right_flank)]
        sequence_id = generate_sequence_id(core_sequence)
        
        # Query to find if the sequence already exists in the database
        found = sequences_collection.find_one({"sequence_data": core_sequence})
        
        if found:
            # If the sequence is found, prepare to update the document
            update_operations = [
                UpdateOne(
                    {"_id": found["_id"]},
                    {
                        "$inc": {"occurrences": 1, f"file_occurrences.{file_path}": 1},
                        "$addToSet": {"source_files": file_path}
                    }
                )
            ]
            # Execute the batch of updates
            sequences_collection.bulk_write(update_operations)
        else:
            # If the sequence is not found, create a new document
            sequence_document = {
                "sequence_id": sequence_id,
                "sequence_data": core_sequence,
                "occurrences": 1,
                "file_occurrences": {file_path: 1},
                "flanking_sequences": {"left": left_flank, "right": right_flank},
                "source_files": [file_path]
            }
            sequences_collection.insert_one(sequence_document)
        
        # Log the processing action
        logging.info(f"Processed sequence from {file_path}")


def process_paired_fastq(forward_path, reverse_path, left_flank, right_flank):
    logging.info(f"Processing pairs: {forward_path} and {reverse_path}")
    try:
        forward_reads = SeqIO.parse(forward_path, 'fastq')
        reverse_reads = SeqIO.parse(reverse_path, 'fastq')
        for forward_record, reverse_record in zip(forward_reads, reverse_reads):
            process_sequence(forward_record.seq, forward_path, left_flank, right_flank)
            process_sequence(reverse_record.seq, reverse_path, left_flank, right_flank)
    except Exception as e:
        logging.error(f"Error processing files {forward_path} and {reverse_path}: {e}")

def process_all_fastq_files(directory, left_flank, right_flank):
    files = os.listdir(directory)
    fastq_files = [f for f in files if f.endswith('.fastq')]
    paired_files = {}
    
    # Group files by prefix assuming naming convention like 'sample1_1.fastq', 'sample1_2.fastq'
    for file in fastq_files:
        prefix = '_'.join(file.split('_')[:-1])
        suffix = file.split('_')[-1][0]  # Get '1' or '2' from '1.fastq' or '2.fastq'

        if prefix not in paired_files:
            paired_files[prefix] = [None, None]
        paired_files[prefix][int(suffix) - 1] = os.path.join(directory, file)

    # Process each pair
    for pair in paired_files.values():
        if pair[0] and pair[1]:  # Ensure both forward and reverse files are found
            process_paired_fastq(pair[0], pair[1], left_flank, right_flank)

def process_sequence_filter(record, file_path, left_flank, right_flank, max_distance, unique_seqs):
    """Process and filter sequences based on edit distance."""
    sequence_str = str(record.seq)
    if sequence_str.startswith(left_flank) and sequence_str.endswith(right_flank):
        core_sequence = sequence_str[len(left_flank):-len(right_flank)]
        if is_sequence_unique(core_sequence, max_distance):
            unique_seqs[record.description] = core_sequence
            insert_unique_sequence_into_db(core_sequence, file_path, left_flank, right_flank)
            
def process_paired_fastq_filter(forward_path, reverse_path, left_flank, right_flank, max_distance, unique_seqs):
    forward_reads = SeqIO.parse(forward_path, 'fastq')
    reverse_reads = SeqIO.parse(reverse_path, 'fastq')
    for forward_record, reverse_record in zip(forward_reads, reverse_reads):
        process_sequence_filter(forward_record, forward_path, left_flank, right_flank, max_distance, unique_seqs)
        process_sequence_filter(reverse_record, reverse_path, left_flank, right_flank, max_distance, unique_seqs)

def export_unique_sequences(unique_seqs):
    with open("unique_sequences.fasta", "w") as fasta_file:
        for header, seq in unique_seqs.items():
            fasta_file.write(f">{header}\n{seq}\n")

def filter_fastq_files(forward_path, reverse_path, left_flank, right_flank, max_distance):
    unique_seqs = {}
    process_paired_fastq_filter(forward_path, reverse_path, left_flank, right_flank, max_distance, unique_seqs)
    export_unique_sequences(unique_seqs)
    print(f"Processed {len(unique_seqs)} unique sequences.")

def levenshtein_distance(s1, s2, max_distance):
    if len(s1) > len(s2):
        s1, s2 = s2, s1

    previous_row = range(len(s1) + 1)
    for i, c2 in enumerate(s2):
        current_row = [i + 1]
        # Check if we can already conclude that the distance will be greater than max_distance
        if min(previous_row) > max_distance:
            return max_distance + 1  # Returning max_distance + 1 to indicate the threshold was exceeded

        for j, c1 in enumerate(s1):
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        
        previous_row = current_row
    
    # Check if the final distance is within the allowed threshold
    if previous_row[-1] > max_distance:
        return max_distance + 1  # Returning max_distance + 1 to indicate the threshold was exceeded
    else:
        return previous_row[-1]

def is_sequence_unique(core_sequence, max_distance):
    """Check if the sequence is unique within a specified edit distance threshold."""
    for seq_doc in sequences_collection.find():
        if levenshtein_distance(core_sequence, seq_doc['sequence_data'], max_distance) <= max_distance:
            return False
    return True        

def insert_unique_sequence_into_db(sequence_data, source_file, left_flank, right_flank):
    """Inserts unique sequence into the database if it's not already present."""
    sequence_id = generate_sequence_id(sequence_data)
    
    sequence_document = {
        "sequence_id": sequence_id,
        "sequence_data": sequence_data,
        "occurrences": 1,
        "flanking_sequences": {"left": left_flank, "right": right_flank},
        "file_occurrences": {source_file: 1},
        "source_files": [source_file]
    }
    sequences_collection.insert_one(sequence_document)
    logging.info(f"Unique sequence added to database from {source_file}")
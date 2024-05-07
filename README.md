### MongoDB Installation
Before running this project, please ensure MongoDB is installed and properly set up:

1. Download MongoDB from [MongoDB Community Server](https://www.mongodb.com/try/download/community).
2. Create a `data/db` directory in your home folder to store the database files.
3. Navigate to the `bin` folder inside your MongoDB installation directory.
4. Start MongoDB with the following command:
```bash
./mongod --dbpath /Users/data
```



### Building the Database
To build the database, run the `build_db.py` script with the necessary arguments:
python build_db.py <Folder_path> <left_flank> <right_flank>
**Example:**
```bash
python build_db.py /Users/mofty/FASTQ_files/FastQ_files AGAA CAAT
```



### Filtering Sequences
To filter sequences, execute the `filter_unique_seq.py` script with the required parameters:
```bash
python filter_unique_seq.py <path_file1> <path_file2> <left_flank> <right_flank> <edit_distance_threshold>
```

To insert unique sequences into the database, include the `--insert` flag:
**Example:**
```bash
python filter_unique_seq.py FASTQ_files/ERR103404_1.fastq FASTQ_files/ERR103404_2.fastq AGAA CAAT 10 --insert
```


### Output Files
The script will generate two output files:
- `unique_sequences.fasta` - Contains all unique sequences.
- `non_unique_sequences.csv` - Lists non-unique sequences, detailing the total number of occurrences and the count of occurrences in each individual file.


## Example Outputs
Please refer to the `examples` folder in this repository to review the output files that have been previously generated
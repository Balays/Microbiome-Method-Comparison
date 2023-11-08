#!/bin/bash

# env_parallel --session

# Variables
INPUT_DIR="$1"
FILE_EXTENSION="$2"
OUTPUT_DIR="$3"
NUM_THREADS="$4"  # Change this to your desired number of threads

#PAIRING_METHOD=$(echo "$3" | tr -d '\n')
#echo "$PAIRING_METHOD" | od -c


# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Define a function that will first merge or concatenate pairs based on PAIRING_METHOD
# then execute the three commands for each processed file
process_file() {
	OUTPUT_DIR=sourmash_NCBI_Toti.v1
	sample1="$1"
	base_name=$(basename "$sample1" "$FILE_EXTENSION")
	#sample2="${sample1/_1$FILE_EXTENSION/_2$FILE_EXTENSION}"
	#processed_output="${OUTPUT_DIR}/${base_name}.processed.fq.gz"
	processed_output=$sample1
	#echo 'merging ' $sample1 ' and ' $sample2 'into: ' $processed_output '...'
	
	#cat "$sample1" "$sample2" > "$processed_output"
	# First command
	sourmash sketch dna "$processed_output" -p k=31,k=51,scaled=1000,abund --name "$base_name" -o "${OUTPUT_DIR}/${base_name}.sig.zip"

	# Second command
	sourmash gather "${OUTPUT_DIR}/${base_name}.sig.zip" \
		/mnt/e/data/databases/sourmash/genbank-2022.03-bacteria-k31.zip \
		/mnt/e/data/databases/sourmash/genbank-2022.03-archaea-k31.zip \
		/mnt/e/data/databases/sourmash/genbank-2022.03-viral-k31.zip \
		/mnt/e/data/databases/sourmash/genbank-2022.03-protozoa-k31.zip \
		/mnt/e/data/databases/sourmash/genbank-2022.03-fungi-k31.zip \
		-k 31 -o "${OUTPUT_DIR}/${base_name}.gather.k31.csv"

	# Third command
	sourmash tax metagenome -g "${OUTPUT_DIR}/${base_name}.gather.k31.csv" -t /mnt/e/data/databases/sourmash/genbank-2022.03-*.lineages.csv.gz -o "${OUTPUT_DIR}/${base_name}.gather.k31" -F kreport

}

export -f process_file


# Use parallel to process the files
find "$INPUT_DIR" -name "*$FILE_EXTENSION" | parallel -j $NUM_THREADS process_file {}


# find "$INPUT_DIR" -name "*_1$FILE_EXTENSION" | env_parallel --env PAIRING_METHOD --env FILE_EXTENSION --env OUTPUT_DIR -j $NUM_THREADS process_file {}

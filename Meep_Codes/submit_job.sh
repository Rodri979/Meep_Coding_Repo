#!/bin/bash

# Run this file once pic_run.sh has been edited. It will run all of the nessesary modes

# Count the number of lines in params.txt
num_lines=$(wc -l < params.txt)

# Export the number of lines as an environment variable
export NUM_LINES=$num_lines

echo "NUM_LINES: $NUM_LINES"

# Submit the job with the correct array size
sbatch --export=NUM_LINES pic_run.sh

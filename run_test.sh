#!/bin/bash

rm -rf test_out/

python seq2HLA.py -1 test_data/input/sample_R1.fastq.gz -2 test_data/input/sample_R2.fastq.gz -o test_out/ -p 12

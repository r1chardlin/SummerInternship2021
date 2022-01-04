# runs one dataset instead of looping through a list of them
# this way all of the files can be run at once

import CloneSeqAnalysisLibrary
from CloneSeqAnalysisLibrary import analyze_clone_seq
import os
from os import path
import sys

colony = sys.argv[1]

reference_file = "/local/storage/rhlin_to_copy/CloneSeqExample/Reference.fa"
attempted_file = "/local/storage/rhlin_to_copy/CloneSeqExample/Mutation_Attempts.txt"

sequencePath = "/local/storage/rhlin_to_copy/CloneSeqExample/"
outputPath = "/local/storage/rhlin/ESP_7_output/"
# outputPath = "/local/storage/rhlin/ESP_7_test/"

sequencing_file = sequencePath + "ESP_7_" + colony + "_trimmed.fastq"
fileName = sequencing_file[-20 - len(colony):-6]
samFileName = fileName + ".sam"
outputFile = fileName + "_Summary.txt"

# print("reference_file: " + reference_file)
# print("attempted_file: " + attempted_file)
# print("sequencePath: " + sequencePath)
# print("outputPath: " + outputPath)
# print("sequencing_file: ",  sequencing_file)
# print("samFileName: ", samFileName)
# print("outputFile: ", outputFile)
# quit()

analyze_clone_seq(sequencing_file, attempted_file, reference_file, samFileName, outputFile, outputPath)
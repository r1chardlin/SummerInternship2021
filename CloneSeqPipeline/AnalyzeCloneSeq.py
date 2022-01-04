# import cloneSeqReal
# from cloneSeqReal import analyze_clone_seq

import CloneSeqAnalysisLibrary
from CloneSeqAnalysisLibrary import analyze_clone_seq
import os
from os import path

reference_file = "/local/storage/rhlin_to_copy/CloneSeqExample/Reference.fa"
attempted_file = "/local/storage/rhlin_to_copy/CloneSeqExample/Mutation_Attempts.txt"

sequencePath = "/local/storage/rhlin_to_copy/CloneSeqExample/"
outputPath = "/local/storage/rhlin/ESP_7_output/"
# outputPath = "/local/storage/rhlin/ESP_7_test/"

sequencingFiles = []
samFileNames = []
outputFiles = []
for file in os.listdir(sequencePath):
    # if file.find("downsampled_1_in_100.fastq") != -1:   # for downsampled files
    if file.endswith("trimmed.fastq"):   # for full dataset, endswith returns boolean value so it was just easier to use
        sequencingFiles.append(sequencePath + file)
        fileName = file[:-6]
        samFileNames.append(fileName + ".sam")
        outputFiles.append(fileName + "_Summary" + ".txt")

# print("reference_file: " + reference_file)
# print("attempted_file: " + attempted_file)
# print("sequencePath: " + sequencePath)
# print("outputPath: " + outputPath)
# print("sequencingFiles: ",  sequencingFiles)
# print("samFileNames: ", samFileNames)
# print("outputFiles: ", outputFiles)
# quit()

for i in range(len(sequencingFiles)):
    # if not path.exists(outputFiles[i]):
    analyze_clone_seq(sequencingFiles[i], attempted_file, reference_file, samFileNames[i], outputFiles[i], outputPath)
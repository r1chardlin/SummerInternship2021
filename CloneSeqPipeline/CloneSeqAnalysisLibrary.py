# revised cloneSeqReal to output more information in summary files
# REMEMBER THAT THE READ COUNTS LOG FILE IS 0-INDEXED AND THE MUTATION STRINGS IN THE SUMMARY FILES ARE 1-INDEXED

import os
from os import error, path
import re
from typing import DefaultDict
import numpy as np
import scipy
from scipy.stats import binom
# import json
import tqdm

cigar_codes = {"M": (True, True), "I": (True, False), "D": (False, True), "N": (False, True), "S": (True, False), 
                "H": (False, False), "P": (False, False), "=": (True, True), "X": (True, True)}
quality_thresh = 30

# Calculates the log space sum of two numbers in log space
# @param x The first addend in log space
# @param y The second addend in log space
# @return The sum of x and y in log space
def log_addition(x, y):
    return max(x, y) + np.log1p(np.exp(-abs(x - y)))

# Parses one row of a sam file
# @param row The row of the sam file
# @param header The list of names for each column in the sam file (or keys in the dict)
# @param header The object types of each column in the sam file
# @return The dict containing all of the parsed data
def parse_sam_row(row, header=["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL"], headerTypes = ["str", "int", "str", "int", "int", "str", "str", "int", "int", "str", "str"]):
    samDict = {}
    colList = row.split("\t")
    for i in range(len(header)):
        try:
            value = colList[i]
        except:
            print(colList, i, header[i])
        if headerTypes[i] == "int":
            value = int(value)
        samDict[header[i]] = value
    return samDict

# Parses the whole sam file using parse_sam_row()
# @param filename The name of the sam file to parse
# @param header The list of names for each column in the sam file (technically useless because it is already the default header for parse_sam_row()) 
# @return A DefaultDict mapping each gene to its corresponding list of parsed rows from the sam file
def parse_sam(filename, header=["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL"]):
    # samList = []
    samFile = open(filename)
    refDict = DefaultDict(list)
    print("\nparsing sam file")
    pbar = tqdm.tqdm(total=54554675, desc="parse sam")
    for row in samFile:
        if row[0] == "@":
            continue
        else:
            samRow = parse_sam_row(row, header)
            if samRow["FLAG"] != 4 and samRow["SEQ"] != "*" and samRow["CIGAR"] != "*":
                refDict[samRow["RNAME"]].append(samRow)
            pbar.update()
    return refDict

# Parses the cigar string for a read
# @param sam_dict the dictionary containing the parsed sam row
# @return a list of tuples with each tuple containing the OP the number and OP of the cigar string
def getSplitCigar(sam_dict):
    split_cigar = re.findall("(\d+)([A-Z])", sam_dict["CIGAR"])
    for i in range(len(split_cigar)):
        tup = split_cigar[i]
        temp = (int(tup[0]), tup[1])
        split_cigar[i] = temp
    return split_cigar

# Aligns a read to its reference
# @param sam_dict
# @param geneSeqDict The dict mapping each gene to its corresponding sequence
# @param cigar_codes A dict mapping each OP to a tuple containing two boolean values that indicating whether OP consumes the reference and/or read, respectively
# @return aligned_read The str containing the aligned read
# @return aligned_ref the str containing the aligned reference
def align(sam_dict, geneSeqDict, cigar_codes = cigar_codes):
    readSeq = sam_dict["SEQ"]
    gene = sam_dict["RNAME"]
    refSeq = geneSeqDict[gene]


    # refFasta = open(ref_fasta)
    # isRef = False
    # for line in refFasta:
    #     if line[0] == ">" and line[1:-1] == refID:
    #         isRef = True
    #     elif isRef:
    #         refSeq = line[:-1]
    #         break
    # refFasta.close()
    
    split_cigar = getSplitCigar(sam_dict)
            
    aligned_read = ""
    aligned_ref = ""
    queryPos = 0
    refPos = sam_dict["POS"] - 1
    
    for i in range(len(split_cigar)):
        consumesQuery = cigar_codes[split_cigar[i][1]][0]
        consumesRef = cigar_codes[split_cigar[i][1]][1]
        if consumesQuery and consumesRef:
            aligned_read += readSeq[queryPos:queryPos + split_cigar[i][0]]
            aligned_ref += refSeq[refPos:refPos + split_cigar[i][0]]
            queryPos += split_cigar[i][0]
            refPos += split_cigar[i][0]
        elif consumesQuery and not consumesRef:
            aligned_read += readSeq[queryPos:queryPos + split_cigar[i][0]]
            lenDifference = len(aligned_read) - len(aligned_ref)
            aligned_ref += "-" * lenDifference
            queryPos += split_cigar[i][0]
        elif consumesRef and not consumesQuery:
            aligned_ref += refSeq[refPos:refPos + split_cigar[i][0]]
            lenDifference = len(aligned_ref) - len(aligned_read)
            aligned_read += "-" * lenDifference
            refPos += split_cigar[i][0]
        
    return aligned_read, aligned_ref

# Creates a dict mapping each position its read counts given an aligned read and reference
# @param sam_dict The dict containing the parsed sam row
# @param aligned_read The str containing the aligned read
# @param aligned_ref The str containing the aligned reference
# @param quality_thresh The int indicating the quality threshold to use
# @param cigar_codes A dict mapping each OP to a tuple containing two boolean values that indicating whether OP consumes the reference and/or read, respectively
# @return A dict mapping each pos of the reference to a dictionary mapping each nucleotide to its count (at that position)
def getSingleReadCounts(sam_dict, aligned_read, aligned_ref, quality_thresh = quality_thresh, cigar_codes = cigar_codes):
    readIndex = 0
    refPos = sam_dict["POS"] - 1
    queryPos = 0
    single_read_counts_filtered = {}
    qualScore = sam_dict["QUAL"]
    
    split_cigar = getSplitCigar(sam_dict)

    for i in range(len(split_cigar)):
        consumesQuery = cigar_codes[split_cigar[i][1]][0]
        consumesRef = cigar_codes[split_cigar[i][1]][1]
        cigNum = split_cigar[i][0]
        if consumesQuery and consumesRef:
            try:
                for j in range(cigNum):
                    templateDict = {"A": 0, "T": 0, "C": 0, "G": 0}
                    if ord(qualScore[queryPos]) - 33 >= quality_thresh and aligned_read[readIndex] != "N":
                        nucl = aligned_read[readIndex] # readSeq[queryPos]
                        templateDict[nucl] += 1
                        single_read_counts_filtered[refPos] = templateDict
                    refPos += 1
                    readIndex += 1
                    queryPos += 1
            except:
                print(aligned_read)
                print(aligned_ref)
                print(readIndex)
                print(split_cigar)
                print("consumesQuery: ", consumesQuery)
                print("consumesRef: ", consumesRef)
                print(cigar_codes)
                print(split_cigar[i][1])
                print(cigNum)
                raise
                print(i)
        elif consumesQuery and not consumesRef:
            readIndex += cigNum
            queryPos += cigNum
        elif consumesRef and not consumesQuery:
            # refPos2 = cigNum
            for j in range(cigNum):
    #             templateDict = {"A": 0, "T": 0, "C": 0, "G": 0}
    #             if ord(qualScore[queryPos]) >= quality_thresh:
    #                 single_read_counts_filtered[refPos] = templateDict
    #             else:
    #                 nucl = aligned_ref[readIndex]
    #                 templateDict[nucl] += 1
    #                 single_read_counts_filtered[refPos2] = templateDict
                refPos += 1
                readIndex += 1
    return single_read_counts_filtered

# Creates a dict mapping each position its read counts using the getSingleReadCounts() and align() for a parsed sam file
# @param sam_lines 
# @param geneSeqDict The dict mapping each gene to its corresponding sequence
# @return a DefaultDict mapping each position to its read counts
def get_read_counts(sam_lines, geneSeqDict, quality_thresh=quality_thresh):
    read_counts = DefaultDict(lambda: {"A": 0, "T": 0, "C": 0, "G": 0})
    for sam_dict in sam_lines:
        aligned_read, aligned_ref = align(sam_dict, geneSeqDict)
        singleReadCounts = getSingleReadCounts(sam_dict, aligned_read, aligned_ref, quality_thresh)
        # singleReadCounts = midmod.get_1H(sam_dict, quality_thresh)
        for pos in singleReadCounts:
            for nucl in singleReadCounts[pos]:
                read_counts[pos][nucl] += singleReadCounts[pos][nucl]
    return read_counts
'''
def parseReadCounts(fileName):
    file = open(fileName)
    readCountsDict = {}
    for line in file:
        if line != "":
            lineSplit = line.strip().split("\t")
            print(lineSplit)
            gene = lineSplit[0]
            pos = int(lineSplit[1])
            singleReadCounts = json.loads(lineSplit[2].replace("'", '"'))
            print(singleReadCounts)
            break
            readCountsDict[gene][pos] = singleReadCounts
'''
# Creates A dict mapping each gene to a dictionary that maps position to read counts using the parseReadCounts() or get_read_counts() function
# @param refDict A DefaultDict mapping each gene to its corresponding list of parsed rows from the sam file
# @param geneSeqDict The dict mapping each gene to its corresponding sequence
# @return A dict mapping each gene to a dictionary that maps position to read counts
def generateReadCounts(refDict, geneSeqDict):
    # readCountsFile = outputPath[:-1] + "_OLD_10/" + (os.path.basename(sequencing_file))[:-6] + "_read_counts.txt"
    # if path.exists(readCountsFile):
    #     return parseReadCounts(readCountsFile)
    # else:
    readCountsDict = {}
    print("\ngenerating read counts")
    for gene in tqdm.tqdm(refDict, desc="read counts"):
        readCounts = get_read_counts(refDict[gene], geneSeqDict)
        readCountsDict[gene] = readCounts
    return readCountsDict
'''
def getGeneSeq(gene, ref_fasta):
    geneSeq = ""
    refFasta = open(ref_fasta)
    isRef = False
    for line in refFasta:
        if line[0] == ">" and line[1:-1] == gene:
            isRef = True
        elif isRef:
            if line[-1] == "\n":
                geneSeq = line[:-1]
            else:
                geneSeq = line
            break
    refFasta.close()
    return geneSeq
'''
'''
def getTotalNucls(readCount):
    totalNucls = 0
    for nucl in readCount:
        totalNucls += readCount[nucl]
    return totalNucls
'''
# Calculates the error rate for sequencing reads
# @param readCounts The dict mapping each gene to a dictionary that maps position to read counts
# @param geneSeqDict The dict mapping each gene to its corresponding sequence
# @return The calculated error rate
def getErrorRate(readCounts, geneSeqDict):
    totalErrors = 0
    totalReads = 0
    # geneSeq = {gene: seq for gene in readCounts}
    print("\ngetting error rate")
    for gene in tqdm.tqdm(readCounts, desc="error rate"):
        # geneSeq = getGeneSeq(gene, refFasta)
        # readCounts = get_read_counts(refDict[gene], refFasta)
        geneSeq = geneSeqDict[gene]
        for pos in readCounts[gene]:
            totalNucls = readCounts[gene][pos]["A"] + readCounts[gene][pos]["T"] + readCounts[gene][pos]["C"] + readCounts[gene][pos]["G"]
            for nucl in readCounts[gene][pos]:
                nuclCount = readCounts[gene][pos][nucl]
                # nuclPercentage = nuclCount / totalNucls
                if nucl != geneSeq[pos]: # and nuclPercentage < 0.01:
                    totalErrors += nuclCount
            totalReads += totalNucls
    return totalErrors / totalReads

# Determines whether attempted mutations were generated and unwanted mutations using a cdf test
# @param gene The gene/ORF to analyze
# @param readCounts The dict that maps position to read counts for the gene
# @param geneSeqDict The dict that maps each gene to its corresponding sequence
# @param mutFile The txt file containing all the attempted mutations
# @param errorRate The calculated error rate for sequencing the reads
# @param outputPath The path to output the cdf log file
# @param colony The colony that is being analyzed
# @return mut_status A dict that maps each attempted mutation to a boolean value indicating whether it was successfully generated
# @return unwantedMuts A list of detected unwanted mutations
# @return clean A boolean indicating if there are any detected unwanted mutations
# @return score The score for each attempted mutation
def getMutants(gene, readCounts, geneSeqDict, mutFile, errorRate, outputPath, colony):
    # geneSeq = ""
    # refFasta = open(ref_fasta)
    # isRef = False
    # for line in refFasta:
    #     if line[0] == ">" and line[1:-1] == gene:
    #         isRef = True
    #     elif isRef:
    #         if line[-1] == "\n":
    #             geneSeq = line[:-1]
    #         else:
    #             geneSeq = line
    #         break
    # refFasta.close()
    
    geneSeq = geneSeqDict[gene]

    # attemptedMutsFile = open(mutFile)
    # attemptedMutsGene = []
    # for line in attemptedMutsFile:
    #     if gene in line:
    #         attemptedMutsGene.append(line.strip())
    # attemptedMutsFile.close()

    attemptedMutsFile = open(mutFile)
    attemptedMuts = []
    for line in attemptedMutsFile:
        attemptedMutStr = line.strip().split("\t")[1]
        # could use ORF ID instead (line.strip().split("\t")[0]), doing it this way just to be safe (gets it from the mutation ID)
        if gene == attemptedMutStr.split("_")[0]:
            attemptedMuts.append(attemptedMutStr)
    attemptedMutsFile.close()
    
    mutDict = {}
    for mutation in attemptedMuts:
        temp = mutation.split("_")
        pos = int(temp[1][1:-1])
        mutNucl = temp[1][-1]
        if pos in mutDict:
            if type(mutDict[pos]) is not list:
                mutDict[pos] = [mutDict[pos]]
            mutDict[pos].append(mutNucl)
        else:
            mutDict[pos] = mutNucl
        # print(temp)
        # print(pos)
        # print(mutNucl)
    # print(gene)
    # print(mutDict)
    
    expectedMutPercentage = (1 / len(attemptedMuts)) - errorRate
    
    foundAttemptedMuts = []
    unwantedMuts = []
    unwantedMutPercentages = []
    attemptedCdfs = {}
    stopPos = len(geneSeq) - 3
    for pos in readCounts:
        totalNucls = totalNucls = readCounts[pos]["A"] + readCounts[pos]["T"] + readCounts[pos]["C"] + readCounts[pos]["G"]
        if totalNucls == 0:
            totalNucls += 999999 * len(attemptedMuts)
        if pos + 1 in mutDict:
            mut = mutDict[pos + 1]
        else:
            mut = "0"   # Nonetype is not iterable
        for nucl in readCounts[pos]:
            # if gene == "1320" and pos == 1107:
            #     print(gene, nucl, readCounts[pos])
            if nucl == geneSeq[pos]:
                continue
            elif nucl in mut:
                nuclCount = readCounts[pos][nucl]
                nuclPercentage = nuclCount / totalNucls
                flag = False
                if totalNucls <= 20:
                    flag = True
                mutStr = gene + "_" + geneSeq[pos] + str(pos + 1) + nucl
                cdf_mut = binom.logcdf(nuclCount, totalNucls, expectedMutPercentage)
                cdf_no_mut = log_addition(binom.logsf(nuclCount, totalNucls, errorRate), binom.logpmf(nuclCount, totalNucls, errorRate))
                cdf_scaled = binom.logcdf((nuclCount / totalNucls) * 100, 100, expectedMutPercentage)
                attemptedCdfs[mutStr] = cdf_mut
                if ((cdf_mut > cdf_no_mut) != (cdf_scaled >= np.log(0.05))) and not flag:
                    cdfFile = open(outputPath + "cdf_log.txt", "a")
                    cdfFile.write(gene + "\t" \
                                + colony + "\t" \
                                + mutStr + "\t" \
                                + str(expectedMutPercentage) + "\t" \
                                + str(cdf_mut) + "\t" \
                                + str(cdf_no_mut) + "\t" \
                                + str(cdf_scaled) + "\n")
                if cdf_mut > cdf_no_mut and not flag:    # if cdf >= 0.05:
                    foundAttemptedMuts.append(mutStr)
                    # print(mutStr)
            elif readCounts[pos][nucl] > 0:
                nuclCount = readCounts[pos][nucl]
                nuclPercentage = nuclCount / totalNucls
                thresholdPerc = 0.15 + 0.0625 * (40 / totalNucls)
                lowReadHighPerc = False
                if (totalNucls <= 90 and nuclPercentage < thresholdPerc) or (pos >= stopPos) or (totalNucls <= 20):   # or (totalNucls == 1 and nuclPercentage == 1.0)
                    lowReadHighPerc = True
                cdf_mut = binom.logcdf(nuclCount, totalNucls, expectedMutPercentage)
                cdf_no_mut = log_addition(binom.logsf(nuclCount, totalNucls, errorRate), binom.logpmf(nuclCount, totalNucls, errorRate))
                # cdf2 = (1 - binom.cdf(nuclCount, totalNucls, 0.05)) + binom.pmf(nuclCount, totalNucls, errorRate)
                if cdf_mut > cdf_no_mut and not lowReadHighPerc: # !=  (cdf2 <= 0.05)):
                    # print(gene, pos, nucl, nuclPercentage)
                    mutStr = gene + "_" + geneSeq[pos] + str(pos + 1) + nucl
                    unwantedMuts.append(mutStr)
                    unwantedMutPercentages.append(nuclPercentage)
                    # print(mutStr, readCounts[pos], cdf_no_mut, cdf_mut, totalNucls, nuclPercentage)

    # if len(unwantedMutPercentages) > 0:
    #     print(unwantedMutPercentages)
    #     print(sum(unwantedMutPercentages))
    #     quit()

    # print(gene, unwantedMutPercentages, sum(unwantedMutPercentages))

    # print("-----------------------------")
    # print(gene)
    # print()
    # print(attemptedMuts)
    # print(foundAttemptedMuts)
    # print(attemptedCdfs)
    # print()
    # print(unwantedMuts)
    # print(unwantedMutPercentages)
    # print(sum(unwantedMutPercentages))

    mut_status = {mut: False for mut in attemptedMuts}
    scores = {}
    for mut in mut_status:
        try:
            scores[mut] = attemptedCdfs[mut] - sum(unwantedMutPercentages)
        except:
            scores[mut] = -5 - sum(unwantedMutPercentages)
        if mut in foundAttemptedMuts:
            mut_status[mut] = True

    # print()
    # print(scores)
    # print("-----------------------------")
    
    clean = True
    if len(unwantedMuts) > 0:
        clean = False

    return mut_status, unwantedMuts, clean, scores   # sum(unwantedMutPercentages)   # return mutants, sum(unwantedMutPercentages)
'''
def getMutResults(mutants, attemptedMuts):     # def getMutResults(gene, mutants, mutFile):
    # attemptedMutsFile = open(mutFile)
    # attemptedMutsGene = []
    # for line in attemptedMutsFile:
    #     attemptedMutStr = line.strip().split("\t")[1]
    #     if gene == attemptedMutStr.split("_")[0]:    # could use ORF ID instead (line.strip().split("\t")[0]), doing it this way just to be safe (gets it from the mutation ID)
    #         attemptedMutsGene.append(attemptedMutStr)

    # attemptedMutsFile.close()
    
    mut_status = {mut: False for mut in attemptedMuts}
    clean = True
    unwantedMuts = []
    for mut in mutants:
        # print(mut)
        if mut in attemptedMuts:
            mut_status[mut] = True
        else:
            clean = False
            unwantedMuts.append(mut)
    
    # print("attempted muts: ", attemptedMutsGene, "\n")
    # print("mut status: ", mut_status, "\n")
    # print("unwanted muts: ", unwantedMuts)

    # otherMutStr = ""
    # for i in range(len(unwantedMuts)):
    #     if i < len(unwantedMuts) - 1:
    #         otherMutStr += unwantedMuts[i] + ";"
    #     else:
    #         otherMutStr += unwantedMuts[i]
        
    return mut_status, clean, unwantedMuts
'''
# Uses all the functions above to generate a summary file for the experiment
# @param sequencing_file The fastq file to use for generating the sam file
# @param attempted_file The txt file that contains all the attempted mutations
# @param reference_file The fasta file that contains all the reference orf sequences
# @param samFileName The name of the sam file to create
# @param output_file The name of the txt to write to
# @param output_path The path to write the output file to
# @return None
def analyze_clone_seq(sequencing_file, attempted_file, reference_file, samFileName, output_file, output_path):
    seqFileBase = os.path.basename(sequencing_file)

    colony = ""
    if seqFileBase[:5] == "ESP_7":
        colony = seqFileBase[6]
    
    if not path.exists(samFileName):
        os.system("bwa index -p refIndex -a is " + reference_file)
        os.system("bwa mem -t 8 refIndex " + sequencing_file + " > " + samFileName)

    refDict = parse_sam(samFileName)

    # cloneSeqParsed = parse_sam(samFileName)
    
    '''
    refDict = {}
    i = 0
    print("filtering reads\n")
    
    pbar = tqdm.tqdm(total=len(cloneSeqParsed))
    while i < len(cloneSeqParsed):
        if cloneSeqParsed[i]["FLAG"] == 4:  # no need to specifically delete the reads that do not align to 
            cloneSeqParsed.pop(i)           # the reference sequence because if a read did not align  
            i -= 1                          # it would have a FLAG of 4
        elif cloneSeqParsed[i]["RNAME"] not in refDict.keys():  
            refDict[cloneSeqParsed[i]["RNAME"]] = [cloneSeqParsed[i]]
        else:
            refDict[cloneSeqParsed[i]["RNAME"]].append(cloneSeqParsed[i])
        i += 1
        pbar.update()
    '''
    # refDict = DefaultDict(list)
    # for row in tqdm.tqdm(cloneSeqParsed):
    #     if row["FLAG"] == 4:
    #         continue
    #     else:
    #         refDict[row["RNAME"]].append(row)

    cloneSeqSummary = open(output_path + output_file, "w")
    logFile = open(output_path + (os.path.basename(sequencing_file))[:-6] + "_read_counts.txt", "w")

    refFasta = open(reference_file)
    geneSeqDict = {}
    gene = ""
    for line in refFasta:
        if line[0] == ">":
            gene = line[1:-1]
        else:
            geneSeqDict[gene] = line.strip()
    refFasta.close()
    
    readCounts = generateReadCounts(refDict, geneSeqDict)
    errorRate = getErrorRate(readCounts, geneSeqDict)
    
    print("\nanalyzing genes")
    for gene in tqdm.tqdm(geneSeqDict, desc="analyze genes"):
        try:
            readCountsGene = readCounts[gene]
        except:
            readCountsGene = {i: {"A": 0, "T": 0, "C": 0, "G": 0} for i in range(len(geneSeqDict[gene]))}
        for pos in readCountsGene:
            logFile.write(gene + "\t" + str(pos) + "\t" + str(readCountsGene[pos]) + "\n")
        mutStatus, otherMuts, clean, scores = getMutants(gene, readCountsGene, geneSeqDict, attempted_file, errorRate, output_path, colony)
        # mutStatus, clean, otherMuts = getMutResults(mutants, attemptedMuts)
        for mut in mutStatus:
            modMutStr = mut + "_" + colony
            score = scores[mut]
            # score = -deduction
            orf = mut.split("_")[0]
            successful = False
            # passedBenchmark = True
            if mutStatus[mut]:
                if clean:
                    successful = True
                # score += 1
            # if score <= 0:
            #     passedBenchmark = False
            cloneSeqSummary.write(orf + "\t" + colony + "\t" + modMutStr + "\t" \
                                + str(mutStatus[mut]) + "\t" + str(otherMuts) + "\t" \
                                + str(clean) + "\t" + str(successful) + "\t" + str(score) + "\n") # + str(passedBenchmark) + "\n")
    cloneSeqSummary.close()
    logFile.close()
    
    return


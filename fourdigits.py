
import hashlib
import os
import sys


import numpy as np

from confidence import calculate_confidence
import config as cfg

"""
Determination of the 4-digit HLA type based on the 2-digit determination by seq2HLA.py
"""

top_alleles = {}

def determine_four_digits_main(infile, bgVec, hla_class, hla_classification, ambiguityfile):
    if hla_classification == "nonclassical":
        result = determine_four_digits_minor_alleles(infile, bgVec, ambiguityfile)
    else:
        dbmhc_file = ""
        dbmhc = []
        if hla_class == 1:
            dbmhc_file = cfg.dbmhc_I
        else:
            dbmhc_file = cfg.dbmhc_II
        dbmhc = create_db_mhc_dict(dbmhc, dbmhc_file)
        result = determine_four_digits(dbmhc, infile, bgVec, ambiguityfile)
    return result
	

def create_db_mhc_dict(dbmhc, dbmhc_file):
    population_dict = {}
    header = True
    with open(dbmhc_file) as inf:
        for line in inf:
            eles = line.split("\t")
            if header:
                #header
                for i in range(1, len(eles)):
                    population_dict[i] = eles[i]
                header = False
            else:
                if ":" in eles[0]:
                    l0_split = eles[0].split(":")
                    fourdigit = "{}:{}".format(l0_split[0], l0_split[1])
                    dbmhc.append(fourdigit)
				
    return dbmhc


def mean_pan_pop_freq(freq_dict):
    pan_list = []
    for entry in freq_dict:
        pan_list.append(float(freq_dict[entry]))
    return np.mean(pan_list)


def determine_four_digits(dbmhc, infile, bgVec, ambiguityfile):    
    result_dict = {}
    result_dict2 = {}
    readcount = []
    ambiguity_dict = {}
    result = []
    allele = ""
    with open(infile) as inf:
        for line in inf:
            eles = line.split(":")
            for i in range(0, len(eles) - 1):
                allele += "{}:".format(eles[i])
            result_dict[allele[0:-1]] = int(eles[len(eles) - 1])
            allele = ""
            readcount.append(int(eles[len(eles) - 1]))
	
    writehandle = open("{}.solutions".format(infile), "w")
    writehandle.write("#Full allele\tp-value\n")
	
    readcount = [x for x in readcount if x != 0]
	
    top = readcount[0]
    readcount_copy = readcount[1:]
    readcount_copy_str = str(readcount_copy)[1:-1].replace(" ","")

    confidence = calculate_confidence(top, readcount)

    repeat = True

    if confidence:
        if 0.0 < confidence < 0.001:
            for item in result_dict:
                if result_dict[item] == top:
                    writehandle.write("{}\t{}\n".format(item, 1 - calculate_outlier(result_dict[item], bgVec)))
                    repeat = False
                    mostProbableAllele = str(item.split(":")[0]) + ":" + str(item.split(":")[1])
                    numbersolutions = 1
			
    percentile = 95
    while repeat:
        if percentile == 0:
            return "no,0"
        percentile -= 5
        perc = int(np.percentile(readcount, percentile))
        if perc == 0:
            perc += 0.1
        count = 0
        for item in result_dict:
            if result_dict[item] >= perc:
                count += 1
		
        #there is only 1 solution:
        if count == 1:
            count = 0
            for item in result_dict:
                if result_dict[item] >= perc:
                    writehandle.write("{}\t{}\n".format(item, 1 - calculate_outlier(result_dict[item], bgVec)))
                    repeat = False
                    item_split = item.split(":")
                    mostProbableAllele = "{}:{}".format(item_split[0], item_split[1])
                    numbersolutions = 1
                #there are more possible solutions
        else:
            count = 0
            #check if allele (4-digits) occurs in dbMHC table. If not, is it a very rare allele in all populations
            for item in result_dict:
                item_split = item.split(":")
                item4digit = "{}:{}".format(item_split[0], item_split[1])
                if result_dict[item] >= perc and item4digit in dbmhc:
                    result.append(item4digit)
                    result_dict2[item] = result_dict[item]
                #only 1 solution (at 4-digit level) left
                if len(set(result)) == 1:
                    #for item in result_dict:
                    for item in result_dict2:
                        item_split = item.split(":")
                        writehandle.write("{}\t{}\n".format(item, 1 - calculate_outlier(result_dict2[item], bgVec)))
                        repeat = False
                        mostProbableAllele = "{}:{}".format(item_split[0], item_split[1])
                        numbersolutions = 1
                    #more than one 2 solution possible (at 4-digit level).
                elif len(set(result)) > 1:
                    numbersolutions = len(set(result))
                    max = 0.0
				
                    for item in result_dict2:
                        item_split = item.split(":")
                        item4digit = "{}:{}".format(item_split[0], item_split[1])
                        if not item4digit in ambiguity_dict:
                            ambiguity_dict[item4digit] = result_dict2[item]
                        elif ambiguity_dict[item4digit] < result_dict2[item]:
                            ambiguity_dict[item4digit] = result_dict2[item]	
                        if result_dict2[item] > max:
                            max = result_dict2[item]
                            mostProbableAllele = item4digit
                        ambiguity_handle = open(ambiguityfile, "a")
                        ambiguity_handle.write("#################################################################################################\n")
                        ambiguity_handle.write("#Ambiguity:\n")
                        ambiguity_handle.write("#Based on the RNA-Seq reads and the dbMHC table, the following 4-digits alleles are possible:\n")
                        for ambi4DigitAllele in ambiguity_dict:
                            ambiguity_handle.write("{}\t{}\n".format(ambi4DigitAllele, 1 - calculate_outlier(ambiguity_dict[ambi4DigitAllele], bgVec)))
				
                        ambiguity_handle.write("#However, by taking into account the read data, the most probable 4-digit allele is:\n")
                        ambiguity_handle.write("{}\n\n".format(mostProbableAllele))
                        ambiguity_handle.close()
                        
                        for item in result_dict2:
                            item_split = item.split(":")
                            item4digit = "{}:{}".format(item_split[0], item_split[1])
                            if result_dict2[item] == max:
                                writehandle.write("{}\t{}\n".format(item, 1 - calculate_outlier(result_dict2[item], bgVec)))
                                repeat = False
    writehandle.close()
    return "{},{}".format(mostProbableAllele, numbersolutions)
	
def calculate_outlier(top, bgVec):
    """Call R-script "commmand_fourdigit.R" to calculate the confidence of the top-scoring allele (four-digits)."""

    # Convert bgVec into md5
    md5sum = hashlib.md5(bgVec.encode("utf-8")).hexdigest()
    if not md5sum in top_alleles:
        top_alleles[md5sum] = {}
    if top in top_alleles[md5sum]:
        return top_alleles[md5sum][top]
    else:
        allele_list = [float(x) for x in bgVec.rstrip(",").split(",")]
        confidence = calculate_confidence(float(top), allele_list)
        top_alleles[md5sum][top] = confidence
        return confidence

#method for determining the 4 digit HLA type for minor HLA alleles, for which no entry in dbMHC is available
def determine_four_digits_minor_alleles(infile, bgVec, ambiguityfile):
	
    result_dict = {}
    result_dict2 = {}
    readcount = []
    ambiguity_dict = {}
    result = []
    allele = ""
    with open(infile) as inf:
        for line in inf:
            l = line.split(":")
            for i in range(0, len(l) - 1):
                allele += "{}:".format(l[i])
            result_dict[allele[0:-1]] = int(l[len(l)-1])
            allele = ""
            readcount.append(int(l[len(l)-1]))
	
    writehandle = open("{}.solutions".format(infile), "w")
    writehandle.write("#Full allele\tp-value\n")
	
    readcount = [x for x in readcount if x != 0]
	
    top = readcount[0]
    readcount_copy = readcount[1:]
    readcount_copy_str = str(readcount_copy)[1:-1].replace(" ","")
	
    confidence = calculate_confidence(top, readcount)
    repeat = True

    if confidence:
        if 0.0 < confidence < 0.001:
            for item in result_dict:
                if result_dict[item] == top:
                    writehandle.write("{}\t{}\n".format(item, 1 - calculate_outlier(result_dict[item], bgVec)))
                    repeat = False
                    item_split = item.split(":")
                    mostProbableAllele = "{}:{}".format(item_split[0], item_split[1])
                    numbersolutions = 1
			
    percentile = 95
    while repeat:
        if percentile == 0:
            return "no,0"
        percentile -= 5
        perc = int(np.percentile(readcount, percentile))
        if perc == 0:
            perc += 0.1
        count = 0
        for item in result_dict:
            if result_dict[item] >= perc:
                count += 1
		
        #there is only 1 solution:
        if count == 1:
            count = 0
            for item in result_dict:
                if result_dict[item] >= perc:
                    writehandle.write("{}\t{}\n".format(item, 1 - calculate_outlier(result_dict[item], bgVec)))
                    repeat = False
                    item_split = item.split(":")
                    mostProbableAllele = "{}:{}".format(item_split[0], item_split[1])
                    numbersolutions = 1
        #there are more possible solutions
        else:
            count = 0
            #check if allele (4-digits) is in the defined percentile range
            for item in result_dict:
                item_split = item.split(":")
                item4digit = "{}:{}".format(item_split[0], item_split[1])
                if result_dict[item] >= perc:
                    result.append(item4digit)
                    result_dict2[item] = result_dict[item]
            #only 1 solution (at 4-digit level) left
            if len(set(result)) == 1:
                #for item in result_dict:
                for item in result_dict2:
                    writehandle.write("{}\t{}\n".format(item, 1 - calculate_outlier(result_dict2[item], bgVec)))
                    repeat = False
                    item_split = item.split(":")
                    mostProbableAllele = "{}:{}".format(item_split[0], item_split[1])
                    numbersolutions = 1
                #more than one 2 solution possible (at 4-digit level).
            elif len(set(result)) > 1:
                numbersolutions = len(set(result))
                max = 0.0
				
                for item in result_dict2:
                    item_split = item.split(":")
                    item4digit = "{}:{}".format(item_split[0], item_split[1])
                    if not item4digit in ambiguity_dict:
                        ambiguity_dict[item4digit] = result_dict2[item]
                    elif ambiguity_dict[item4digit] < result_dict2[item]:
                        ambiguity_dict[item4digit] = result_dict2[item]	
                    if result_dict2[item] > max:
                        max = result_dict2[item]
                        #if meanPanPopFreq(dbmhc_prob[item4digit]) > max:
                        #	max = meanPanPopFreq(dbmhc_prob[item4digit])
                        mostProbableAllele = item4digit
                        #if np.mean(dbmhc_prob[item4digit]) > max:
                        #	max = np.mean(dbmhc_prob[item4digit])
                        #	mostProbableAllele=item4digit
                        ambiguity_handle = open(ambiguityfile, "a")
                        ambiguity_handle.write("#################################################################################################\n")
                        ambiguity_handle.write("#Ambiguity:\n")
                        ambiguity_handle.write("#Based on the RNA-Seq reads and the dbMHC table, the following 4-digits alleles are possible:\n")
                        for ambi4DigitAllele in ambiguity_dict:
                            ambiguity_handle.write("{}\t{}\n".format(ambi4DigitAllele, 1 - calculate_outlier(ambiguity_dict[ambi4DigitAllele], bgVec)))
                            
                        ambiguity_handle.write("#However, by taking into account the read data, the most probable 4-digit allele is:\n")
                        ambiguity_handle.write("{}\n\n".format(mostProbableAllele))
                        ambiguity_handle.close()
				
                        for item in result_dict2:
                            item_split = item.split(":")
                            item4digit = "{}:{}".format(item_split[0], item_split[1])
                            #if np.mean(dbmhc_prob[item4digit]) == max:
                            #if meanPanPopFreq(dbmhc_prob[item4digit]) == max:
                            if result_dict2[item] == max:
                                writehandle.write("{}\t{}\n".format(item, 1 - calculate_outlier(result_dict2[item], bgVec)))
                                repeat = False
    writehandle.close()			
    return "{},{}".format(mostProbableAllele, numbersolutions)

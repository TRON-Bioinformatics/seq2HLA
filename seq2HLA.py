#!/usr/bin/env python

"""
Seq2HLA is an in-silico method, written in python and R, which takes standard RNA-Seq sequence reads in fastq format as input, uses a bowtie index comprising all HLA alleles and outputs the most likely HLA class I and class II genotypes (in 4 digit resolution), a p-value for each call, and the expression of each class.
"""

__author__ = "Sebastian Boegel"
__copyright__ = "Copyright (c) 2023 Bioinformatics Group at TRON - Translational Oncology at the Medical Center of the Johannes Gutenberg-University Mainz gGmbH"
__license__ = "MIT"
__maintainer__ = "Patrick Sorn"
__mail__ = "patrick.sorn@tron-mainz.de"


from argparse import ArgumentParser
import gzip
import linecache
import logging
from operator import itemgetter
import os
import subprocess
import sys

# External modules
from Bio import SeqIO
import numpy as np

# Internal modules
import config as cfg

from determine_hla import calculate_confidence
from determine_hla import determine_four_digits_main
from version import version


def is_gzipped(infile):
    gzipped = False
    f = gzip.open(infile, "r")
    try:
        first_line = f.readline()
        gzipped = True
    except:
        gzipped = False
    return gzipped


def get_read_len(fq_file, gzipped):
    cmd = None

    if gzipped:
        cmd = "zcat {} | sed '2q;d' | wc -L".format(fq_file)
    else:
        cmd = "sed '2q;d' {} | wc -L".format(fq_file)
            
    process_wc = subprocess.Popen(["bash", "-c", cmd], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    read_len = process_wc.communicate()[0]
    return read_len



def get_P(locus, finaloutput):
    """
    Return the p value of a prediction, which is stored in the intermediate textfile.
    """
    locus_dict = {
        "A": 2, "DQA1": 2, "E": 2,
        "B": 3, "DQB1": 3, "F": 3,
        "C": 4, "DRB1": 4, "G": 4,
        "H": 5, "DRA": 5,
        "J": 6, "DPA1": 6,
        "K": 7, "DPB1": 7,
        "L": 8,
        "P": 9,
        "V": 10
    }
    return linecache.getline(finaloutput, locus_dict[locus]).split("\t")[4][0:-1]



def get_best_P(infile):
    """
    In case of ambiguous typings, the allele(s) with the best p-value
    (which is actually the smallest one) is reported and thus the min(p) is returned.
    """
    p = []
    with open(infile) as inf:
        for line in inf:
            if not line[0] == "#":
                if line.split("\t")[1][0:-1]=="NA":
                    return "NA"
                p.append(float(line.split("\t")[1][0:-1]))
    try:
        return min(p)
    except ValueError:
        return "NA"



def get_best_allele_per_group(readspergroup, readcount):
    """
    For each allele, to which at least 1 read map, find the allele which has the most reads within a group (2-digit-level, e.g. A*02") and save
    #i) the key of this allele as ambassador for the group => maxallelepergroup holds for each group the top-scoring allele
    #ii) the number of reads mapping to the top-scoring allele => readspergroup
    """

    top_counts = {}
    top_alleles = {}
    for group in readspergroup:
        top_counts[group] = readspergroup[group]
    for key in readcount:
        if readcount[key] > 0:
            group = key.split(":")[0]
            if top_counts[group] <= readcount[key]:
                top_counts[group] = readcount[key]
                top_alleles[group] = key
    return (top_counts, top_alleles)



def parse_mapping(hla_dict, sam, readcount):
    """Open sam-file and count mappings for each allele."""
    with open(sam) as samhandle:
        for line in samhandle:
            if line[0] != '@':
                hla = line.split('\t')[2]
                readcount[hla_dict[hla]] += 1
    return readcount



class Pipeline(object):
    def __init__(self, working_dir, fq1, fq2, trim3, threads):
        self.working_dir = os.path.abspath(working_dir.rstrip("/"))
        if not os.path.exists(self.working_dir):
            os.makedirs(self.working_dir)
        self.fq1 = fq1
        self.fq2 = fq2
        self.trim3 = trim3
        self.threads = threads
        self.gzipped = is_gzipped(self.fq1)

        
        self.logger = logging.getLogger("main")
        self.logger.handlers = []
        self.logger.setLevel(logging.DEBUG)

        formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")

        streamHandler = logging.StreamHandler(sys.stdout)
        streamHandler.setLevel(logging.DEBUG)
        streamHandler.setFormatter(formatter)

        fileHandler = logging.FileHandler(os.path.join(self.working_dir, "run.log"), "w+")
        fileHandler.setLevel(logging.DEBUG)
        fileHandler.setFormatter(formatter)

        #self.logger = logging.getLogger("main")
        #self.logger.addHandler(streamHandler)
        
        self.logger.addHandler(streamHandler)
        self.logger.addHandler(fileHandler)



    def run(self):
        """This function starts predicting the HLA type for the input reads."""
        script_call = "python {} -1 {} -2 {} -3 {} -o {} -p {}".format(
            os.path.realpath(__file__),
            self.fq1,
            self.fq2,
            self.trim3,
            self.working_dir,
            self.threads
        )

        self.logger.info("Starting seq2HLA pipeline version {}!".format(version))
        self.logger.info("CMD: {}".format(script_call))
        self.logger.info("Estimating mismatch ratio from read length.")
        
        read_len = get_read_len(self.fq1, self.gzipped)
        mismatch_ratio = round(int(read_len) / 50.0)

        mapopt = "-p {} -a -v {}".format(self.threads, mismatch_ratio)

        self.logger.info("The read length of your input fastq was determined to be {}, so {} mismatches will be allowed and {} threads will be used by bowtie.".format(read_len, mismatch_ratio, self.threads))


	# call HLA typing for Class I, classical (A,B.C)
        self.call_HLA(
            cfg.bowtiebuild_class_I, 
            cfg.fasta_class_I, 
            mapopt, 
            cfg.class_I_list, 
            1, 
            "classical", 
            cfg.len_dict_I
        )

	# call HLA typing for Class I, non-classical (E,G....)
        self.call_HLA(
            cfg.bowtiebuild_class_I_nonclass, 
            cfg.fasta_class_I_nonclass, 
            mapopt, 
            cfg.class_I_nonclass_list, 
            1, 
            "nonclassical", 
            cfg.len_dict_I_nonclass
        )

        # call HLA typing for Class II
        self.call_HLA(
            cfg.bowtiebuild_class_II,
            cfg.fasta_class_II,
            mapopt,
            cfg.class_II_list,
            2,
            "classical",
            cfg.len_dict_II
        )


    def call_HLA(self, bowtiebuild, hla1fasta, mapopt, locus_list, hla_class, hla_classification, length_dict):
        """This method calls HLA types for class I."""
        twodigits1 = {}
        fourdigits1 = {}
        fourDigit_solutions1 = {}
        twodigits2 = {}
        fourdigits2 = {}
        fourDigit_solutions2 = {}
        fourDigit_solutions3 = {}
        twodigits3 = {}
        finalAlleles = {}

        if hla_class == 1:
            self.logger.info("----------HLA class I------------")
            if hla_classification == "classical":
                self.logger.info(">---classical HLA alleles---")
            else:
                self.logger.info(">---nonclassical HLA alleles---")
        else:
            self.logger.info("----------HLA class II------------")
        prefix = "hla_{}_{}".format(hla_class, hla_classification)

        
        fq1_2 = os.path.join(self.working_dir, "{}_iteration_2_1.fq".format(prefix))
        fq2_2 = os.path.join(self.working_dir, "{}_iteration_2_2.fq".format(prefix))
        sam1 = os.path.join(self.working_dir, "{}_iteration_1.sam".format(prefix))
        sam2 = os.path.join(self.working_dir, "{}_iteration_2.sam".format(prefix))
        sam3 = os.path.join(self.working_dir, "{}_iteration_3.sam".format(prefix))
        aligned = os.path.join(self.working_dir, "{}.aligned".format(prefix))
        aligned_1 = os.path.join(self.working_dir, "{}_1.aligned".format(prefix))
        aligned_2 = os.path.join(self.working_dir, "{}_2.aligned".format(prefix))
        bowtielog = os.path.join(self.working_dir, "{}.bowtielog".format(prefix))
        output1 = os.path.join(self.working_dir, "{}.digitalhaplotype1".format(prefix))
        output2 = os.path.join(self.working_dir, "{}.digitalhaplotype2".format(prefix))
        output3 = os.path.join(self.working_dir, "{}.digitalhaplotype3".format(prefix))
        medianfile = os.path.join(self.working_dir, "{}.digitalhaplotype1".format(prefix))
        ambiguityfile = os.path.join(self.working_dir, "hla.ambiguity")
        expressionfile = os.path.join(self.working_dir, "{}.expression".format(prefix))
        finaloutput = os.path.join(self.working_dir, "{}.HLAgenotype2digits".format(prefix))
        finaloutput4digit = os.path.join(self.working_dir, "{}.HLAgenotype4digits".format(prefix))
        
        self.logger.info("Starting 1st iteration:")
        self.logger.info("Mapping reads...")
        self.mapping(sam1, aligned, bowtielog, self.fq1, self.fq2, bowtiebuild, 1, mapopt)
        
        medians = {}
        for locus in locus_list:
            medians[locus] = 0
        medianflag = False

        self.logger.info("Calculation of first digital haplotype...")
        (hla_dict, readcount, readspergroup) = self.create_ref_dict(hla1fasta, locus_list)
        readcount = parse_mapping(hla_dict, sam1, readcount)
        fourDigitString1 = self.predict_HLA(sam1, medians, output1, medianflag, locus_list, hla_class, readcount, readspergroup, hla_classification, ambiguityfile)
        self.logger.info("1st iteration done.\nNow removing reads that mapped to the three top-scoring groups .......")
        
        try:
            self.remove_reads(fq1_2, fq2_2, aligned_1, aligned_2, self.create_remove_list(sam1, output1, hla_dict, locus_list))
        except IOError:
            self.logger.debug("Nothing to remove\n")
	
        self.logger.info("Starting 2nd iteration:")
        self.logger.info("Mapping reads...")
        medians = {}
        
        self.mapping(sam2, aligned, bowtielog, fq1_2, fq2_2, bowtiebuild, 2, mapopt)

        row = 2
        for locus in locus_list:
            medians[locus] = linecache.getline(medianfile, row).split('\t', 3)[2]
            row += 1
        medianflag = True

        #Calculation of second digital haplototype
        self.logger.info("Calculating second digital haplotype...")
        (hla_dict, readcount, readspergroup) = self.create_ref_dict(hla1fasta, locus_list)
        readcount = parse_mapping(hla_dict, sam2, readcount)
        fourDigitString2 = self.predict_HLA(sam2, medians, output2, medianflag, locus_list, hla_class, readcount, readspergroup, hla_classification, ambiguityfile)
        self.logger.info("2nd iteration done.")
        self.report_HLA_genotype(output1, output2, finaloutput, locus_list)
        self.logger.info("Calculating locus-specific expression...")
        try:
            self.expression(aligned_1, sam1, expressionfile, bowtielog, output1, locus_list, length_dict, hla_dict)
        except IOError:
            for locus in locus_list:
                self.logger.debug("{}: 0 RPKM".format(locus))
	
        #-----3rd iteration in case of at least one homozygous call-----------
        f1 = fourDigitString1.split(",")
        f2 = fourDigitString2.split(",")
        alleleTwoDigit_index = 0
        alleleFourDigit_index = 1
        numberSolutions_index = 2
        for locus in locus_list:
            twodigits1[locus] = f1[alleleTwoDigit_index]
            fourdigits1[locus] = f1[alleleFourDigit_index]
            if twodigits1[locus] == "no":
                fourDigit_solutions1[locus] = f1[numberSolutions_index]
            else:
                fourDigit_solutions1[locus] = int(f1[numberSolutions_index])
            fourDigit_solutions2[locus] = f2[numberSolutions_index]
            twodigits2[locus] = f2[alleleTwoDigit_index]
            fourdigits2[locus] = f2[alleleFourDigit_index]
            alleleTwoDigit_index += 3
            alleleFourDigit_index += 3
            numberSolutions_index += 3
	
        #try-catch block to prevent IO-Error in case of no expression
        try:
            if "no" in twodigits2.values():
                self.remove_reads(fq1_2, fq2_2, aligned_1, aligned_2, self.create_remove_list_four_digits(sam1, hla_dict, fourdigits1))
                
                self.mapping(sam3, aligned, bowtielog, fq1_2, fq2_2, bowtiebuild, 2, mapopt)
                readcount = parse_mapping(hla_dict, sam3, readcount)

                fourDigitString3 = self.predict_HLA(sam3, medians, output3, medianflag, locus_list, hla_class, readcount, readspergroup, hla_classification, ambiguityfile)
                f3 = fourDigitString3.split(",")
                alleleTwoDigit_index = 0
                numberSolutions_index = 2
                for locus in locus_list:
                    twodigits3[locus] = f3[alleleTwoDigit_index]
                    fourDigit_solutions3[locus] = f3[numberSolutions_index]
                    alleleTwoDigit_index += 3
                    numberSolutions_index += 3
            #1st allele--------------------------------
            for locus in locus_list:
                allele1 = fourdigits1[locus]
                if fourDigit_solutions1[locus]:
                    if int(fourDigit_solutions1[locus]) > 1:
                        allele1 += "'"
                    finalAlleles[locus] = "{}\t{}\t".format(allele1, get_best_P("{}.4digits{}1.solutions".format(sam1, locus)))
                else:
                    finalAlleles[locus] = "no\tNA\t"
	
            #2nd allele--------------------------------	
            alleleFourDigit_index = 1
            for locus in locus_list:
                allele2 = twodigits2[locus]
                if allele2 == "no":
                    allele3 = twodigits3[locus]
                    if allele3 == "no" or fourDigit_solutions1[locus] == 1 or not allele3.split(":")[0] == fourdigits1[locus].split(":"):
                        if not fourDigit_solutions1[locus]:
                            finalAlleles[locus] += "no\tNA"
                        else:
                            finalAlleles[locus] += "{}\t{}".format(fourdigits1[locus], get_P(locus, finaloutput))
                    else:
                        allele3 = f3[alleleFourDigit_index]
                        if int(f3[alleleFourDigit_index+1]) > 1:
                            allele3 += "'"
                        finalAlleles[locus] += "{}\t{}".format(allele3, get_best_P("{}.4digits{}2.solutions".format(sam3, locus)))
                        alleleFourDigit_index += 3
                else:
                    allele2 = fourdigits2[locus]
                    if int(fourDigit_solutions2[locus]) > 1:
                        allele2 += "'"
                    finalAlleles[locus] += "{}\t{}".format(allele2, get_best_P("{}.4digits{}2.solutions".format(sam2, locus)))
                    alleleFourDigit_index += 3
	
            #output the 4-digit type (stdout and to file)
            self.report_HLA_four_digit_genotype(finalAlleles, finaloutput4digit, locus_list)
        except IOError:
            #no expression"
            self.logger.info("no expression")
            for locus in locus_list:
                finalAlleles[locus] = "no\tNA\tno\tNA"
            self.report_HLA_four_digit_genotype(finalAlleles, finaloutput4digit, locus_list)

        self.cleanup()


    def cleanup(self):
        """Remove intermediate files."""
        
        files_to_keep = [
            "hla_1_classical_iteration_1.sam",
            "hla_1_classical.bowtielog",
            "hla_1_classical.HLAgenotype2digits",
            "hla_1_classical.HLAgenotype4digits",
            "hla_1_nonclassical.bowtielog",
            "hla_1_nonclassical.HLAgenotype2digits",
            "hla_1_nonclassical.HLAgenotype4digits",
            "hla_2_classical.bowtielog",
            "hla_2_classical.HLAgenotype2digits",
            "hla_2_classical.HLAgenotype4digits",
            "hla.ambiguity",
            "hla_1_classical.expression",
            "hla_1_nonclassical.expression",
            "hla_2_classical.expression",
            "run.log"
        ]

        for filename in os.listdir(self.working_dir):
            if not filename in files_to_keep:
                file_path = os.path.join(self.working_dir, filename)
                self.logger.debug("Removing {}.".format(file_path))
                os.remove(file_path)


    def mapping(self, sam, aligned_file, bowtielog, fq1, fq2, bowtiebuild, iteration, mapopt):
        """Performs the bowtie mapping for the 2 iterations using the given parameters"""
        mapping_cmd = ""
        if os.path.exists(sam):
            self.logger.info("Skipping mapping step as already done...")
            return
        if iteration == 1:
            mapping_cmd = "(bowtie -3 {} -S {} --al {} -x {} -1 {} -2 {} |  awk -F \'\\t\' '$3 != \"*\"{{ print $0 }}' > {}) 2> {}".format(self.trim3, mapopt, aligned_file, bowtiebuild, fq1, fq2, sam, bowtielog)
        elif iteration == 2:
            mapping_cmd = "bowtie -3 {} -S {} -x {} -1 {} -2 {} {}".format(self.trim3, mapopt, bowtiebuild, fq1, fq2, sam)
        #execute bowtie
        self.logger.info("CMD: {}".format(mapping_cmd))
        process_mapping = subprocess.Popen(['bash', '-c', mapping_cmd], stdout=subprocess.PIPE)
        out, err = process_mapping.communicate()

    
    def create_ref_dict(self, hlafasta, class1_list):
        """
        Create dictionary 'hla_dict', that contains all IMGT/HLA-allele names as keys
        and the allele name (e.g. A*02:01:01) as value:
        dictionary 'readcount' is initialized with 0 for each allele
        dictionary 'readspergroup' is initialized with 0 for each group (2digit, e.g. A*01)
        dictionary 'allelesPerLocus' stores the number of alleles per locus.
        """
        
        hla_dict = {}
        allelesPerLocus = {}
        readcount = {}
        readspergroup = {}
	
        for locus in class1_list:
            allelesPerLocus[locus] = 0

        with open(hlafasta) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                l = record.description.split(' ')
                hlapseudoname = l[0]
                hlaallele = l[1]
                locus = hlaallele.split('*')[0]
                if locus in class1_list:
                    hla_dict[hlapseudoname] = hlaallele
                    readcount[hlaallele] = 0
                    readspergroup[hlaallele.split(":")[0]] = 0
                    allelesPerLocus[hlaallele.split('*')[0]] += 1
        return (hla_dict, readcount, readspergroup)


    def predict_HLA(self, sam, medians, output, medianflag, locus_list, hla_class, readcount, readspergroup, hla_classification, ambiguityfile):
        """Predict HLA type"""
        maxAlleles = {}
        maxkey = {}
        allele_dict = {}
        fourdigits_dict = {}	
        alleleVector = {}
        alleleCount = {}
        readspergroup_dict_list = {}
        readspergroup_dict = {}
        fourdigits_dict = {}
        fourdigits_sorted_dict = {}

        for locus in locus_list:
            readspergroup_dict_list[locus] = []
            readspergroup_dict[locus] = {}
            fourdigits_dict[locus] = {}
            fourdigits_sorted_dict[locus] = {}
            alleleVector[locus] = []


        (top_counts, top_alleles) = get_best_allele_per_group(readspergroup, readcount)


        with open("{}.readspergrouphandle".format(sam), "a") as outf:
            for group in top_counts:
                outf.write("{}\t{}\n".format(group, top_counts[group]))
        #consider all A-,B-, and C-groups seperately
        #readspergroup<A|B|C>list = list of all reads mapping to the top-scoring groups minus the decision threshold (which is 0 in the first iteration)
        #readspergroup<A|B|C> = contains the same entries as the list, but the entries are uniquely accessible via the group-key (e.g. B*27)
        for key in top_counts:
            locus = key.split("*")[0]
            readspergroup_dict_list[locus].append(top_counts[key] - float(medians[locus]))
            readspergroup_dict[locus][key] = top_counts[key] - float(medians[locus])
	
        #Determine top-scoring group of the whole locus (A,B,C) and store it
        #maxkey<A,B,C> = group (e.g. A*02) with the most reads
        #It can be that, e.g. in cancer cells A whole locus is lost. For that reason it is checked if the 
        #number of mapping reads of the top-scoring group maxkey<A,B,C> is > 0, otherwise "no" ist reported for this locus
        for locus in locus_list:
            if len(readspergroup_dict_list[locus]) > 0:
                maxkey[locus] = max(readspergroup_dict[locus], key=lambda a:readspergroup_dict[locus].get(a))
                if readspergroup_dict[locus][maxkey[locus]] > 0:
                    maxAlleles[locus] = top_alleles[maxkey[locus]]
                    allele_dict[locus] = maxkey[locus]
                else:
                    allele_dict[locus] = "no"
            else:
                allele_dict[locus] = "no"
	
        for key in readcount:
            for locus in locus_list:	
                if key.split(":")[0] == maxkey[locus]:
                    fourdigits_dict[locus][key] = readcount[key] 
	
        for locus in locus_list:
            fourdigits_sorted_dict[locus] = sorted(fourdigits_dict[locus].items(), key=itemgetter(1), reverse=True)
	
            if medianflag:
                iteration = 2
            else:
                iteration = 1
	
            with open("{}.4digits{}{}".format(sam, locus, iteration), "w") as outf:
                for key in fourdigits_sorted_dict[locus]:
                    outf.write("{}: {}\n".format(key[0], key[1]))

            readspergrouplocus = sorted(readspergroup_dict[locus].items(), key=itemgetter(1), reverse=True)

            if medianflag:
                #in the 2nd iteration: 
                #1.) DO NOT remove the top-scoring group from the set <a,b,c>vec as this is more strict when calculating the probability of the top scoring group being an outlier
                #The strings <a,b,c>vec are used by the R-script to calculate the probability of the top scoring group being an outlier

                for key in top_alleles:
                    if key.split("*")[0] == locus:
                        alleleVector[locus].append(readcount[top_alleles[key]])
                #2.) Add the decision thresholds to the sets <a,b,c>vec, as this enables measuring the distance of the top-scoring group to this distance
                if allele_dict[locus] == "no":
                    alleleVector[locus].append(int(float(medians[locus])))
                #3.) In case of, e.g. a loss of whole HLA-locus (A,B,C), <a,b,c>vec only contain median[<0|1|2>].
                #To avoid errors in the R-Script, set to 0
                if not alleleVector[locus] or alleleVector[locus][0] == int(float(medians[locus])):
                    alleleVector[locus] = []
                    alleleCount[locus] = 0
                else:
                    alleleCount[locus] = readcount[top_alleles[maxkey[locus]]]

            else:
                #in the 1st iteration: remove the top-scoring group from the set <a,b,c>vec as this increases the certainty when calculating the probability of the top scoring group being an outlier
                for key in top_alleles:
                    if key.split("*")[0] == locus:
                        if key != maxkey[locus]:
                            alleleVector[locus].append(readcount[top_alleles[key]])
                #2.) DO NOT add the decision thresholds to the sets <a,b,c>vec
                if not alleleVector[locus]:
                    alleleCount[locus] = 0
                else:
                    alleleCount[locus] = readcount[top_alleles[maxkey[locus]]]


        confidence_vals = []
        for locus in locus_list:
            allele_list = [float(x) for x in alleleVector[locus]]
            confidence_val = calculate_confidence(float(alleleCount[locus]), allele_list)
            confidence_vals.append((maxkey[locus], confidence_val))

        fourDigitString = ""
        if medianflag:
            iteration = 2
        else:
            iteration = 1

        self.logger.info("Determining 4 digits HLA type...")
        for locus in locus_list:
            if not allele_dict[locus] == "no":
                fourDigitString += "{},{},".format(allele_dict[locus], determine_four_digits_main("{}.4digits{}{}".format(sam, locus, iteration), alleleVector[locus], hla_class, hla_classification, ambiguityfile))
            else:
                fourDigitString += "{},,,".format(allele_dict[locus])

        self.pred2file(confidence_vals, readspergroup_dict, output, allele_dict, hla_class, locus_list)
        return fourDigitString
	

    def pred2file(self, entries, readspergroup_dict, output, allele_dict, HLAclass, locus_list):
        """Write digital haplotype into file"""
        self.logger.info("Writing prediction to file (file={}).".format(output))
        with open(output, 'w') as outf:
            outf.write("HLA\tHLA1\tmedian-Value\talternative\tp-Value\n")
            index = 0
            for locus in locus_list:
                outf.write("{}\t{}\t".format(locus, allele_dict[locus]))
                #compute the decision threshold (median) for homozygosity vs. heterozygosity for the second iteration
                outf.write("{}\t".format(np.median(list(readspergroup_dict[locus].values()))))
                outf.write("{}\t{}\n".format(entries[index][0], entries[index][1]))
                index += 1

        self.logger.info("The digital haplotype is written into {}".format(output))


    def create_remove_list(self, sam, output, hla_dict, locus_list):
        """
        Open mapping file and all read ids to the list 'removeList',
        which map to one of the three groups in 'alleles'
        """
        removeList = {}
        alleles = []
        line = 2
        for locus in locus_list:
            alleles.append(linecache.getline(output, line).split('\t', 2)[1])
            line += 1

        with open(sam) as samhandle:
            for line in samhandle:
                if line[0] != '@':
                    illuminaid = line.split("\t")[0]
                    hlapseudoname = line.split("\t")[2]
                    if hla_dict[hlapseudoname].split(':')[0] in alleles:
                        removeList[illuminaid] = 1
        return removeList


    def create_remove_list_four_digits(self, sam, hla_dict, fourdigits1_dict):
        """
        Open mapping file and all read ids to the list 'removeList',
        which map to one of the three four-digit alleles in 'alleles'.
        """
        removeList = {}

        with open(sam) as samhandle:
            for line in samhandle:
                if line[0] != '@':
                    illuminaid = line.split("\t")[0]
                    hlapseudoname = line.split("\t")[2]
                    fourdigits = hla_dict[hlapseudoname].split(':')[0]+":"+hla_dict[hlapseudoname].split(':')[1]
                    if fourdigits in fourdigits1_dict:
                        removeList[illuminaid] = 1
        return removeList


    def remove_reads(self, fq1, fq2, al1, al2, removeList):
        """
        Remove reads that mapped to the three top-scoring alleles
        and write the remaining reads into two new read files.
        """
        #r1 and r2, which are the input of bowtie in the 2nd iteration
        r1 = open(fq1, "w")
        r2 = open(fq2, "w")
        #open the 2 files, that contain the reads that mapped in the 1st iteration
        aligned_handle1 = open(al1, "r")
        aligned_handle2 = open(al2, "r")
	
        #One read entry consists of 4 lines: header, seq, "+", qualities.
        for record in SeqIO.parse(aligned_handle1, "fastq"):
            #find exact id, which also appears in the mapping file
            illuminaid = record.id.split('/')[0].split(' ')[0]
            if not illuminaid in removeList:
                SeqIO.write(record, r1, "fastq")
	
        for record in SeqIO.parse(aligned_handle2, "fastq"):
            #find exact id, which also appears in the mapping file
            illuminaid = record.id.split('/')[0].split(' ')[0]
            if not illuminaid in removeList:
                SeqIO.write(record, r2, "fastq")

                
    def report_HLA_genotype(self, output1, output2, finaloutput, locus_list):
        """Write the final prediction (both digital haplotypes) to file and stdout."""
        filehandle1 = open(output1, "r").readlines()[1:len(locus_list)+1]
        filehandle2 = open(output2, "r").readlines()[1:len(locus_list)+1]
        with open(finaloutput, "w") as outfile:
            outfile.write("#Locus\tAllele 1\tConfidence\tAllele 2\tConfidence\n")
	
            for i in range(len(filehandle1)):
                filehandle1[i] = filehandle1[i][0:-1]
            for i in range(len(filehandle2)):
                filehandle2[i] = filehandle2[i][0:-1]
            self.logger.info("-----------2 digit typing results-------------")
            self.logger.info("#Locus\tAllele 1\tConfidence\tAllele 2\tConfidence")
            line = 0
            for locus in locus_list:
                allele1 = filehandle1[line].split('\t', 2)[1]
                allele1_score = filehandle1[line].split('\t')[4]
                allele2 = filehandle2[line].split('\t', 2)[1]
                if allele2 == "no":
                    allele2 = "hoz("+filehandle2[line].split('\t')[3]+")"
                allele2_score = filehandle2[line].split('\t')[4]
                line += 1
                #write complete HLA genotype to file
                outfile.write("{}\t{}\t{}\t{}\t{}\n".format(
                    locus,
                    allele1,
                    allele1_score,
                    allele2,
                    allele2_score
                ))

                #.. and print it to STDOUT
                self.logger.info("{}\t{}\t{}\t{}\t{}".format(
                    locus,
                    allele1,
                    allele1_score,
                    allele2,
                    allele2_score
                ))


    def report_HLA_four_digit_genotype(self, finalAlleles, finaloutput4digit, locus_list):

        """Write complete HLA genotype at four-digit-level to file."""
        
        with open(finaloutput4digit, "w") as outfile:
            outfile.write("#Locus\tAllele 1\tConfidence\tAllele 2\tConfidence\n")
            #.. and print it to STDOUT
            self.logger.info("-----------4 digit typing results-------------")
            self.logger.info("#Locus\tAllele 1\tConfidence\tAllele 2\tConfidence")
            for locus in locus_list:
                outfile.write("{}\t{}\n".format(locus, finalAlleles[locus]))
                self.logger.info("{}\t{}".format(locus, finalAlleles[locus]))


    def expression(self, aligned, sam, output, logfile, alleles_in, locus_list, length_dict, hla_dict):
        """Calculate locus-specific expression."""
        totalreads = float(linecache.getline(logfile, 1).split(':')[1])

        alleles = []
        line = 2
        for locus in locus_list:
            alleles.append(linecache.getline(alleles_in, line).split('\t', 2)[1])
            alleles.append(linecache.getline(alleles_in, line).split('\t')[3])
            line += 1
	
        #create read dictionary
        reads = {}
        with open(aligned) as aligned_handle1:
            for record in SeqIO.parse(aligned_handle1, "fastq"):
                illuminaid = record.id.split('/')[0].split(' ')[0]#find exact id, which also appears in the mapping file
                reads[illuminaid] = {}
                for locus in locus_list:
                    reads[illuminaid][locus] = 0
        with open(sam) as samhandle:
            for line in samhandle:
                if line[0] != '@':
                    illuminaid = line.split("\t")[0].split('/')[0].split(' ')[0]
                    hlapseudoname = line.split("\t")[2]
                    if hla_dict[hlapseudoname].split(':')[0] in alleles:
                        reads[illuminaid][hla_dict[hlapseudoname].split('*')[0]] += 1
        count = {}
        for locus in locus_list:
            count[locus] = 0
        for key in reads:
            n = 0
            for locus in reads[key]:
                if reads[key][locus] > 0:
                    n += 1
            for locus in reads[key]:
                if reads[key][locus] > 0:
                    count[locus] += float(1.0 / float(n))
	
        #Calculate RPKM and print expression values for each locus to stdout
        with open(output, 'w') as outf:
            for locus in count:
                rpkm = float((1000.0 / length_dict[locus])) * float((1000000.0 / totalreads)) * count[locus]
                self.logger.info("{}: {} RPKM".format(locus, rpkm))
                outf.write("{}: {} RPKM\n".format(locus, rpkm))


def main():
    parser = ArgumentParser(description="Predict HLA type of paired-end FASTQs")
    parser.add_argument("-1", "--fq1", dest="fq1",
                        help="Specify first FASTQ file. Gzipped input supported.",
                        required=True)
    parser.add_argument("-2", "--fq2", dest="fq2",
                        help="Specify second FASTQ file. Gzipped input supported.",
                        required=True)
    parser.add_argument("-o", "--working_dir", dest="working_dir",
                        help="Specify working dir. Results will be saved into this folder",
                        default=".")
    parser.add_argument("-p", "--threads", dest="threads", type=int,
                        help="Specify number of threads to be used by bowtie",
                        default=1)
    parser.add_argument("-3", "--trim3", dest="trim3", type=int,
                        help="Specify trimming cutoff for the low quality end. Bowtie option: -3 <int> trims <int> bases from the low quality 3' end of each read. Default: 0",
                        default=0)
    parser.add_argument("-V", "--version", action="version",
                        version="seq2HLA version {}".format(version),
                        help="show the version number and exit")


    args = parser.parse_args()

    pipe = Pipeline(args.working_dir, args.fq1, args.fq2, args.trim3, args.threads)

    pipe.run()



if __name__ == "__main__":
    main()


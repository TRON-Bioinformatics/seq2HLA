import sys, os
import numpy as np

"""
Determination of the 4-digit HLA type based on the 2-digit determination by seq2HLA.py
"""


#hoz=1


def determine4digits_main(file,bgVec,runName,HLAclass,hla_classification):
	if hla_classification=="nonclassical":
		result=determine4digitsMinorAlleles(file, bgVec, runName)
	else:
		dbmhc_file=""
		dbmhc=[]
		if HLAclass==1:
			dbmhc_file=os.path.abspath(os.path.dirname(sys.argv[0]))+"/HLAI.dbmhc"
		else:
			dbmhc_file=os.path.abspath(os.path.dirname(sys.argv[0]))+"/HLAII.dbmhc"
		dbmhc=createdbMHCdict(dbmhc,dbmhc_file)
		result=determine4digits(dbmhc, file, bgVec, runName)
	return result
	

def createdbMHCdict(dbmhc,dbmhc_file):
	population_dict={}
	linenumber=0
	for line in open(dbmhc_file,"r"):
		l=line.split("\t")
		if linenumber==0:
			#header
			for i in range(1,len(l)):
				population_dict[i]=l[i]
			linenumber=1
		else:
			if ":" in l[0]:
				fourdigit=str(l[0].split(":")[0])+":"+str(l[0].split(":")[1])
				dbmhc.append(fourdigit)
				
	return dbmhc

def meanPanPopFreq(freq_dict):
	pan_list=[]
	for entry in freq_dict:
		pan_list.append(float(freq_dict[entry]))
	return np.mean(pan_list)
	
def determine4digits(dbmhc, file, bgVec,runName):
	
	result_dict={}
	result_dict2={}
	readcount=[]
	ambiguity_dict={}
	result=[]
	allele=""
	for line in open(file,"r"):
		l=line.split(":")
		for i in range(0,len(l)-1):
			allele+=l[i]+":"
		result_dict[allele[0:-1]]=int(l[len(l)-1])
		allele=""
		readcount.append(int(l[len(l)-1]))
	
	writehandle=open(file+".solutions","w")
	writehandle.write("#Full allele\tp-value\n")
	
	readcount=[x for x in readcount if x != 0]
	
	top=readcount[0]
	readcount_copy=readcount[1:]
	readcount_copy_str=str(readcount_copy)[1:-1].replace(" ","")
	
	#call R-script "commmand_fourdigit.R" to calculate the confidence of the top-scoring allele (four-digits)
	routput=os.popen("R --vanilla < "+os.path.abspath(os.path.dirname(sys.argv[0]))+"/command_fourdigit.R --args "+str(top)+" "+str(top)+","+readcount_copy_str).read()
	parseOutput=routput.split("\n")

	entries = []
	repeat=1
	for entry in parseOutput:
		if entry[0:3]=="[1]":
			entries.append(str(entry[4:len(entry)]))
	#print entries[0]
	if not entries[0]=="NA":
		if float(entries[0]) <0.001 and float(entries[0]) > 0.0:
			for item in result_dict:
				if result_dict[item]==top:
					writehandle.write(item+"\t"+str(calculateOutlier(result_dict[item],bgVec))+"\n")
					repeat=0
					mostProbableAllele=str(item.split(":")[0])+":"+str(item.split(":")[1])
					numbersolutions=1
			
	percentile=95
	while repeat==1:
		if percentile==0:
			return "no,0"
		percentile-=5
		perc=int(np.percentile(readcount,percentile))
		if perc==0:
			perc+=0.1
		count=0
		for item in result_dict:
			if result_dict[item]>=perc:
				count+=1
		
		#there is only 1 solution:
		if count==1:
			count=0
			for item in result_dict:
				if result_dict[item]>=perc:
					writehandle.write(item+"\t"+str(calculateOutlier(result_dict[item],bgVec))+"\n")
					repeat=0
					mostProbableAllele=str(item.split(":")[0])+":"+str(item.split(":")[1])
					numbersolutions=1
		#there are more possible solutions
		else:
			count=0
			#check if allele (4-digits) occurs in dbMHC table. If not, is it a very rare allele in all populations
			for item in result_dict:
				item4digit=str(item.split(":")[0])+":"+str(item.split(":")[1])
				if result_dict[item]>=perc and item4digit in dbmhc:
					result.append(item4digit)
					result_dict2[item]=result_dict[item]
			#only 1 solution (at 4-digit level) left
			if len(set(result))==1:
				#for item in result_dict:
				for item in result_dict2:
					writehandle.write(item+"\t"+str(calculateOutlier(result_dict2[item],bgVec))+"\n")
					repeat=0
					mostProbableAllele=str(item.split(":")[0])+":"+str(item.split(":")[1])
					numbersolutions=1
			#more than one 2 solution possible (at 4-digit level).
			elif len(set(result))>1:
				numbersolutions=len(set(result))
				max=0.0
				
				for item in result_dict2:
					item4digit=str(item.split(":")[0])+":"+str(item.split(":")[1])
					if not item4digit in ambiguity_dict:
						ambiguity_dict[item4digit]=result_dict2[item]
					elif ambiguity_dict[item4digit]<result_dict2[item]:
						ambiguity_dict[item4digit]=result_dict2[item]	
					if result_dict2[item] > max:
						max=result_dict2[item]
						mostProbableAllele=item4digit
				ambiguity_handle=open(runName+".ambiguity", "a")
				ambiguity_handle.write("#################################################################################################\n")
				ambiguity_handle.write("#Ambiguity:\n")
				ambiguity_handle.write("#Based on the RNA-Seq reads and the dbMHC table, the following 4-digits alleles are possible:\n")
				for ambi4DigitAllele in ambiguity_dict:
					ambiguity_handle.write(ambi4DigitAllele+"\t"+str(calculateOutlier(ambiguity_dict[ambi4DigitAllele],bgVec))+"\n")
				
				ambiguity_handle.write("#However, by taking into account the read data, the most probable 4-digit allele is:\n")
				ambiguity_handle.write(mostProbableAllele+"\n\n")
				ambiguity_handle.close
				
				for item in result_dict2:
					item4digit=str(item.split(":")[0])+":"+str(item.split(":")[1])
					if result_dict2[item]==max:
						writehandle.write(item+"\t"+str(calculateOutlier(result_dict2[item],bgVec))+"\n")
						repeat=0
	writehandle.close()			
	return mostProbableAllele+","+str(numbersolutions)
	
def calculateOutlier(top,bgVec):
	#call R-script "commmand_fourdigit.R" to calculate the confidence of the top-scoring allele (four-digits)
	routput=os.popen("R --vanilla < "+os.path.abspath(os.path.dirname(sys.argv[0]))+"/command_fourdigit.R --args "+str(top)+" "+bgVec).read()
	parseOutput=routput.split("\n")
	
	entries = []
	repeat=1
	for entry in parseOutput:
		if entry[0:3]=="[1]":
			entries.append(str(entry[4:len(entry)]))
	return entries[0]

#method for determining the 4 digit HLA type for minor HLA alleles, for which no entry in dbMHC is available
def determine4digitsMinorAlleles(file, bgVec,runName):
	
	result_dict={}
	result_dict2={}
	readcount=[]
	ambiguity_dict={}
	result=[]
	allele=""
	for line in open(file,"r"):
		l=line.split(":")
		for i in range(0,len(l)-1):
			allele+=l[i]+":"
		result_dict[allele[0:-1]]=int(l[len(l)-1])
		allele=""
		readcount.append(int(l[len(l)-1]))
	
	writehandle=open(file+".solutions","w")
	writehandle.write("#Full allele\tp-value\n")
	
	readcount=[x for x in readcount if x != 0]
	
	top=readcount[0]
	readcount_copy=readcount[1:]
	readcount_copy_str=str(readcount_copy)[1:-1].replace(" ","")
	
	#call R-script "commmand_fourdigit.R" to calculate the confidence of the top-scoring allele (four-digits)
	routput=os.popen("R --vanilla < "+os.path.abspath(os.path.dirname(sys.argv[0]))+"/command_fourdigit.R --args "+str(top)+" "+str(top)+","+readcount_copy_str).read()
	parseOutput=routput.split("\n")

	entries = []
	repeat=1
	for entry in parseOutput:
		if entry[0:3]=="[1]":
			entries.append(str(entry[4:len(entry)]))
	#print entries[0]
	if not entries[0]=="NA":
		if float(entries[0]) <0.001 and float(entries[0]) > 0.0:
			for item in result_dict:
				if result_dict[item]==top:
					writehandle.write(item+"\t"+str(calculateOutlier(result_dict[item],bgVec))+"\n")
					repeat=0
					mostProbableAllele=str(item.split(":")[0])+":"+str(item.split(":")[1])
					numbersolutions=1
			
	percentile=95
	while repeat==1:
		if percentile==0:
			return "no,0"
		percentile-=5
		perc=int(np.percentile(readcount,percentile))
		if perc==0:
			perc+=0.1
		count=0
		for item in result_dict:
			if result_dict[item]>=perc:
				count+=1
		
		#there is only 1 solution:
		if count==1:
			count=0
			for item in result_dict:
				if result_dict[item]>=perc:
					writehandle.write(item+"\t"+str(calculateOutlier(result_dict[item],bgVec))+"\n")
					repeat=0
					mostProbableAllele=str(item.split(":")[0])+":"+str(item.split(":")[1])
					numbersolutions=1
		#there are more possible solutions
		else:
			count=0
			#check if allele (4-digits) is in the defined percentile range
			for item in result_dict:
				item4digit=str(item.split(":")[0])+":"+str(item.split(":")[1])
				if result_dict[item]>=perc:
					result.append(item4digit)
					result_dict2[item]=result_dict[item]
			#only 1 solution (at 4-digit level) left
			if len(set(result))==1:
				#for item in result_dict:
				for item in result_dict2:
					writehandle.write(item+"\t"+str(calculateOutlier(result_dict2[item],bgVec))+"\n")
					repeat=0
					mostProbableAllele=str(item.split(":")[0])+":"+str(item.split(":")[1])
					numbersolutions=1
			#more than one 2 solution possible (at 4-digit level).
			elif len(set(result))>1:
				numbersolutions=len(set(result))
				max=0.0
				
				for item in result_dict2:
					item4digit=str(item.split(":")[0])+":"+str(item.split(":")[1])
					if not item4digit in ambiguity_dict:
						ambiguity_dict[item4digit]=result_dict2[item]
					elif ambiguity_dict[item4digit]<result_dict2[item]:
						ambiguity_dict[item4digit]=result_dict2[item]	
					if result_dict2[item] > max:
						max=result_dict2[item]
					#if meanPanPopFreq(dbmhc_prob[item4digit]) > max:
					#	max = meanPanPopFreq(dbmhc_prob[item4digit])
						mostProbableAllele=item4digit
					#if np.mean(dbmhc_prob[item4digit]) > max:
					#	max = np.mean(dbmhc_prob[item4digit])
					#	mostProbableAllele=item4digit
				ambiguity_handle=open(runName+".ambiguity", "a")
				ambiguity_handle.write("#################################################################################################\n")
				ambiguity_handle.write("#Ambiguity:\n")
				ambiguity_handle.write("#Based on the RNA-Seq reads and the dbMHC table, the following 4-digits alleles are possible:\n")
				for ambi4DigitAllele in ambiguity_dict:
					ambiguity_handle.write(ambi4DigitAllele+"\t"+str(calculateOutlier(ambiguity_dict[ambi4DigitAllele],bgVec))+"\n")
				
				ambiguity_handle.write("#However, by taking into account the read data, the most probable 4-digit allele is:\n")
				ambiguity_handle.write(mostProbableAllele+"\n\n")
				ambiguity_handle.close
				
				for item in result_dict2:
					item4digit=str(item.split(":")[0])+":"+str(item.split(":")[1])
					#if np.mean(dbmhc_prob[item4digit]) == max:
					#if meanPanPopFreq(dbmhc_prob[item4digit]) == max:
					if result_dict2[item]==max:
						writehandle.write(item+"\t"+str(calculateOutlier(result_dict2[item],bgVec))+"\n")
						repeat=0
	writehandle.close()			
	return mostProbableAllele+","+str(numbersolutions)	

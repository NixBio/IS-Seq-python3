import sys
import os
import re
import seaborn as sns
import pandas as pd
import numpy as np
from multiprocessing import Pool
from scipy import stats
import matplotlib.pyplot as plt
import pysam


def processBam(bam):
	'''
	Processing function: calls pool of worker functions
	to extract from a bam file two definitions of the edit distances to the reference genome scaled by read length
	Returned in a pandas DataFrame
	'''
	samfile = pysam.AlignmentFile(bam, "rb")
	if not samfile.has_index():
		pysam.index(bam)
		samfile = pysam.AlignmentFile(bam, "rb")  # Need to reload the samfile after creating index
	chromosomes = samfile.references
	datadf = pd.DataFrame()
	pool = Pool(processes=12)
	params = zip([bam]*len(chromosomes), chromosomes)
	try:
		output = [results for results in pool.imap(extractFromBam, params)]
	except KeyboardInterrupt:
		print("Terminating worker threads")
		pool.terminate()
		pool.join()
		sys.exit()
	datadf["editDistancesNM"] = np.array([x for y in [elem[0] for elem in output] for x in y])
	datadf["editDistancesMD"] = np.array([x for y in [elem[1] for elem in output] for x in y])
	return datadf

def processBam1(bam):
	'''
	Processing function: calls pool of worker functions
	to extract from a bam file two definitions of the edit distances to the reference genome scaled by read length
	Returned in a pandas DataFrame
	'''
	samfile = pysam.AlignmentFile(bam, "rb")
	if not samfile.has_index():
		pysam.index(bam)
		samfile = pysam.AlignmentFile(bam, "rb")  # Need to reload the samfile after creating index
	chromosomes = samfile.references
	datadf = pd.DataFrame()
	pool = Pool(processes=12)
	params = zip([bam]*len(chromosomes), chromosomes)
	try:
		output = [results for results in pool.imap(extractFromBam1, params)]
	except KeyboardInterrupt:
		print("Terminating worker threads")
		pool.terminate()
		pool.join()
		sys.exit()
	datadf["readName"] = np.array([x for y in [elem[0] for elem in output] for x in y])
	datadf["editDistancesNM"] = np.array([x for y in [elem[1] for elem in output] for x in y])
	datadf["readLengthNM"] = np.array([x for y in [elem[2] for elem in output] for x in y])
	datadf["editDistancesMD"] = np.array([x for y in [elem[3] for elem in output] for x in y])
	datadf["readLengthMD"] = np.array([x for y in [elem[4] for elem in output] for x in y])
	return datadf

def extractFromBam(params):
	'''
	Worker function per chromosome
	loop over a bam file and create tuple with lists containing metrics:
	two definitions of the edit distances to the reference genome scaled by aligned read length
	'''
	bam, chromosome = params
	samfile = pysam.AlignmentFile(bam, "rb")
	editDistancesNM = []
	editDistancesMD = []
	for read in samfile.fetch(reference=chromosome, multiple_iterators=True):
		editDistancesNM.append(read.get_tag("NM")/read.query_alignment_length)
		editDistancesMD.append(
			(sum([len(item) for item in re.split('[0-9^]', read.get_tag("MD"))]) +  # Parse MD string to get mismatches/deletions
			sum([item[1] for item in read.cigartuples if item[0] == 1]))  # Parse cigar to get insertions
			/read.query_alignment_length)
	return (editDistancesNM, editDistancesMD)

def extractFromBam1(params):
	'''
	Worker function per chromosome
	loop over a bam file and create tuple with lists containing metrics:
	two definitions of the edit distances to the reference genome scaled by aligned read length
	'''
	bam, chromosome = params
	samfile = pysam.AlignmentFile(bam, "rb")
	readName = []
	editDistancesNM = []
	readLengthNM = []
	editDistancesMD = []
	readLengthMD = []
	for read in samfile.fetch(reference=chromosome, multiple_iterators=True):
		readName.append(read.query_name)
		editDistancesNM.append(read.get_tag("NM"))
		readLengthNM.append(read.query_alignment_length)
		editDistancesMD.append(
			(sum([len(item) for item in re.split('[0-9^]', read.get_tag("MD"))]) +  # Parse MD string to get mismatches/deletions
			sum([item[1] for item in read.cigartuples if item[0] == 1])))  # Parse cigar to get insertions
		readLengthMD.append(read.query_alignment_length)
	return (readName,editDistancesNM,readLengthNM,editDistancesMD,readLengthMD)

def makePlot(datadf):
	try:
		plot = sns.jointplot(
        	x='editDistancesNM',
	        y='editDistancesMD',
	        data=datadf,
	        kind="kde",
	        color="#4CB391",
	        stat_func=stats.pearsonr,
	        space=0,
	        size=10)
		plot.savefig('EditDistancesCompared_kde.png', format='png', dpi=1000)
	except:  # throws an error with perfect correlation! 
		pass
	plot = sns.jointplot(
		x='editDistancesNM',
		y='editDistancesMD',
		data=datadf,
		kind="scatter",
		color="#4CB391",
		stat_func=stats.pearsonr,
		space=0,
		joint_kws={"s": 1},
		size=10)
	plot.savefig('EditDistancesCompared_scatter.png', format='png', dpi=1000)


df = processBam(sys.argv[1])
print(df)

fileName=os.path.basename(sys.argv[1])

print(fileName)

sampleName=fileName.rsplit('.', maxsplit=1)[0]

print(fileName)
fileName='/home/local_scratch/ISseqOutput/230221M086400005000000000KWN5YRerun/'+sampleName+'_edit_distance.csv'

df.to_csv(fileName) 
#makePlot(df)

df1 = processBam1(sys.argv[1])

print(df1)

fileName1='/home/local_scratch/ISseqOutput/230221M086400005000000000KWN5YRerun/'+sampleName+'_edit_distance_1.csv'
print(fileName1)

df1.to_csv(fileName1) 


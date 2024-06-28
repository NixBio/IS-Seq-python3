# This script takes as input the "final_parse" file of each sample and applies 7 bp window rule to merge IS within one sample 
# by considering LTR direction(strand: +/-)

# Input: POOL-ISA-AVRO-6-Preclin_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_final_parse_filterNo.txt or 
# POOL-ISA-AVRO-6-Preclin_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_final_parse_filter60.txt
# add as input the 
# Input dir is Cutadapt

# Output: Dir: Diversity out

##################### Todo:#################################
## Properly hand over the Index information!!! currently handling is worng

# *final_parse_filter*_report.txt
# *filter*_NonGrouped.txt
# *final_parse_filter*_grouped_IS.txt

import sys
import operator
import re
import os
import subprocess
import numpy as np
from itertools import groupby
import fnmatch

def unique(seq):
    # order preserving
    noDupes = []
    [noDupes.append(i) for i in seq if not noDupes.count(i)]
    return noDupes


## orig
# CutAdapt/TCRseq-ISA-testseq-1_FP-ISA-i5-B2-LTR_FP-ISA-i7-B2-LC_final_parse_filter60.txt
myblast=open(sys.argv[1]) 
# /reference/hg38/hg38.genome.num.txt
my_legend=open(sys.argv[2])
# ISout, summary per chomosome location..
dwdt=sys.argv[3]

## troubleshooting input
#myblast=open("/out/ISSeqout_20240626/CutAdapt/"TCRseq-ISA-testseq-1_FP-ISA-i5-B2-LTR_FP-ISA-i7-B2-LC_final_parse_filter60.txt")
#my_legend=open("/reference/hg38/hg38.genome.num.txt")
#dwdt="/out/ISSeqout_20240626/CutAdapt/ISout"


os.chdir(dwdt)

baseDir = os.path.dirname(sys.argv[1])
baseFile = os.path.basename(sys.argv[1])
outDir = os.path.join(baseDir,"DiversityOut")

if not os.path.exists(outDir):
    os.makedirs(outDir)

## Define output files
my_report=open(os.path.join(outDir,baseFile[:-4]+'_report.txt'),'w')
grouped_IS=open(os.path.join(outDir,baseFile[:-4]+'_grouped_IS.txt'),'w')
NonGrouped=open(os.path.join(outDir,baseFile[:-4]+'_NonGrouped.txt'),'w')


## Define variables

list_pos=[] # all input file paths
list_strand=[]
list_diff=[]
list_item1=[]
list_item3=[]
final_list=[]
mydict={}
myleg={} # chromosome information
myData={}
my_pos={}
dist_cutoff = 7



#List all the file containing the reads per IS
for f in os.listdir(dwdt):
    if fnmatch.fnmatch(f,'*LTR*LC*filter*_chr*.txt'):
        list_pos.append(f)

# this are all basenames of the input files
list_posbase = [os.path.basename(x) for x in list_pos]
#list_posbase = [x.rstrip(".txt") for x in list_posbase]

#Read the input file (final_parse) and store the infos about strand, chr and pos
for l in myblast.readlines():
    l= l.rstrip()
    if l[0]=='#':
        continue
    l= l.split('\t')
    tmp= l[4], l[5]
    try:
        mydict[tmp].append(l[0])
    except:
        mydict[tmp]=[l[0]]
#
for l in my_legend.readlines():
    l= l.rstrip()
    if l[0]=='#':
        continue
    l= l.split('\t')
    try:
        myleg[l[1]].append(l[0])
    except:
        myleg[l[1]]=[l[0]]



#Defines the correspondence between the input file and the files into the listi
# The following code assumes there is a "_" between the different functional parts of the file name, in addition it assumes that all barcodes have to be named in a certain fashion: ID.1, ID.2, while we have ...B2..,D2,C2

# Input: POOL-ISA-AVRO-6-Preclin_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_final_parse_filterNo.txt or 
# fields: for below: POOL-ISA-AVRO-6-Preclin, FB-P5-Rd1-LTR, 9, FB-P7-Rd2-LC, 9, final, parse, filterNo, txt 
# best practice, there is a dictionary of Adaptercombination per subsample.
# 
baseFile=baseFile.rstrip(".txt")
baseFileS=re.split('\_', baseFile)
if len(baseFileS) > 6:
    print("there are underscores used in the samplename and or Barcode names!, this is forbidden")
    sys.exit(1)

baseFileLB=re.split('\.|\_', baseFile)
if len(baseFileLB)==8:
    # for barcode containing a [.] in the name followedby an increasing number 
    #       4: "9"
    # 7: filterNo
    print(baseFileLB[2]+'\t'+baseFileLB[4]+'\t'+baseFileLB[7]+'\n')
    for element in list_posbase:
        element=re.split('\.|\_', element)
        #Finds the correct LTR and LC barcode
        if element[2]==baseFileLB[2] and element[4]==baseFileLB[4]:
            #Finds the correct filter: filterNo, filter30, filter45, filter60
            if element[5]==baseFileLB[7]:
                print(element[2]+'\t'+element[4]+'\t'+element[5]+'__'+baseFileLB[2]+'\t'+baseFileLB[4]+'\t'+baseFileLB[7]+'\n')
                for chr,pos in mydict:
                    #if the IS falls into a standard chr
                    if len(element)==9:
                        if chr == element[6] and pos == element[7]:
                            if str(set(mydict[chr,pos]))[1:9] == "'R1_rev'":
                                for name in myleg:
                                    if element[6]==name:
                                        for code in myleg[name]:
                                            #print chr, pos, code, '-', str(len(mydict[chr,pos]))
                                            list_strand.append([chr, pos, code, '-', str(len(mydict[chr,pos]))])
                            else:
                                for name in myleg:
                                    if element[6]==name:
                                        for code in myleg[name]:
                                            #print chr, pos, code, '-', str(len(mydict[chr,pos]))
                                            list_strand.append([chr, pos, code, '+', str(len(mydict[chr,pos]))])
                    #if the IS falls into a chr_something
                    if len(element)==10:
                        if chr == (element[6]+'_'+element[7]) and pos == element[8]:
                            if str(set(mydict[chr,pos]))[1:9] == "'R1_rev'":
                                for name in myleg:
                                    if (element[6]+'_'+element[7])==name:
                                        for code in myleg[name]:
                                            #print chr, pos, code, '-', str(len(mydict[chr,pos]))
                                            list_strand.append([chr, pos, code, '-', str(len(mydict[chr,pos]))])
                            else:
                                for name in myleg:
                                    if (element[6]+'_'+element[7])==name:
                                        for code in myleg[name]:
                                            #print chr, pos, code, '+', str(len(mydict[chr,pos]))
                                            list_strand.append([chr, pos, code, '+', str(len(mydict[chr,pos]))])
                    #if the IS falls into a chr_something_something
                    if len(element)==11:
                        if chr == (element[6]+'_'+element[7]+'_'+element[8]) and pos == element[9]:
                            if str(set(mydict[chr,pos]))[1:9] == "'R1_rev'":
                                for name in myleg:
                                    if (element[6]+'_'+element[7]+'_'+element[8])==name:
                                        for code in myleg[name]:
                                            #print chr, pos, code, '-', str(len(mydict[chr,pos]))
                                            list_strand.append([chr, pos, code, '-', str(len(mydict[chr,pos]))])
                            else:
                                for name in myleg:
                                    if (element[6]+'_'+element[7]+'_'+element[8])==name:
                                        for code in myleg[name]:
                                            #print chr, pos, code, '+', str(len(mydict[chr,pos]))
                                            list_strand.append([chr, pos, code, '+', str(len(mydict[chr,pos]))])
############################
# For barcodes not containing a [.] in the name, differences are in the whole name
########
elif len(baseFileS)==6:
    print(baseFileS[1]+'\t'+baseFileS[2]+'\t'+baseFileS[5]+'\n')
    for element in list_posbase:
        element=re.split('\.|\_', element)
        #Finds the correct LTR and LC barcode
        if element[1]==baseFileS[1] and element[2]==baseFileS[2]:
            #Finds the correct filter: filterNo, filter30, filter45, filter60
            if element[3]==baseFileS[5]:
                #print(element[2]+'\t'+element[4]+'\t'+element[5]+'__'+baseFile[2]+'\t'+baseFile[4]+'\t'+baseFile[7]+'\n')
                print(element[1]+'\t'+element[2]+'\t'+element[3]+'__'+baseFileS[1]+'\t'+baseFileS[2]+'\t'+baseFileS[5]+'\n')
                for chr,pos in mydict:
                    #if the IS falls into a standard chr
                    if len(element)==7:
                        if chr == element[4] and pos == element[5]:
                            if str(set(mydict[chr,pos]))[1:9] == "'R1_rev'":
                                for name in myleg:
                                    if element[4]==name:
                                        for code in myleg[name]:
                                            #print chr, pos, code, '-', str(len(mydict[chr,pos]))
                                            list_strand.append([chr, pos, code, '-', str(len(mydict[chr,pos]))])
                            else:
                                for name in myleg:
                                    if element[4]==name:
                                        for code in myleg[name]:
                                            #print chr, pos, code, '-', str(len(mydict[chr,pos]))
                                            list_strand.append([chr, pos, code, '+', str(len(mydict[chr,pos]))])
                    #if the IS falls into a chr_something
                    if len(element)==8:
                        if chr == (element[4]+'_'+element[5]) and pos == element[6]:
                            if str(set(mydict[chr,pos]))[1:9] == "'R1_rev'":
                                for name in myleg:
                                    if (element[4]+'_'+element[5])==name:
                                        for code in myleg[name]:
                                            #print chr, pos, code, '-', str(len(mydict[chr,pos]))
                                            list_strand.append([chr, pos, code, '-', str(len(mydict[chr,pos]))])
                            else:
                                for name in myleg:
                                    if (element[4]+'_'+element[5])==name:
                                        for code in myleg[name]:
                                            #print chr, pos, code, '+', str(len(mydict[chr,pos]))
                                            list_strand.append([chr, pos, code, '+', str(len(mydict[chr,pos]))])
                    #if the IS falls into a chr_something_something
                    if len(element)==9:
                        if chr == (element[4]+'_'+element[5]+'_'+element[6]) and pos == element[7]:
                            if str(set(mydict[chr,pos]))[1:9] == "'R1_rev'":
                                for name in myleg:
                                    if (element[4]+'_'+element[5]+'_'+element[6])==name:
                                        for code in myleg[name]:
                                            #print chr, pos, code, '-', str(len(mydict[chr,pos]))
                                            list_strand.append([chr, pos, code, '-', str(len(mydict[chr,pos]))])
                            else:
                                for name in myleg:
                                    if (element[4]+'_'+element[5]+'_'+element[6])==name:
                                        for code in myleg[name]:
                                            #print chr, pos, code, '+', str(len(mydict[chr,pos]))
                                            list_strand.append([chr, pos, code, '+', str(len(mydict[chr,pos]))])
#Sorts the list by chr and pos
for item in sorted(list_strand, key=lambda x: (int(x[2]), int(x[1]))):
    try:
        myData[item[2]].append(int(item[1]))
    except:
        myData[item[2]]=[int(item[1])]
#
for chr in myData:
    if len(myData[chr])>=2:
        count=0
        for count in range(0,(len(myData[chr])-1)):
            next_count=count+1
            if count==next_count:
                continue
            else:
                #compare the i element of this list with the next one
                myDiff= int(sorted(myData[chr])[count]) - int(sorted(myData[chr])[next_count])
                if abs(myDiff)<= dist_cutoff:
                    #if the distance between two element of this list is lower than the cut_off, make another list
                    list_diff.append([chr, myData[chr][count], myData[chr][next_count]])
#
for item in sorted(list_strand, key=lambda x: (int(x[2]), int(x[1]))):
    control= True
    for r in sorted(list_diff):
        if int(item[2])==int(r[0]) and int(item[1])==int(r[1]):
            list_item1.append([item[0], item[1], item[2], item[3], int(item[4])])
            control= False
        if int(item[2])==int(r[0]) and int(item[1])==int(r[2]):
            list_item1.append([item[0], item[1], item[2], item[3], int(item[4])])
            control= False
    if control:
        list_item3.append([item[0], item[1], item[2], item[3], int(item[4])])
#
unique_list=unique(list_item1)
for i in unique_list:
    tmp=i[0],i[1],i[2],i[3],i[4]
    try:
        my_pos[int(i[2])].append(tmp)
    except:
        my_pos[int(i[2])]=[tmp]

for a in my_pos:
    b=my_pos[a]
    d=sorted(b, key=lambda x: (int(x[2]), int(x[1])))
    diff = [int(d[i+1][1])-int(d[i][1]) for i in range(len(d)-1)]
    avg = sum(diff) / len(diff)
    c=[[d[0][0]]]
    m=[[d[0][1]]]
    n=[[d[0][4]]]
    strand=[[d[0][3]]]
    for i in range(0,len(d)-1):
        if int(d[i+1][1])-int(d[i][1]) <= dist_cutoff and d[i+1][3]==d[i][3]:
            c[-1].append(d[i+1][0])
            m[-1].append(d[i+1][1])
            n[-1].append(d[i+1][4])
            strand[-1].append(d[i+1][3])
        else:
            c.append([d[i+1][0]])
            m.append([d[i+1][1]])
            n.append([d[i+1][4]])
            strand.append([d[i+1][3]])
    c=np.array(c,dtype=object)
    m=np.array(m,dtype=object)
    n=np.array(n,dtype=object)
    strand=np.array(strand,dtype=object)
    #Check if this script is trying to merge different integration sites together.
    #Opposite strand will not be part of the same integration site and you will get the warnining message.
    for segno in strand:
        for contatore in range(0,len(segno)-1):
            if segno[contatore+1]==segno[contatore]:
                continue
            else:
                print('Warninig: The merged sites contain opposite strands.')
                print('Please check the report and modify the grouped_IS file, if necessary.')
    print_array= np.dstack((c,m,n,strand))
    my_array= np.vstack(([c],[m],[n],[strand]))
    #Print a report in which the reads of the same integration sites are grouped together
    for i in my_array[[2],:]:
        count=0
        for x in i:
            x= np.array(x,dtype=np.int64)
            reads_count= sum(x)
            l= x.argmax()
            final_list.append([str(my_array[0][count][l]), str(my_array[1][count][l]), str(d[0][2]), str(my_array[3][count][l]), str(reads_count)])
            if len(my_array[0][count])>=2:
                print('>'+str(my_array[0][count][l])+'_'+str(my_array[1][count][l])+'_'+str(d[0][2])+'_'+str(my_array[3][count][l])+'_'+str(reads_count))
                my_report.write('>'+str(my_array[0][count][l])+'_'+str(my_array[1][count][l])+'_'+str(d[0][2])+'_'+str(my_array[3][count][l])+'_'+str(reads_count)+'\n')
                for valori in range(len(my_array[0][count])):
                    my_report.write(str(my_array[0][count][valori])+'\t'+str(my_array[1][count][valori])+'\t'+str(d[0][2])+'\t'+str(my_array[3][count][valori])+'\t'+str(my_array[2][count][valori])+'\n')
                    print(my_array[0][count][valori], my_array[1][count][valori], d[0][2], my_array[3][count][valori], my_array[2][count][valori])
            count+=1
#
for c in list_item3:
    NonGrouped.write(str(c[0])+'\t'+str(c[1])+'\t'+str(c[2])+'\t'+str(c[3])+'\t'+str(c[4])+'\n')
    final_list.append([str(c[0]), str(c[1]), str(c[2]), str(c[3]), str(c[4])])

for y in sorted(final_list, key=lambda x: (int(x[2]), int(x[1]))):
    grouped_IS.write('\t'.join(y)+'\n')
#
myblast.close()
my_legend.close()
my_report.close()
#

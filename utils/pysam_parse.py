# This script uses pysam to parse bam alignment and select only those reads showing an exact match in the 3 nt following the end of the LTR, which ensures the insertion site is properly detected. 
# This filter only applies to the read that is first in pair (R1) because these where carring the LTR before alignment.

# input:  R1_R2_Barcode_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_aligned_mem_allFilter_rehead.bam
# output: R1_R2_Barcode_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_aligned_mem_allFilter_rehead_exact3nt.bam

import pysam
import sys

samfile = pysam.AlignmentFile(sys.argv[1],'rb',check_sq=False)
out_reads = pysam.AlignmentFile(sys.argv[1][:-4]+'_exact3nt.bam', 'wb', template=samfile)
for read in samfile.fetch():
    #if the read is aligned in reverse
    if read.is_reverse:
        #if the read is first in pair (R1)
        if read.is_read1:
            for i in read.cigar:
                if read.cigar.index(i) == (len(read.cigar)-1):
                    if i[0] == 0 and i[1] >= 3:
                        #if the first 3 nt are all an exact match in MD tag
                        try:
                            if int(read.get_tag('MD')[(len(read.get_tag('MD'))-1)])>=3:
                                out_reads.write(read)
                                break
                            else:
                                if int(read.get_tag('MD')[(len(read.get_tag('MD'))-2)])>=0:
                                    out_reads.write(read)
                                    break
                        except:
                            continue
        #if the read is the second in pair (R2) there is no need for filtering
        else:
            out_reads.write(read)
    #if the read is aligned in forward
    else:
        #if the read is first in pair (R1)
        if read.is_read1:
            for i in read.cigar:
                #if the first 3 nt are a match in cigar
                if read.cigar.index(i) == 0:
                    if i[0] == 0 and i[1] >= 3:
                        #if the first 3 nt are all an exact match in MD tag
                        try:
                            if int(read.get_tag('MD')[0])>=3:
                                out_reads.write(read)
                                break
                            else:
                                if int(read.get_tag('MD')[1])>=0:
                                    out_reads.write(read)
                                    break
                        except:
                            continue
        #if the read is the second in pair (R2) there is no need for filtering
        else:
            out_reads.write(read)
out_reads.close()
samfile.close()

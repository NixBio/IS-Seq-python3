
nohup docker run --rm -v path/to/IS-Seq-python3/utilsRefData/IsSeq:/out aiminy/isseq:2.4 Rscript /usr/src/IS-Seq-python3/R/makeREFIndex1.R -i https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/GRCh38.primary_assembly.genome.fa.gz -g https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.annotation.gtf.gz -r https://hgdownload.soe.ucsc.edu/goldenPath/currentGenomes/Homo_sapiens/database/rmsk.txt.gz -m https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/chromInfo.txt.gz -o /out/hg38/GRCh38.primary_assembly.genome.fa > path/to/ISseqOutput/log/logMakeHg38.txt 2>&1 &

wait

nohup docker run --rm -v path/to/IS-Seq-python3/data:/in --rm -v  path/to/IS-Seq-python3/sample_research:/in1 --rm -v path/to/IS-Seq-python3/utilsRefData/IsSeq:/in2 --rm -v path/to/ISseqOutput:/out aiminy/isseq:2.4 python /usr/src/IS-Seq-python3/ISpipelineFv3_test.py -1 /in/simulationUp_R1.fq.gz -2 /in/simulationUp_R2.fq.gz -s POOL-ISA-AVRO-6-Preclin -o /out -t DEMO -r /in1/20210121_AssociationFIle_POOL6_Preclinical.csv -u /in2 -p /usr/src/IS-Seq-python3/utils -a read -c nothing -q 30 > path/to/ISseqOutput/log/logDEMO_read.txt 2>&1 &

wait

nohup docker run --rm -v path/to/IS-Seq-python3/data:/in --rm -v  path/to/IS-Seq-python3/sample_research:/in1 --rm -v path/to/IS-Seq-python3/utilsRefData/IsSeq:/in2 --rm -v path/to/ISseqOutput:/out aiminy/isseq:2.4 python /usr/src/IS-Seq-python3/ISpipelineFv3_test.py -1 /in/simulationUp_R1.fq.gz -2 /in/simulationUp_R2.fq.gz -s POOL-ISA-AVRO-6-Preclin -o /out -t DEMO -r /in1/20210121_AssociationFIle_POOL6_Preclinical.csv -u /in2 -p /usr/src/IS-Seq-python3/utils -a umi -c nothing -q 30 > path/to/ISseqOutput/log/logDEMO_umi.txt 2>&1 &

wait

nohup docker run --rm -v path/to/IS-Seq-python3/data:/in --rm -v  path/to/IS-Seq-python3/sample_research:/in1 --rm -v path/to/IS-Seq-python3/utilsRefData/IsSeq:/in2 --rm -v path/to/ISseqOutput:/out aiminy/isseq:2.4 python /usr/src/IS-Seq-python3/ISpipelineFv3_test.py -1 /in/simulationUp_R1.fq.gz -2 /in/simulationUp_R2.fq.gz -s POOL-ISA-AVRO-6-Preclin -o /out -t DEMO -r /in1/20210121_AssociationFIle_POOL6_Preclinical.csv -u /in2 -p /usr/src/IS-Seq-python3/utils -a fragment -c nothing -q 30 > path/to/ISseqOutput/log/logDEMO_Frag.txt 2>&1 &
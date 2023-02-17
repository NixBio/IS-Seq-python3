# docker build --no-cache -f $HOME/IS-Seq/IS-Seq-python3/Docker/Dockerfile -t aiminy/isseq:2.5 .

From ubuntu:focal

MAINTAINER Aimin Yan <aimin.at.work@gmail.com>

ENV DEBIAN_FRONTEND=noninteractive


RUN apt-get update && apt-get -y upgrade && \
	apt-get install -y lsb-release build-essential wget curl git default-jdk ant && \
	apt-get	install -y libssl-dev libtbb2 libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev libcurl3-dev libncurses5 && \
	apt-get clean && apt-get purge && \
	rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*


WORKDIR /

RUN wget https://github.com/compgen-io/ngsutilsj/releases/download/ngsutilsj-0.4.17/ngsutilsj
RUN chmod a+x /ngsutilsj 

ENV PATH "/":${PATH}
RUN echo "${PATH}" >> .bashrc

# Will copy from existing Docker image
COPY --from=continuumio/miniconda3:4.12.0 /opt/conda /opt/conda

ENV PATH=/opt/conda/bin:$PATH

# install the tools from bioconda
RUN set -ex && \
    conda config --set always_yes yes --set changeps1 no && \
    conda info -a && \
    conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda install --quiet --freeze-installed -c main conda-pack && \
    conda install -c bioconda fastx_toolkit && \
    conda install -c bioconda seqtk && \ 
    conda install -c bioconda -y ea-utils && \
    conda install -c bioconda -y bwa-mem2 && \
    conda install -c bioconda -y blat && \
    conda install -c bioconda -y starcode && \
    conda install -c bioconda samtools && \
    conda install -c bioconda art

# install flexbar

WORKDIR /

RUN wget https://github.com/seqan/flexbar/releases/download/v3.5.0/flexbar-3.5.0-linux.tar.gz && \
        tar xvf flexbar-3.5.0-linux.tar.gz && \
        rm flexbar-3.5.0-linux.tar.gz

ENV PATH "/flexbar-3.5.0-linux":${PATH}
RUN echo "${PATH}" >> .bashrc

# install bedtools picard-tools cutadapt ncbi-blast+(for using blastn)

RUN apt-get update && apt-get -y upgrade && apt install -y bedtools picard-tools cutadapt ncbi-blast+ 

RUN pip install biopython numpy more_itertools Bio

# install umi_tools
RUN wget https://github.com/CGATOxford/UMI-tools/archive/refs/tags/1.1.2.tar.gz && \
    tar xvf 1.1.2.tar.gz && \
    rm 1.1.2.tar.gz && \
    cd UMI-tools-1.1.2 && \
    python setup.py install --user
ENV PATH=$HOME/root/.local/bin/:$PATH

RUN echo "${PATH}" >> .bashrc

WORKDIR /usr/src

#IS-Seq

RUN git clone https://github.com/aiminy/IS-Seq-python3

RUN apt update -qq && \
    apt install -y --no-install-recommends software-properties-common dirmngr && \
    wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc && \
    add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" && \
    apt install -y --no-install-recommends r-base

RUN Rscript -e "install.packages(c('reshape'),lib='/usr/lib/R/site-library')"

RUN Rscript -e "install.packages(c('reshape2'),lib='/usr/lib/R/site-library')"

RUN Rscript -e "install.packages(c('stringr','ggplot2','plyr','optparse','curl','IRanges','data.table','igraph','tidysq'),lib='/usr/lib/R/site-library')"

RUN apt-get update && \
    apt-get install libxml2-dev

RUN Rscript -e "install.packages(c('XML'),lib='/usr/lib/R/site-library')"

RUN Rscript -e "install.packages(c('restfulr'),lib='/usr/lib/R/site-library')"

RUN Rscript -e "install.packages(c('sonicLength'),lib='/usr/lib/R/site-library')"

RUN Rscript -e "install.packages(c('tidyverse'),lib='/usr/lib/R/site-library')"

RUN Rscript -e 'if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager",lib="/usr/lib/R/site-library")' && R -e 'BiocManager::install(c("rtracklayer","GenomicRanges","BSgenome.Hsapiens.UCSC.hg38","ShortRead"),lib="/usr/lib/R/site-library")'

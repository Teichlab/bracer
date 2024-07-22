# run as:
#   docker build -t bracer .

#start off with a plain Debian
FROM debian:trixie-20230703-slim

#basic setup stuff, including bowtie2 for Bracer, libcurl4-openssl-dev r-base libxml2-dev for Alakazam
RUN apt-get update && apt-get -y upgrade && \
apt-get -y install \
bowtie2 \
build-essential \
cmake \
curl \
default-jre \
git \
graphviz \
jellyfish \
libcairo2-dev \
libcurl4-openssl-dev \
libfreetype6-dev \
libgirepository1.0-dev \
libxml2-dev \
pkg-config \
python3-dev \
python3-pip \
python3-venv \
r-base \
salmon \
samtools  \
unzip \
wget \
zlib1g-dev \
&& \
apt-get clean && \
rm -rf /var/lib/apt/lists/*

#Trinity - depends on zlib1g-dev and openjdk-8-jre installed previously
RUN wget https://github.com/trinityrnaseq/trinityrnaseq/releases/download/Trinity-v2.15.1/trinityrnaseq-v2.15.1.FULL.tar.gz 
RUN tar -xf trinityrnaseq-v2.15.1.FULL.tar.gz
# Currently trinity won't compile without a #include <string> in this file
RUN sed -i '1s;^;#include <string>\n;' /trinityrnaseq-v2.15.1/trinity-plugins/bamsifter/sift_bam_max_cov.cpp
RUN cd /trinityrnaseq-v2.15.1 && make

#IgBLAST
RUN wget https://ftp.ncbi.nih.gov/blast/executables/igblast/release/1.21.0/ncbi-igblast-1.21.0-x64-linux.tar.gz && \
tar -xf ncbi-igblast-1.21.0-x64-linux.tar.gz && \
rm ncbi-igblast-1.21.0-x64-linux.tar.gz


COPY docker_helper_files/gencode_parse.py /bracer/docker_helper_files/gencode_parse.py

#aligners - kallisto
RUN wget https://github.com/pachterlab/kallisto/releases/download/v0.48.0/kallisto_linux-v0.48.0.tar.gz && tar -xzvf kallisto_linux-v0.48.0.tar.gz && rm kallisto_linux-v0.48.0.tar.gz

#obtaining the transcript sequences, no need for kallisto/salmon indices
RUN mkdir GRCh38 && cd GRCh38 && wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.transcripts.fa.gz && \
gunzip gencode.v43.transcripts.fa.gz && python3 /bracer/docker_helper_files/gencode_parse.py gencode.v43.transcripts.fa && rm gencode.v43.transcripts.fa

RUN mkdir GRCm38 && cd GRCm38 && wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M32/gencode.vM32.transcripts.fa.gz && \
gunzip gencode.vM32.transcripts.fa.gz && python3 /bracer/docker_helper_files/gencode_parse.py gencode.vM32.transcripts.fa && rm gencode.vM32.transcripts.fa

#regular BLAST
RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.14.0/ncbi-blast-2.14.0+-x64-linux.tar.gz \
&&  tar -xzvf ncbi-blast-2.14.0+-x64-linux.tar.gz && rm ncbi-blast-2.14.0+-x64-linux.tar.gz

#phylip
RUN wget https://phylipweb.github.io/phylip/download/phylip-3.697.tar.gz && tar -xzvf phylip-3.697.tar.gz && rm phylip-3.697.tar.gz
RUN cd phylip-3.697/src && sed -i 's/^CFLAGS =/CFLAGS = -fcommon/g' Makefile.unx && make -f Makefile.unx install

#Trim Galore! plus its dependency FastqC
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip && unzip fastqc_v0.12.1.zip && rm fastqc_v0.12.1.zip
RUN chmod 755 /FastQC/fastqc
RUN ln -s /FastQC/fastqc /usr/local/bin/fastqc
RUN curl -fsSL  https://github.com/FelixKrueger/TrimGalore/archive/refs/tags/0.6.10.tar.gz -o trim_galore.tar.gz
RUN tar xvzf trim_galore.tar.gz && mv TrimGalore-0.6.10/trim_galore /usr/bin

#R dependencies
RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install(c('GenomicAlignments', 'Biostrings', 'IRanges'))"
RUN R -e "install.packages(c('remotes', 'ggplot2'), repos='http://cran.us.r-project.org')"
RUN R -e "require(remotes);install_version('alakazam', version = '1.2.0', repos='http://cran.us.r-project.org')"


#Bowtie 2 as needs version 2.5.1 due to a bug in 2.5.0
RUN wget https://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.5.1/bowtie2-2.5.1-linux-x86_64.zip && unzip bowtie2-2.5.1-linux-x86_64.zip && rm bowtie2-2.5.1-linux-x86_64.zip

#bracer proper, no need to reposition resources as config will now know where this lives
COPY ./docker_helper_files/requirements_stable.txt ./bracer/docker_helper_files/
WORKDIR /bracer
RUN python3 -m venv venv
ENV VIRTUAL_ENV=/bracer/venv
ENV PATH=$VIRTUAL_ENV/bin:$PATH
RUN pip install --upgrade pip --break-system-packages
RUN pip3 install -r docker_helper_files/requirements_stable.txt --break-system-packages
COPY . /bracer
RUN python3 setup.py install
WORKDIR /

#placing a preconfigured bracer.conf in ~/.bracerrc
RUN cp /bracer/docker_helper_files/docker_bracer.conf /home/.bracerrc

ENV BRACER_CONF=/bracer/docker_helper_files/docker_bracer.conf
ENV IGDATA=/ncbi-igblast-1.21.0/

#this is a bracer container, so let's point it at bracer and set -e
ENTRYPOINT ["bash", "/bracer/docker_helper_files/docker_wrapper.sh"]

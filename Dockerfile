# run as:
#   docker build -t bracer .

#start off with a plain Debian
FROM debian:latest

#basic setup stuff, including bowtie2
RUN apt-get update && apt-get -y upgrade
RUN apt-get -y install wget curl unzip build-essential zlib1g-dev git python3 python3-pip bowtie2 openjdk-8-jre

#Trinity - depends on zlib1g-dev and openjdk-8-jre installed previously
RUN wget https://github.com/trinityrnaseq/trinityrnaseq/archive/Trinity-v2.4.0.zip
RUN unzip Trinity-v2.4.0.zip && rm Trinity-v2.4.0.zip
RUN cd /trinityrnaseq-Trinity-v2.4.0 && make

#IgBLAST, plus the setup of its super weird internal_data thing. don't ask. just needs to happen
#and then on top of that, the environmental variable thing facilitates the creation of a shell wrapper. fun
RUN wget ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/1.4.0/ncbi-igblast-1.4.0-x64-linux.tar.gz
RUN tar -xzvf ncbi-igblast-1.4.0-x64-linux.tar.gz && rm ncbi-igblast-1.4.0-x64-linux.tar.gz
RUN cd /ncbi-igblast-1.4.0/bin/ && wget -r ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/internal_data && \
	mv ftp.ncbi.nih.gov/blast/executables/igblast/release/internal_data . && rm -r ftp.ncbi.nih.gov

#aligners - kallisto and salmon
RUN wget https://github.com/pachterlab/kallisto/releases/download/v0.43.1/kallisto_linux-v0.43.1.tar.gz
RUN tar -xzvf kallisto_linux-v0.43.1.tar.gz && rm kallisto_linux-v0.43.1.tar.gz
RUN wget https://github.com/COMBINE-lab/salmon/releases/download/v0.8.2/Salmon-0.8.2_linux_x86_64.tar.gz
RUN tar -xzvf Salmon-0.8.2_linux_x86_64.tar.gz && rm Salmon-0.8.2_linux_x86_64.tar.gz

#graphviz, along with its sea of dependencies that otherwise trip up the dpkg -i
RUN apt-get -y install libgd3 libgts-0.7-5 liblasi0 libltdl7 freeglut3 libglade2-0 libglu1-mesa libglu1 libgtkglext1 libxaw7
RUN wget http://www.graphviz.org/pub/graphviz/stable/ubuntu/ub13.10/x86_64/libgraphviz4_2.38.0-1~saucy_amd64.deb
RUN dpkg -i libgraphviz4_2.38.0-1~saucy_amd64.deb && apt-get -y -f install
RUN wget http://www.graphviz.org/pub/graphviz/stable/ubuntu/ub13.10/x86_64/graphviz_2.38.0-1~saucy_amd64.deb
RUN dpkg -i graphviz_2.38.0-1~saucy_amd64.deb && apt-get -y -f install
RUN rm libgraphviz4_2.38.0-1~saucy_amd64.deb && rm graphviz_2.38.0-1~saucy_amd64.deb

#regular BLAST
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.6.0/ncbi-blast-2.6.0+-x64-linux.tar.gz
RUN tar -xzvf ncbi-blast-2.6.0+-x64-linux.tar.gz && rm ncbi-blast-2.6.0+-x64-linux.tar.gz

#phylip
RUN wget http://evolution.gs.washington.edu/phylip/download/phylip-3.696.tar.gz
RUN tar -xzvf phylip-3.696.tar.gz && rm phylip-3.696.tar.gz
RUN cd phylip-3.696/src && make -f Makefile.unx install

#Trim Galore! plus its dependency FastqC
RUN wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip
RUN unzip fastqc_v0.11.5.zip && rm fastqc_v0.11.5.zip
RUN chmod 755 /FastQC/fastqc
RUN ln -s /FastQC/fastqc /usr/local/bin/fastqc
RUN curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.4.3.tar.gz -o trim_galore.tar.gz
RUN tar xvzf trim_galore.tar.gz && mv TrimGalore-0.4.3/trim_galore /usr/bin

#R dependencies. libxml2-dev is a ghost dependency of an alakazam dependency not mentioned by the install crash
RUN apt-get -y install r-base libxml2-dev
RUN R -e "install.packages(c('alakazam', 'ggplot2'), repos='http://cran.us.r-project.org')"

#bracer proper, no need to reposition resources as config will now know where this lives
COPY . /bracer
RUN cd /bracer && pip3 install -r docker_helper_files/requirements_stable.txt && python3 setup.py install

#obtaining the transcript sequences, no need for kallisto/salmon indices
RUN mkdir GRCh38 && cd GRCh38 && wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_27/gencode.v27.transcripts.fa.gz && \
	gunzip gencode.v27.transcripts.fa.gz && python3 /bracer/docker_helper_files/gencode_parse.py gencode.v27.transcripts.fa && rm gencode.v27.transcripts.fa
RUN mkdir GRCm38 && cd GRCm38 && wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M15/gencode.vM15.transcripts.fa.gz && \
	gunzip gencode.vM15.transcripts.fa.gz && python3 /bracer/docker_helper_files/gencode_parse.py gencode.vM15.transcripts.fa && rm gencode.vM15.transcripts.fa

#placing a preconfigured bracer.conf in ~/.bracerrc
RUN cp /bracer/docker_helper_files/docker_bracer.conf ~/.bracerrc

#this is a bracer container, so let's point it at a bracer wrapper that sets the silly IgBLAST environment variable thing
ENTRYPOINT ["bash", "/bracer/docker_helper_files/docker_wrapper.sh"]

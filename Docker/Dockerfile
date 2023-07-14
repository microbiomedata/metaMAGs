FROM continuumio/miniconda3:4.12.0
## python 3.8.5
MAINTAINER chienchi@lanl.gov

LABEL version="0.5"
LABEL software="NMDC_MAG_Pipeline"
LABEL tags="bioinformatics"

# Tell systemd we're running in Docker!

ENV container docker

ADD * /opt/conda/bin/
RUN cp opt/conda/lib/libgomp.so.1 /lib

RUN \
    conda install -c bioconda -y biopython=1.78 metabat2=2.15 hmmer=3.3.2 prodigal=2.6.3 pplacer=1.1.alpha19 fasttree=2.1.11 samtools=1.6 numpy==1.23.5 && \
    conda install -y sqlite=3.39.2 && \
    pip3 install matplotlib==3.7.1 pysam==0.20.0 checkm-genome==1.2.1 gtdbtk==2.1.1 
   # conda install -c bioconda/label/cf201901 -y metabat2 && \

ENV GTDBTK_DATA_PATH /databases/GTDBTK_DB
ENV CHECKM_DATA_PATH /databases/CheckM_DB

WORKDIR /data

CMD ["/bin/bash"]


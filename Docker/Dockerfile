FROM continuumio/miniconda3:4.9.2
## python 3.8.5
MAINTAINER chienchi@lanl.gov

LABEL version="0.1.6"
LABEL software="NMDC_MAG_Pipeline"
LABEL tags="bioinformatics"

# Tell systemd we're running in Docker!

ENV container docker

RUN apt update --allow-releaseinfo-change -y | true && apt install --no-install-recommends -y vim procps gzip zip time libgomp1 samtools && \
    conda config --add channels defaults && \ 
    conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    mkdir -p /data && \
    mkdir -p /databases/GTDBTK_DB && \ 
    mkdir -p /checkM_DB && cd /checkM_DB && \
    wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz && \
    tar xvzf checkm_data_2015_01_16.tar.gz && \
    rm checkm_data_2015_01_16.tar.gz && \
    conda install -c bioconda -y \
       biopython=1.78 \
       boost=1.70 \
       boost-cpp=1.70 \
       fastani=1.32 \
       fasttree=2.1.10 \
       hmmer=3.1b2 \
       sqlite=3.33.0 \
       pymysql=1.0.2 \
       requests=2.27.1 \
       checkm-genome=1.1.9 \
       gtdbtk=1.3.0 && \
    checkm data setRoot /checkM_DB  2>/dev/null && \
    apt autoremove -y && apt clean && \
    conda clean -y -a && \
    cd /root && rm -rf .cache /tmp/* /var/lib/apt/lists/* /var/tmp/* 
   # conda install -c bioconda/label/cf201901 -y metabat2 && \

ADD *.py /opt/conda/bin/
ADD metabat2_bin/* /opt/conda/bin/
ADD init.sh /

ENV GTDBTK_DATA_PATH /databases/GTDBTK_DB

WORKDIR /data

ENTRYPOINT ["/init.sh"]

CMD ["/bin/bash"]


FROM continuumio/miniconda3:latest
MAINTAINER chienchi@lanl.gov

LABEL version="0.1.2"
LABEL software="NMDC_MAG_Pipeline"
LABEL tags="bioinformatics"

# Tell systemd we're running in Docker!

ENV container docker

ADD *.py /opt/conda/bin/
ADD metabat2_bin/* /opt/conda/bin/
ADD init.sh /
RUN apt update && apt install -y vim procps gzip time libgomp1 && \
    conda config --add channels defaults && \ 
    conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    mkdir -p /data && \
    mkdir -p /databases/GTDBTK_DB && \ 
    mkdir -p /checkM_DB && cd /checkM_DB && \
    wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz && \
    tar xvzf checkm_data_2015_01_16.tar.gz && \
    rm checkm_data_2015_01_16.tar.gz && \
    conda install -c bioconda -y samtools biopython fastani fasttree sqlite  pymysql requests checkm-genome gtdbtk=1.2.0 && \
    checkm data setRoot /checkM_DB  2>/dev/null && \
    apt autoremove -y && apt clean && \
    conda clean -y -a && \
    cd /root && rm -rf .cache /tmp/* /var/lib/apt/lists/* /var/tmp/* 
   # conda install -c bioconda/label/cf201901 -y metabat2 && \

ENV GTDBTK_DATA_PATH /databases/GTDBTK_DB

WORKDIR /data

ENTRYPOINT ["/init.sh"]

CMD ["/bin/bash"]


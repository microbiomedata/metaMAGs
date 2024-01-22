FROM mambaorg/micromamba:latest
## python 3.9.12
LABEL version="0.5"
LABEL software="NMDC_MAG_Pipeline"
LABEL tags="bioinformatics"
LABEL maintainer="Chienchi Lo"

# Tell systemd we're running in Docker!

ENV container docker
ENV DEBIAN_FRONTEND=noninteractive
ENV LANG=en_US.UTF-8
ENV LC_ALL=C

USER root

RUN \
    micromamba create -n mags -c conda-forge -c bioconda -c dnachun -c defaults -y python=3.9.12 biopython=1.78 metabat2=2.15 hmmer=3.3.2 prodigal=2.6.3 pplacer=1.1.alpha19 fasttree=2.1.11 mash=2.3 samtools=1.6 numpy=1.23.5 sqlite=3.39.2 google-cloud-sdk gcsfuse && \
    eval "$(micromamba shell hook --shell bash)" && \ 
    micromamba activate mags && \
    pip3 install matplotlib==3.7.1 pysam==0.20.0 checkm-genome==1.2.1 gtdbtk==2.1.1  && \
    micromamba clean -y -a && pip cache purge

ADD * /opt/conda/envs/mags/bin/
RUN cp /opt/conda/envs/mags/lib/libgomp.so.1 /lib

ENV PATH="${PATH}:/opt/conda/envs/mags/bin"
ENV GTDBTK_DATA_PATH /databases/GTDBTK_DB
ENV CHECKM_DATA_PATH /databases/CheckM_DB

WORKDIR /data

CMD ["/bin/bash"]


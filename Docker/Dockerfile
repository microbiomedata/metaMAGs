FROM doejgi/mbin:v1.0.1
## python 3.9.12

LABEL maintainer="chienchi@lanl.gov"
LABEL version="1.0.5"
LABEL software="NMDC_MAG_Pipeline"
LABEL tags="bioinformatics"

# Tell systemd we're running in Docker!

ENV container="docker"


ADD Docker/* /opt/conda/bin/


ENV GTDBTK_DATA_PATH=/databases/GTDBTK_DB
ENV CHECKM_DATA_PATH=/databases/CheckM_DB
ENV EUKCC2_DB=/databases/EUKCC2_DB

WORKDIR /data

CMD ["/bin/bash"]


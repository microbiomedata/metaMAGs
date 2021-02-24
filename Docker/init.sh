#!/usr/bin/env bash

checkm data setRoot /databases/checkM_DB  2>/dev/null 
chmod 777 /opt/conda/lib/python3.7/site-packages/checkm/DATA_CONFIG
chmod -R 777 /databases/checkM_DB 

export GTDBTK_DATA_PATH=/databases/GTDBTK_DB

exec "$@"



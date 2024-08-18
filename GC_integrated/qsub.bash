#!/bin/bash
#! -S /bin/bash

module use /usr/local/package/modulefiles/
module load R/4.4.0
module load python/3.12.0

OUTPUT_DIR="Rout"

RFILE_NAME1="14DEG_target_gene"
RFILE_NAME2=""
RFILE_NAME3=""
RFILE_NAME4=""
RFILE_NAME5=""

if [[ -n "$RFILE_NAME1" ]]; then
    R CMD BATCH --no-save --no-restore src/${RFILE_NAME1}.R ${OUTPUT_DIR}/${RFILE_NAME1}.Rout
fi

if [[ -n "$RFILE_NAME2" ]]; then
    R CMD BATCH --no-save --no-restore src/${RFILE_NAME2}.R ${OUTPUT_DIR}/${RFILE_NAME2}.Rout
fi

if [[ -n "$RFILE_NAME3" ]]; then
    R CMD BATCH --no-save --no-restore src/${RFILE_NAME3}.R ${OUTPUT_DIR}/${RFILE_NAME3}.Rout
fi

if [[ -n "$RFILE_NAME4" ]]; then
    R CMD BATCH --no-save --no-restore src/${RFILE_NAME4}.R ${OUTPUT_DIR}/${RFILE_NAME4}.Rout
fi

if [[ -n "$RFILE_NAME5" ]]; then
    R CMD BATCH --no-save --no-restore src/${RFILE_NAME5}.R ${OUTPUT_DIR}/${RFILE_NAME5}.Rout
fi
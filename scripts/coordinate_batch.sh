#!/bin/bash

$FILE_TO_SUBMIT=$1
$AFTER_JOB="after_job.sh"

jid1=$(sbatch --mem=12g --cpus-per-task=4 $FILE_TO_SUBMIT)
jid2=$(sbatch --dependency=afterany:$jid1 $AFTER_JOB)
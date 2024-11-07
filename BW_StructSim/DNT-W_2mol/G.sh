#!/bin/sh

#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -pe gau 12
#$ -q all.q

module load gaussian/g16
export GAUSS_SCRDIR=/scr/$JOB_ID
mkdir /scr/$JOB_ID

g16 xxx.gjf
rm -rf /scr/$JOB_ID


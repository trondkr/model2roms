#!/bin/bash
#
#  Give the job a name
#PBS -N "model2roms_WS4KM"
#
#  Specify the project the job belongs to
#PBS -A nn9297k
#PBS -q normal
#PBS -l mppwidth=1,walltime=06:10:00
#PBS -l mppmem=1000MB

#
#  Send me an email on  a=abort, b=begin, e=end
#PBS -m abe
#
#  Use this email address (check that it is correct):
#PBS -M trond.kristiansen@niva.no
#
#  Write the standard output of the job to file 'mpijob.out' (optional)
#PBS -o  model2roms_A20.out
#
#  Write the standard error of the job to file 'mpijob.err' (optional)
#PBS -e  model2roms_A20.err
#

#  Make sure I am in the correct directory
cd /work/shared/nn9297k/model2roms
module load python

export MPLCONFIGDIR=${pwd}
export TMP=`pwd`
export PYTHON_EGG_CACHE=/work/shared/nn9297k/model2roms

aprun -B python main.py > output.log
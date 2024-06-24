#!/bin/bash
#JSUB -n 35
##JSUB -app g16
##JSUB -m gpu12
#JSUB -q gpu
##JSUB -gpgpu 1
#JSUB -e error.%J
#JSUB -o output.%J
#JSUB -J y61500
export g16root=~/YWQ
export GAUSS_SCRDIR=$g16root/g16/tmp
source $g16root/g16/bsd/g16.profile
export PATH=$PATH:$g16root/g16

icc=0
for i in {1..200}
do
for inf in ${i}.gjf
do
((icc++))
echo RUNNING $icc

time g16 < ${inf} > ${inf//gjf/log}
rename fort.7 ${inf//gjf/pun} fort.7
echo $icc has finished

done
done

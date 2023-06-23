#!/bin/bash

usage="$0 input window\n"

if [[ $# -ne 2 ]] || [[ $1 == "-h" ]] || [[ $1 == "--help" ]]
then
    echo -e $usage
    exit
fi

awk 'BEGIN{FS=OS="\t";i=n=total=0}{if(i<window-1){total+=$2;i+=1}else{print (n,total/window);n+=1;total=i=0}}' window=$2 $1 > ${1}_sum_$2.csv

#!/bin/bash -l

echo "../divine.py -q ./Pfeiffer.hpo -v ./Pfeiffer.vcf -o ./Pfeiffer -e 1 -c ../../../config/filterconf_dp10.txt -k 0 --reuse"

../divine.py -q ./Pfeiffer.hpo -v ./Pfeiffer.vcf -o ./Pfeiffer -e 1 -c ../../../config/filterconf_dp10.txt -k 0 --reuse

#!/bin/bash -l

echo "../divine.py -v ./Pfeiffer.vcf -o ./Pfeiffer_noHpo -e 1 -c ../../../config/filterconf_dp10.txt --reuse -k 0"

../divine.py -v ./Pfeiffer.vcf -o ./Pfeiffer_noHpo -e 1 -c ../../../config/filterconf_dp10.txt --reuse -k 0

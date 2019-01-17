#!/bin/bash -l

echo "../divine.py -q ./Pfeiffer.hpo -v ./Pfeiffer.vcf -o ./Pfeiffer -e 1 -c ../../../config/divine.conf -k 1 --reuse"

../divine.py -q ./Pfeiffer.hpo -v ./Pfeiffer.vcf -o ./Pfeiffer -e 1 -c ../../../config/divine.conf -k 1 --reuse


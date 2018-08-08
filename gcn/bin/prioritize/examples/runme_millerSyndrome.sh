#!/bin/bash -l

echo "../divine.py -q ./MillerSyndrome.hpo -v ./MillerSyndrome.vcf -k 0 --reuse -o ./MillerSyndrome"

../divine.py -q ./MillerSyndrome.hpo -v ./MillerSyndrome.vcf -k 0 --reuse -o ./MillerSyndrome

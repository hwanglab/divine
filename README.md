![Divine_workflow](documents/oxfj_bioinf2018_hong/slide_figures/divine_logo_small.PNG)

# *Divine* : Prioritizing Genes for Rare Mendelian Disease in Whole Exome Sequencing Data

Divine is designed to make a daily-routine molecular diagnosis with high-throughput whole exome sequencing data more efficient. Using both patient phenotypic information and genetic variants, Divine that integrates patients’ phenotype(s) and WES data with 30 prior biological knowledge (e.g., human phenotype ontology, gene ontology, pathway database, protein-protein interaction networks, etc.) to prioritize potential disease-causing genes.

![Divine_workflow](documents/oxfj_bioinf2018_hong/slide_figures/divine_v0.7.png)

## Website

https://github.com/hwanglab/divine
# Tutorial
A [tutorial](https://github.com/hwanglab/divine/blob/master/documents/tutorial/divine_tutorial.md) is available, from the installation to case studies.

# Input
- a text file containing Human Phenotype Ontology (HPO) IDs (e.g., HP:0002307) that describe patient's clinical features

- or (and) VCF file

# Algorithm
- Refer to [Divine: Prioritizing Genes for Rare Mendelian Disease in Whole Exome Sequencing Data", Hong. *et. al.*](https://www.biorxiv.org/content/early/2018/08/21/396655)

# Output
- When only HPO IDs are given, Divine generates a prioritized gene list and an inferred disease list.

- If VCF file (or with HPO IDs), Divine also generates an annotated variant table with a ranking score in Microsoft Excel format.

# Developer note
- Divine is in active development. Please check the above website for newer available versions. Please contact the developer at hongc2{at}ccf.org to report any problems or for additional help.

# Setup

## Prerequisite:
- Linux, 4 GB RAM, and 130 GB available hard disk space
- python 2.7+
- make sure that 'pip', 'wget' (with Internet connection), 'grep', 'awk', and 'sort' command in $PATH
- locate your shell configuration file (e.g., $HOME/.bash_profile, $HOME/.profile, or $HOME/.cshrc)
- It does not require a super user (root) account.

## Python modules to be installed
Divine requires the following modules but, during the setup process, the modules will be installed automatically if necessary.

- *`fastsemsim-0.9.4`: https://sites.google.com/site/fastsemsim
- *`hgvs`: https://github.com/counsyl/hgvs
- *`hpo_similarity`: https://github.com/jeremymcrae/hpo_similarity
- `ConfigParser, backports-abc, html5lib==0.999999999, backports.ssl-match-hostname, certifi,decorator, matplotlib, networkx, nose, numpy, pandas, pygr, pyparsing, pysam, python-dateutil, pytz,scipy, scikit-learn, singledispatch, six, tornado, xlwt, dill, ped_parser`

- [*]: the python modules are already included in the Divine package.

## Install

download divine source codes from github
```
$ git clone https://github.com/hwanglab/divine.git
```

### Option 1: fresh install

It requires downloading 21 GB database files and so be patient!
```
$ setup.py --install
```

### Optional 2: upgrade Divine database/resource/examples/3rd-party python modules (when you previously installed Divine)
```
$ setup.py --install --update_db
```

## Configuration
- Get updated environment variables (PATH, DIVINE, and PYTHONPATH) at the end of the installation setup message. Then,

- Add the environment variables into your shell configuration (e.g., $HOME/.bash_profile, $HOME/.profile, or $HOME/.cshrc)

- Make sure that the new configuration is applied by
  - logoff and login your account, or
  - running
  ```
  $ source $HOME/your_shell_config_file
  ```

## Uninstall
### Optional 1: only uninstall python modules of dependency
```
$ setup.py --uninstall
```

### Optional 2: also, uninstall resource files
```
$ setup.py --uninstall --remove_db
```

# Usage
## Input: a text file containing HPO IDs

First, visit either [http://compbio.charite.de/phenomizer](http://compbio.charite.de/phenomizer), [https://hpo.jax.org/app](https://hpo.jax.org/app), or [https://mseqdr.org/search_phenotype.php](https://mseqdr.org/search_phenotype.php). Enter the patient phenotype terms, description, or specific keywords you think important. Get the best matching HPO IDs. Paste the HPO IDs in the format of HP:XXXXXXX (e.g., HP:0002307) into a text file line by line and save it as a text file (e.g., P0001.hpo)

For example, an HPO file looks like
```
$ cat P0001.hpo
#my_patient_ID
HP:0002307
HP:0000639
HP:0001252
HP:0100543
HP:0002120
HP:0000708
HP:0001344
HP:0008872
HP:0000510
HP:0001513
HP:0006979
HP:0000752
HP:0012469
HP:0000577
HP:0001010
HP:0006887
HP:0002650
HP:0005469
HP:0002312
HP:0010808
HP:0002136
HP:0200085
HP:0002311
```

## Input: a VCF file
Any VCF file following the standard format (e.g., [https://samtools.github.io/hts-specs/VCFv4.2.pdf](https://samtools.github.io/hts-specs/VCFv4.2.pdf)). Refer to by GATK, samtools, or, freebayes to convert FASTQ to VCF file.

## Output

### Known disease matching by patient phenotypes
`hpo_to_diseases.tsv`: From an input HPO file, Divine prioritize which disease the patient likely has. The output format is
```
$ head -n3 hpo_to_disease.tsv
#disease_ID    genes    score[FunSimMax]
OMIM:101600    FGFR1,FGFR2    0.000911782
OMIM:101200    FGFR2    0.000674322
:
:
```
- disease_ID: OMIM, ORPHANET, or DECIPHER disease ID
- genes: gene symbol
- score[FunSimMax]: patient phenotype matching score to the disease

### Gene rank
- gene: gene symbol
- predicted_score: finalized score after heat diffusion
- seed_score: a score combined by phenotype matching score and genetic damage score
- gt_dmg_score: a genetic pathogenicity score assessed in VCF file
- pheno_score: score[FunSimMax] 
 
```
head -n6 gene_rank.tsv
#gene   predicted_score seed_score      gt_dmg_score    pheno_score     contain_known_pathogenic
FGFR2   0.00925436      1.60806e-06     0.00545116      0.00117607      NA
EVC     0.00634702      1.08469e-06     0.00981397      0.000439834     NA
VPS13B  0.00477042      8.25496e-07     0.00322986      0.00102016      NA
CCT4    0.0045629       7.92439e-07     0.00648582      0.000487018     NA
ZMYND11|ENSP00000452959 0.00405376      6.66455e-07     0.00745641      0.000356124     NA
```
 
### Gene enrichment from known diseases matched by patient phenotypes
```bash
head -n6 diseases_rank.tsv
```
- disease_ID: disease ID
- disease_description: disease description
- inheritance: 
        - inheritance model, autosomal recessive or dominant assoc_pheno_genes (^:mutated, *:known_pathogenic)known to be associated with the disease
        - ^: Does this gene have a mutated variant in the patient VCF file?
        - *: if yes, then variant is reported to a known pathogenic in either ClinVar or HGMD?
- num_of_assoc_pheno_genes: the number of genes associated with the disease
- num_of_gt_dmg_genes: the number of genes containing variants
- pheno_match_score: A phenotype matching score to the disease
- avg_combined_dmg_score: average combined pathogenic score
- max_combined_dmg_score: maximum combined pathogenic score
- avg_harmonic_score: average finalized score yi
- max_harmonic_score: max finalized score yi
- external_genes_of_interest(kegg-ppi_or_GO_enriched[harmonic_score]): mutated gene enriched by KEGG pathway, PPI, or Gene Ontology database
- PPI-KEGG_pathway_desc: PPI and KEGG pathway description


### Annotated VCF files
- `divine.vcf`
    A VCF file annotated by Varant where we improve the original method significantly. Refer to Varant website (http://compbio.berkeley.edu/proj/varant).

- `vfilter.vcf`
    A subset of the VCF file above (i.e., `divine.vcf`). This VCF file contains variants after filtering out all variant calls that are either located in an intergenic region or its MAF is high or in repeated region.

### Prioritized genes and Microsoft Excel file
- `divine.xls`
   An excel table generated by Varant. Divine annotates a ranking score per gene so that you can sort the variant list by the ranking score (refer to 'Comments_on_Genes', an interpreted column in the 2nd tab, VCF).

### Log files

- check log files under a subdirectory `logs` in an output directory

## Cases

- When HPO file is only given,

```
../divine.py -q dir_to_the_hpo/P0001.hpo -o dir_to_output/P0001
```

- When VCF file is only given,

```
../divine.py -v dir_to_the_vcf/P0002.vcf -o dir_to_output/P0002
```

- When both HPO file and VCF file are given,

```
../divine.py -q dir_to_the_hpo/P0003.hpo -v dir_to_the_vcf/P0003.vcf -o dir_to_output/P0003
```

- When both HPO file and VCF file are given,

```
../divine.py -q dir_to_the_hpo/P0003.hpo -v dir_to_the_vcf/P0003.vcf -o dir_to_output/P0003
```

## Examples
We include 4 to 5 demo samples in the resource package, 
	
	```
	cd $DIVINE/gcn/bin/prioritize/examples
	./runme_angelman.sh #when only HPO data is available
	./runme_pfeisffer_noHpo.sh #when only VCF is available 
	./runme_pfeisffer.sh #when both HPO and VCF are available
	./runme_millerSyndrome.sh #when both HPO and VCF are available
	./runme_trio.sh #analyze family samples (PED file requires and sample ID should be matched with the ones in VCF file)
	```

# Help
```
usage: divine.py [-h] [-q HPO_QUERY] [-v VCF] [-o OUT_DIR] [-c VCF_FILTER_CFG]
                 [-f PED] [-p PROBAND_ID] [-d EXP_TAG] [-i INDEL_FIDEL]
                 [-K TOP_K_DISEASE] [-r GO_SEED_K] [-e REF_EXON_ONLY]
                 [-C CADD] [-j COSMIC] [-D DBLINK] [-H HGMD] [-k VKNOWN]
                 [-t CAPKIT] [--reuse]

Divine (v0.1.2) [author:hongc2@ccf.org]

optional arguments:
  -h, --help            show this help message and exit
  -q HPO_QUERY, --hpo HPO_QUERY
                        Input patient HPO file. A file contains HPO IDs (e.g.,
                        HP:0002307), one entry per line. Refer to
                        http://compbio.charite.de/phenomizer or
                        https://mseqdr.org/search_phenotype.php
  -v VCF, --vcf VCF     input vcf file
  -o OUT_DIR, --out_dir OUT_DIR
                        output directory without white space. If not exist,
                        the directory will be created.
  -c VCF_FILTER_CFG, --vcf_filter_cfg VCF_FILTER_CFG
                        vcf filter configuration file [None]
  -f PED, --family_fn PED
                        family pedigree file [None]
  -p PROBAND_ID, --proband_id PROBAND_ID
                        proband sample ID [None]
  -d EXP_TAG, --exp_tag EXP_TAG
                        specify experiment tag without white space. The tag
                        will be contained in the output file name.[None]
  -i INDEL_FIDEL, --indel INDEL_FIDEL
                        the level of fidelity of indell call in VCF, [1]:low
                        (e.g., samtools), 2:high (GATK haplotype caller)
  -K TOP_K_DISEASE      focus on top-K disease associated with the input
                        phenotypes [0], set 0 to consider all
  -r GO_SEED_K, --go_seed_k GO_SEED_K
                        the number of top-k diseases for GO enrichment [3];
                        set to 0 to disable
  -e REF_EXON_ONLY, --ref_exon_only REF_EXON_ONLY
                        the annotation process only runs on RefSeq coding
                        regions 0:No, [1]:Yes
  -C CADD, --cadd CADD  use CADD prediction score, 0:No, [1]:Yes
  -j COSMIC, --cosmic COSMIC
                        enable COSMIC, [0]:No, 1:Yes
  -D DBLINK, --dblink DBLINK
                        enable dblink, [0]:No, 1:Yes
  -H HGMD, --hgmd HGMD  enable HGMD (requires a license), [0]:No, 1:Yes
  -k VKNOWN, --vknown VKNOWN
                        apply variant-level pathogenic annotation (e.g.,
                        either ClinVar or HGMD) to prioritization strategy,
                        0:No, [1]:Yes
  -t CAPKIT             capture kit symbol [None],SureSelect_V6,SeqCapEZ_Exome
  --reuse               Reuse previous annotation file (divine.vcf) if it is
                        available [False]
```

# FAQ

- Q.1: I have a VCF file generated from either WGS or WES dataset, and the number of variants is a lot! It seems that Divine is slow and I want to have a result as soon as possible.

- A.1: Depending on your hardware specification and the number of variants in an input VCF file, the computational time varies. In my computer setting (e.g., Intel Core 2 Duo @ 2.93GHz), it takes 30 min to handle 37,000 variants. The computational bottleneck is directly associated with hard disk drive I/O. HDD connected to USB 3.0 or SDD (solid-state drive) can be helpful. Also, we are actively working on Divine to speed up the process. For now, Divine is optimized to analyze WES or targeted panel samples but not WGS. As a default option, Divine focuses on NCBI RefGene exon regions with +/-20bp flanking. However, be aware that it may not detect a pathogenic variant in an intergenic region or up/downstream. For a VCF file containing a small number of variants, we suggest `-e 0` instead. 

- Q.2: I don't like a default filtering scheme used in Divine. I want my own filtering strategy (e.g., include a specific flag in FILTER; not to use ExAC; 0.03 for MAF cutoff).

- A.2: open $DIVINE/gcn/config/filterconf.txt and edit the configuration file. For example, set `excl=LowQual` to filter out all LowQual in the VCF file. Specify the option with `-c` in Divine option. Later, we will provide more comprehensive instruction for this file. The default configuration is

```
[fltr]
excl=LowDP:LowGQ:LowGQX:LowQual:IndelGap:SnpGap
#incl=PASS
[infoflag]
excl=DB137
[infoval]
kgaf=yes
espaf=yes
exacaf=yes
splice_dist=20
hgmd_filter=2
regulome=no
#min_depth=10
[reg]
incl=CodingExonic:NonCodingExonic:CodingIntronic:NonCodingIntronic
[freq]
incl=0.01
[freq_cli]
incl=0.05
[gid]
min=0.1
```

- Q.3: I purchase HGMD professional license and how can I use this feature?

- A.3: Contact to hongc2@ccf.org

- Q.4: Previously, I ran Divine on a patient sample dataset which took so long. Now, I want to prioritize genes with a different setting (e.g., filtering condition or a different option in `divine.py`) and I am sure that I didn't update Divine/GCN database. How can I make the repeat analysis faster?

- A.4: Unless Divine database is changed, or purchase HGMD license newly, or divine.vcf created previously is truncated, you can reuse the previous annotated VCF file using `--reuse`. If you want to compare a new result with the previous one, try

```
$ $DIVINE/gcn/bin/prioritize/divine.py -q dir_to_the_hpo/P0005.hpo \
   -v dir_to_the_vcf/P0005.vcf --reuse -d rev -o dir_to_the_output/P0005
```

- Q.5: For the same input and configuration, will Divine generate the same output?

- A.5: Yes, Divine is deterministic. There is no randomness in the analysis. Also, Divine generates all log files so that you can audit/trace the previous experimental results and database maintenance at any time.

- Q.6: I think database or annotation Divine provides is outdated. I would like to keep all the latest database so that I can improve my diagnosis result at best.

- A.6: Divine uses Varant as an annotation framework. We will frequently update the database. Each resource package will be archived. Furthermore, we plan to provide a stream pipeline to update annotation database.

- Q.7: I know Divine is designed for germline or constitutional disease samples. Can I use Divine for somatic mutation (cancer sample) analysis?

- A.7: We will work on this feature as well.

- Q.8: I have a patient dataset but I concern if Divine collects the patient dataset without my permission and send it somewhere.

- A.8: This standalone package never collects any user information or input dataset the user works. But, we need your feedback and bug report!

- Q.9: Does Divine support analysis for trio samples (e.g., proband, mother, and father)?

- A.9: Yes, if VCF file contains multiple samples, then trio analysis can be done in `divine.xls`. The next release of Divine will provide more systematic trio analysis.

- Q.10: I don't have a Linux computer. Can Divine run on MacOS or Windows?

- A.10: Not now but we are currently working on the compatibility issue.

- Q.11: What is the memory requirement?

- A.11: In Varant annotation step, memory usage peaks around 1.5GB.

- Q.12: My company wants to use Divine.

- A.12: Divine is free for academic or research purpose. Contact to hongc2@ccf.org for commercial use.

- Q.13: I follow this README file and help. My Divine run fails. Can you help me?

- A.13: Three sample scripts with datasets in `$DIVINE/gcn/bin/prioritize/examples` are available. Use the examples as a template. All three scripts should work if your installation were successful. If the problem is not resolved, send me a log file (e.g., `divine_err.log`) by running,

```
$ your_divine_command 2>&1 | tee divine_err.log
```

- Q.14: Can Divine detect a gene previously never known to be associated with a specific disease?

- A.14: This is a challenging task, but it is essential. Currently, Divine uses gene ontology enrichment and also it uses a system biology approach (heat diffusion model) over STRING protein functional network.

# Change Log
- v.0.1.1 (June 15 2016)
    - Original release
- v.0.1.2 (Aug 8 2018)
    - Minor fix
    - protein domain pathogenicity
    - enhanced the new gene-to-disease association discovery
    
# License
GNU GENERAL PUBLIC LICENSE
https://www.gnu.org/licenses/gpl-3.0.en.html

# Disclaimer
Not intended for direct clinical application. Divine suggests an order of genes to be inspected so that it can make molecular diagnosis efficient. The validation is the responsibility of the user. Neither Divine developer nor any software module integrated is responsible for clinical actions that may result from the use of this software. By using this tool, the user assumes all responsibility for any information that may be generated.

# Reference
- OMIM: Online Mendelian Inheritance in Man, OMIM®. McKusick-Nathans Institute of Genetic Medicine, Johns Hopkins University (Baltimore, MD), {date}. World Wide Web URL: http://omim.org/
- Sebastian Köhler, Sandra C Doelken, Christopher J. Mungall, Sebastian Bauer, Helen V. Firth, et al. "The Human Phenotype Ontology project: linking molecular biology and disease through phenotype data", Nucl. Acids Res. (1 January 2014) 42 (D1): D966-D974 doi:10.1093/nar/gkt1026
- HGMD: https://www.qiagenbioinformatics.com/files/user_manuals/HGMD_User_Manual_2015.2.pdf
- CADD: http://cadd.gs.washington.edu/
- Varant: http://compbio.berkeley.edu/proj/varant/publication.html
- Divine: Prioritizing Genes for Rare Mendelian Disease in Whole Exome Sequencing Data, Changjin Hong and TaeHyun Hwang (Oxford Bioinformatics Application Note)

# Contact
- Changjin Hong, Ph.D (hongc2{at}ccf.org)

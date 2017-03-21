# CNV and SV papers in NGS
Relevant studies with Structual Variants and Copy Number Variants in NGS (Genome, Exome and Target Sequencing) pipelines.

## Background CNV and SV
* [2016 Frequency and Complexity of De Novo Structural Mutation in Autism](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4833290/)

## WES

## WGS
* [FindSV](https://github.com/J35P312/FindSV): FindSV is a structural variation pipeline written in nextflow and python. FindSV performs variant calling using TIDDIT and CNVnator, and Manta.
* [SvABA](https://github.com/walaj/svaba): Structural variation and indel analysis by assembly.
* [2012 DELLY: structural variant discovery by integrated paired-end and split-read analysis](https://academic.oup.com/bioinformatics/article/28/18/i333/245403/DELLY-structural-variant-discovery-by-integrated) - [github](https://github.com/dellytools/delly): Delly2 was the best sv caller in the [DREAM challenge]()
* [2014 LUMPY: a probabilistic framework for structural variant discovery](http://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-6-r84) - [github](https://github.com/arq5x/lumpy-sv)
* [2015 Wham: Identifying Structural Variants of Biological Consequence](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004572) - [github](https://github.com/zeeev/wham) - [mergeSVcallers](https://github.com/zeeev/mergeSVcallers)
* [2016 Manta: rapid detection of structural variants and indels for germline and cancer sequencing applications](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btv710) - [github](https://github.com/Illumina/manta)
* [2017 SV2: Accurate Structural Variation Genotyping and De Novo Mutation Detection](http://biorxiv.org/content/early/2017/03/17/113498) - [github](https://github.com/dantaki/SV2)
* [2017 GRIDSS: sensitive and specific genomic rearrangement detection using positional de Bruijn graph assembly](http://biorxiv.org/content/early/2017/02/21/110387) - [github](https://github.com/PapenfussLab/gridss)
## TS 

## Annotation
* [2015 SpeedSeq: ultra-fast personal genome analysis and interpretation](http://www.nature.com/nmeth/journal/v12/n10/full/nmeth.3505.html) - SVTyper - [github](https://github.com/hall-lab/svtyper)
* [2016 SVScore: an impact prediction tool for structural variation](https://academic.oup.com/bioinformatics/article/doi/10.1093/bioinformatics/btw789/2748212/SVScore-an-impact-prediction-tool-for-structural) - [github](https://github.com/lganel/SVScore)

## Visualization
* [CNVplot](https://github.com/dantaki/CNVplot): Plot CNV data with a genome viewer in R.
* [cnvgram](https://github.com/cc2qe/cnvgram): Draw CNV diagrams.
* [Stupid Simple Structural Variant View](https://github.com/ryanlayer/svv): A two-step process that can help visualize the coverage near a variant across multiple BAMs.
* [CNView](https://github.com/RCollins13/CNView): Visualization, quantitation, and annotation of CNVs from population-scale whole-genome sequencing data.
* [2015 svviz: a read viewer for validating structural variants](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btv478) - [github](https://github.com/svviz/svviz)
* [2016 Ribbon: Visualizing complex genome alignments and structural variation](http://biorxiv.org/content/early/2016/10/20/082123) - [github](https://github.com/MariaNattestad/Ribbon)

## Others
* [SVsim](https://github.com/GregoryFaust/SVsim): a tool that generates synthetic Structural Variant calls as benchmarks to test/evaluate SV calling pipelines.
* [Structural Variant Catalog](https://github.com/tobiasrausch/svcatalog): A repository for human genetic structural variants (SVs) discovered by Delly in the 1000 Genomes cohort of samples.
* [SVDB](https://github.com/J35P312/SVDB): SVDB is a toolkit for constructing and querying structural variant databases. The databases are constructed using the output vcf files from structural variant callers such as TIDDIT, Manta, Fermikit or Delly. [The thousand genomes structural variant calls may also be used as a database](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/).

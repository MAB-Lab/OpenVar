# OpenVar: annotation of genetic variants using deep genome annotation
OpenVar is a deep genome annotation tool. OpenVar currently supports usual annotations (Ensembl and NCBI RefSeq), as well as the OpenProt deep open reading frame (ORF) annotation (www.openprot.org). However, the package can be used with any genome annotation if supplemented with the adequate SnpEff database (see http://pcingola.github.io/SnpEff/se_buildingdb/ for information on how to build a custom SnpEff database).


## Installation of OpenVar

To install the OpenVar package, you will need to clone this repository on your computer. On the main page of this repository, click the **Clone** button. Click on **HTTPS** to clone with https, or **SSH** to clone with ssh, and copy the link. 
On your computer, go to the desired directory and type:
    <code>git clone [the link of the repository]</code>

For more information on how to clone a github repository, see https://docs.github.com/en/github/creating-cloning-and-archiving-repositories/cloning-a-repository-from-github/cloning-a-repository


## Using OpenVar

Depending on your architecture, you may want to edit the paths for the genome files at the beginning of the openvar.py file. All the necessary data for OpenVar to run are within the **data** folder of this repository.

Once you have edited the necessary paths, you can start using OpenVar in python. Below are basic commands to get you started, but see below for a detailed description of each class.

First, you will need to import openvar:
    <code>from OpenVar.openvar import *</code>

Then, you can create a SeqStudy object with the following command:
    <code>vcf = SeqStudy(data_dir = 'path/to/data/directory', 
  file_name = 'filename.vcf', 
  results_dir = 'path/to/results/directory', 
  study_name = 'studyname', 
  specie = 'human', 
  genome_version = 'hg38', 
  annotation = 'OP_Ensembl', 
  picard_path = 'path/to/picard/directory')</code>
  
Then, create an OpenVar object with the following command:
    <code>opv = OpenVar(snpeff_path = 'path/to/snpeff/', 
  vcf = vcf)</code>

Then, to run OpenVar with a classical annotation (Ensembl or NCBI RefSeq), use the following command:
    <code>opv.run_snpeff(vcf.single_vcf_path, vcf.annotation)</code>
If using the OpenProt annotation, use the following command:
    <code>opv.run_snpeff_parallel_pipe()</code>

Then, generate a report with the following command:
    <code>opvr = OPVReport(opv)</code>
    <code>opvr.aggregate_annotated_vcf()</code>
    <code>opvr.write_tabular()</code>
    
If using a classical annotation (Ensembl or NCBI RefSeq), you can generate report statistics and figures with the following command:
    <code>opvr.compute_summary_stats()</code>
If using the OpenProt annotation, use the following command:
    <code>opvr.compute_chrom_gene_level_stats(write_summary_pkl = True)</code>

The last command will generate a pickle object identical to the one generated when using the OpenVar web-based application. This allow you to quickly go back to previous analyses and see general statistics as presented on the Results page of the OpenVar web-based application.
To load the pickle object, simply run the following command:
    <code>import pickle
  summary_path = 'path/to/summary/pickle'
  pickle.load(open(summary_path, 'rb'))</code>


If you have any question regarding OpenVar, don't hesitate to contact us: https://openprot.org/p/ng/contactUs

### Input file format
An example input file can be found [here](https://github.com/MAB-Lab/OpenVar_WebApp/blob/main/vcf_cosmic_HEY2.vcf)
The expected input format is a Variant Call Format ([VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf). It is the de facto standard file format for genomic variants. Other formats should be converted to a VCF format, below are a few examples to run on a shell.
#### BED files
BED input files can be converted using [PLINK](https://www.cog-genomics.org/plink/) with the following command:
<code>plink --bfile [filename] --recode vcf --out [vcf name]</code>
#### dbSNP identifiers
In order to produce a VCF input from a list of dbSNP identifiers, download the VCF file containing all variants within dbSNP [here](https://ftp.ncbi.gov/snp/organisms/human_9606/VCF/)
Then use the following command: <code>grep -wFf dbsnp_id_list.txt my_vcf.vcf > /path_to_output_folder</code>

<div align="center">
    <h1>Custom 16S database for Kraken</h1>
    <br />
    <h2>Script which combines text files with 16S data to multi-FASTA file, compatible with Kraken</h2>
    <br />
</div>


## General information
* **Author:** Kaitlin Weber
* **Commissioned by:** Rijksinstituut voor Volksgezondheid en Milieu (RIVM)

## About this script

This python script combines two in-house databases, into one multi-FASTA file. In the multi-FASTA file, by giving an email (linked to NCBI account) as input and the scientific name from the files, the taxonomy ID is placed in the header. The taxonomy ID is the reason this database can be made compatible with Kraken 2.

This script was designed for creating a database for the purpose of using it with [ExId16S](https://github.com/Kaitlinweber/exid16s). 

## Installation

1. Clone the repository:

```
git clone https://github.com/Kaitlinweber/combine_16S_database
```

2. Enter the directory with the pipeline and install the conda environment:

```
cd combine_16S_database
conda env install -f envs/database_env.yaml
```

### Required parameters

* ```-i, --input```  Kraken summary kreport from ExId16S results, based on the custom made database.
* ```-w, --wgs``` Kraken summary kreport obtained from Kraken2 and Bracken with WGS data.
* ```-s, --sample_list``` Excel sheet with 16S rDNA Sanger sequencing results, with the sample number, genus name and species name in a seperate column. 
* ```-o, --output``` Pathway to output directory, if directory does not exists, directory will be created


### Optional parameters

* ```-si, --silva``` For the comparison of the created database with the SILVA database, the pathway to the SILVA database Kraken kreports can be given here


### The base command to run this script 

```
python  -i [path/to/input/dir/ExId16S] -w [path/to/input/dir/WGS] -s [path/to/input/dir/16S_sequencing] -o [path/to/output/dir] 
```

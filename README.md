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

* ```-db, --database```  Path to directory, which can have multiple subdirectories with 16S rDNA information. The files need to contain ORGANSIM/DEFENITION (for extraction scientific name), ORIGIN (for extraction sequence) (format derived from Genbank format, but in text file format)
* ```-e, --email``` Email which is linked to a NCBI account, this is for obtaining the taxonomy ID, which is required for building a Kraken 2 database.
* ```-o, --output``` Pathway to output directory, if directory does not exists, directory will be created



### The base command to run this script 

```
python combine_16S_database.py -db [path/to/input/dir/16Sdata] -o [email from NCBI] -o [path/to/output/dir] 
```

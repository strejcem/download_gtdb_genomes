# download_gtdb_genomes
Downloads NCBI assemblies from GenBank based on GTDB taxonomy labels.
## Requirements:
### Packages:
```bash
conda create -n download_gtsb_genomes r-tidyverse r-argparser ncbi-datasets-cli>=15.11.0
conda activate download_gtsb_genomes
```
### Database:
```bash
wget https://data.gtdb.ecogenomic.org/releases/latest/ar53_metadata.tar.gz
wget https://data.gtdb.ecogenomic.org/releases/latest/bac120_metadata.tar.gz
wget https://data.gtdb.ecogenomic.org/releases/latest/VERSION.txt

tar xzf ar53_metadata.tar.gz
tar xzf bac120_metadata.tar.gz
cat *_metadata*.tsv | gz > gtdb_metadata.tsv.gz

```
### Example of usage:
```bash
Rscript download_gtdb_genomes.R -h
```

Download five random genomes of each `class` in both Archaeal and Bacterial domains. Filter the genomes to only representatives, high quality, rename the downloaded genomes to simplified names, i.e. 'accession.fna', write a file that maps all the scaffolds/contigs to 'c2g.tsv', and make the whole procedure reproducible by setting `--seed` to '123'.
```bash
Rscript download_gtdb_genomes.R --metadata gtdb_metadata.tsv.gz \
                                --out_dir gtdb_genomes \
                                --representatives TRUE \
                                --gtdb_tag d__Archaea,d__Bacteria \
                                --num_genomes 5 \
                                --tax_level class \
                                --contigs2genomes c2g.tsv
                                --seed 123 \
                                --mimag high \
                                --rename TRUE
```

To get only file of accessions, use `-j/--just_acc_list`.

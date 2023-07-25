#!/usr/bin/env Rscript

library(argparser, quietly = TRUE)

p <- arg_parser("Version 1.21, by Michal Strejcek @ UCT Prague
                Downloads NCBI assemblies from RefSeq or GenBank based on GTDB taxonomy labels.
		
                At least one of Arc/Bac taxonomy/metadata file needs to be specified.
                Taxonomy/metadata files can be accessed at https://data.gtdb.ecogenomic.org/releases/latest/
                If metadata file is used, you can filter for gtdb_representative genomes only, i.e. genomes used in GTDB reference tree.

		Depends on:
	       		'conda isntall r-tidyverse r-argparser ncbi-genome-download'
	                'wget https://data.gtdb.ecogenomic.org/releases/latest/ar53_metadata.tar.gz'
	                'wget https://data.gtdb.ecogenomic.org/releases/latest/bac120_metadata.tar.gz'
			'wget https://data.gtdb.ecogenomic.org/releases/latest/VERSION.txt'
		")
p <- add_argument(p, "--arc_tax", help = "File location for 'ar122_taxonomy_rXX.tsv(.gz)'", default = NA)
p <- add_argument(p, "--bac_tax", help = "File location for 'bac120_taxonomy_rXX.tsv(.gz)'", default = NA)
p <- add_argument(p, "--arc_meta", help = "File location for 'ar122_metadata_rXX.tsv(.gz)'", default = NA)
p <- add_argument(p, "--bac_meta", help = "File location for 'bac120_metadata_rXX.tsv(.gz)'", default = NA)
p <- add_argument(p, "--out_dir", help = "Output directory", default = NA)
p <- add_argument(p, "--tag_GTDB", help = "Comma-separated list of gtdb taxonomy
                  e.g. 'p__Desantisbacteria,g__Escherichia,g__Bacilli'", default = NA)
p <- add_argument(p, "--representative", help = "Download only GTDB representative genomes.", flag = TRUE)
p <- add_argument(p, "--keep_intermediate", help = "Keep intermediate files.", flag = TRUE)
argv <- parse_args(p)

suppressPackageStartupMessages(library(tidyverse))

out_dir <- file.path(argv$out_dir)

if(is.na(argv$tag_GTDB)) {
  stop("Taxonomy tags must be specified!")
}

if(is.na(out_dir)) {
  stop("Output directory must be specified!")
}

at <- if_else(is.na(argv$arc_tax), FALSE, TRUE)
am <- if_else(is.na(argv$arc_meta), FALSE, TRUE)
bt <- if_else(is.na(argv$bac_tax), FALSE, TRUE)
bm <- if_else(is.na(argv$bac_meta), FALSE, TRUE)
arc <- any(at, am)
bac <- any(bt, bm)

if(!any(at, am, bt ,bm)) {
  stop("At least one of taxonomy/metadata files must be specified!")
}

if(all(at, am) | all(at, bm) | all(bt, bm) | all(bt, am)) {
  stop("Please use either taxonomy or metadata file(s) only!")
}

to_pick <- str_split(argv$tag_GTDB, pattern = ",")[[1]]

#prepare GTDB taxonomy database
if(am) {
  gtdb_arc <-read_tsv(file.path(argv$arc_meta),
                      col_types = cols_only(
                        accession = 'c',
                        gtdb_taxonomy = 'c',
                        gtdb_representative = 'l')
                      ) %>% 
    separate(
      gtdb_taxonomy,
      into = c(
        "Domain",
        "Phylum",
        "Class",
        "Order",
        "Family",
        "Genus",
        "Species"
      ),
      sep = ";"
    )
} else if(at) {
  gtdb_arc <-
    read_tsv(file.path(argv$arc_tax),
             col_names = c("accession", "tax"),
             col_types = cols()) %>%
    separate(
      tax,
      into = c(
        "Domain",
        "Phylum",
        "Class",
        "Order",
        "Family",
        "Genus",
        "Species"
      ),
      sep = ";"
    )
} else {
  gtdb_arc <- NULL
}

if(bm) {
  gtdb_bac <-read_tsv(file.path(argv$bac_meta),
                      col_types = cols_only(
                        accession = 'c',
                        gtdb_taxonomy = 'c',
                        gtdb_representative = 'l')
                      ) %>% 
    separate(
      gtdb_taxonomy,
      into = c(
        "Domain",
        "Phylum",
        "Class",
        "Order",
        "Family",
        "Genus",
        "Species"
      ),
      sep = ";"
    )
} else if(bt) {
  gtdb_bac <-
    read_tsv(file.path(argv$bac_tax),
             col_names = c("accession", "tax"),
             col_types = cols()) %>%
    separate(
      tax,
      into = c(
        "Domain",
        "Phylum",
        "Class",
        "Order",
        "Family",
        "Genus",
        "Species"
      ),
      sep = ";"
    )
} else {
  gtdb_bac <- NULL
}

gtdb_tax <- bind_rows(gtdb_arc, gtdb_bac) %>% 
  mutate(
    section = case_when(str_detect(accession, "^GB_") ~ "GB",
                        str_detect(accession, "^RS_") ~ "RS"),
    accession = str_remove(accession, "^(GB|RS)_")
  )

tax_selection <- gtdb_tax %>%
  {if(argv$representative) filter(., gtdb_representative) else .} %>%
  select(accession, section, Domain, Phylum, Class, Order, Family, Genus, Species) %>% 
  pivot_longer(names_to = "taxonomy", cols = -c(accession, section)) %>%
  filter(value %in% to_pick) %>%
  distinct(accession, section)

num_acc <- nrow(tax_selection)
message("Found ", num_acc, " assembly accessions.")

if (num_acc == 0) {
  stop("No accession to download!")
}

# filter what is GB and what RS
GB <- nrow(filter(tax_selection, section == "GB"))
RS <- nrow(filter(tax_selection, section == "RS"))

NCBI_groups <- case_when(
  arc & bac ~ "archaea,bacteria",
  !arc & bac ~ "bacteria",
  arc & !bac ~ "archaea"
)
rng <- sample(10000:99999, 1)
tmp_dir <-file.path(str_c("tmp_", rng))
dir.create(tmp_dir)

#download GenBank
if (GB > 0) {
  message("Downloading ", GB, " genomes from GenBank")
  GB_file <- file.path(tmp_dir, "accession_to_download_GB.txt")
  tax_selection %>%
    filter(section == "GB") %>%
    select(accession) %>%
    write_tsv(file = GB_file, col_names = FALSE)
  
  ncbi_genome_download_arguments <- str_c(sep = " ",
    "-A", GB_file, "-s genbank -F fasta -p 4 -o", tmp_dir, "--no-cache", "--debug", NCBI_groups
  )
  system2("ncbi-genome-download", args = ncbi_genome_download_arguments)
}

#download RefSeq
if (RS > 0) {
  message("Downloading ", RS, " genomes from RefSeq")
  RS_file <- file.path(tmp_dir, "accession_to_download_RS.txt")
  tax_selection %>%
    filter(section == "RS") %>%
    select(accession) %>%
    write_tsv(file = RS_file, col_names = FALSE)
  
  ncbi_genome_download_arguments <- str_c(sep = " ",
    "-A", RS_file, "-s refseq -F fasta -p 4 -o", tmp_dir, "--no-cache", "--debug", NCBI_groups
  )
  
  system2("ncbi-genome-download", args = ncbi_genome_download_arguments)
}

#move genome files to a new directory
dir.create(out_dir, recursive = TRUE)
ncbi_genomes <-
  list.files(
    tmp_dir,
    pattern = "*fna.gz",
    full.names = TRUE,
    recursive = TRUE
  )
only_accession <-
  str_extract(basename(ncbi_genomes), "(GCA|GCF)_[:digit:]+\\.[:digit:]+")
new_names <- str_c(only_accession, ".fna.gz")
file.copy(
  from = file.path(ncbi_genomes),
  to = file.path(out_dir, new_names),
  overwrite = TRUE
)

if(!argv$keep_intermediate) {
  unlink(tmp_dir, recursive = TRUE)
}

#check that everything went ok
downloaded_files <-
  list.files(out_dir, pattern = "*fna.gz", full.names = FALSE)
downloaded_accession <- str_remove(downloaded_files, ".fna.gz")
missing_genomes <- tax_selection %>%
  filter(!accession %in% downloaded_accession) %>% 
  pull(accession)


if (length(missing_genomes) == 0) {
  message("All genomes succesfully downloaded!")
} else {
  message("Missing these genomes:")
  gtdb_tax %>% 
    filter(accession %in% missing_genomes) %>% 
    select(accession, Domain, Phylum, Genus) %>% 
    print(n=20)
}

#!/usr/bin/env Rscript


library(argparser, quietly = TRUE)


p <- arg_parser(hide.opts = TRUE,
  "Version 1.3, by Michal Strejcek @ UCT Prague
  Downloads NCBI assemblies from RefSeq or GenBank based on GTDB taxonomy labels.

  At least one of Arc/Bac taxonomy/metadata file needs to be specified.
  Taxonomy/metadata files can be accessed at https://data.gtdb.ecogenomic.org/releases/latest/
  If metadata file is used, you can filter for gtdb_representative genomes only, i.e. genomes used in GTDB reference tree.

  Depends on:
  'conda install r-tidyverse r-argparser ncbi-genome-download'
  'wget https://data.gtdb.ecogenomic.org/releases/latest/ar53_metadata.tar.gz'
  'wget https://data.gtdb.ecogenomic.org/releases/latest/bac120_metadata.tar.gz'
  'wget https://data.gtdb.ecogenomic.org/releases/latest/VERSION.txt'
"
)
p <- add_argument(p, "--arc_tax", help = "File location for 'ar122_taxonomy_rXX.tsv(.gz)'", default = NA)
p <- add_argument(p, "--bac_tax", help = "File location for 'bac120_taxonomy_rXX.tsv(.gz)'", default = NA)
p <- add_argument(p, "--arc_meta", help = "File location for 'ar122_metadata_rXX.tsv(.gz)'", default = NA)
p <- add_argument(p, "--bac_meta", help = "File location for 'bac120_metadata_rXX.tsv(.gz)'", default = NA)
p <- add_argument(p, "--out_dir", help = "Output directory", default = NA)
p <- add_argument(p, "--gtdb_tag", help = "Comma-separated list of gtdb taxonomy
                  e.g. 'p__Desantisbacteria,g__Escherichia,g__Bacilli'", default = NA)
p <- add_argument(p, "-F", help = "[fasta|protein-fasta|...] Argument for 'ncbi-genome-download',
                  see its help for more options.", default = "fasta")
p <- add_argument(p, "--representatives", help = "[TRUE|FALSE] Download only GTDB representative
                  or all genomes.", default = TRUE)
p <- add_argument(p, "--keep_intermediate", help = "Keep intermediate files. Includes accession lists
                  and donwloaded files", defau = TRUE)
p <- add_argument(p, "--num_genomes", help = "[integer] Pick specific number of random genomes
                  per the '--gtdb_tag' entry.", default = NA, type = "integer")
p <- add_argument(p, "--tax_level", help = "One of [phylum, class, order, family,
                  genus, species]. To be used with '--num_genomes'. Move the 
                  '--num_genomes' selection to a specified taxnomic rank within the '--gtdb_tag'
                  entries. Eg. '--gtdb_tag f__Pseudomonadaceae,c__Heimdallarchaeia --num_genomes 3
                  --tax_level genus' will download three random genomes of all genera within
                  Pseudomonadaceae family and  Heimdallarchaeia class.", default = NA)
p <- add_argument(p, "--just_list", help = "Creates the accession list but no download.", flag = TRUE)
p <- add_argument(p, "--seed", help = "Seed number for reproducible random sampling", type = "integer", default = NA)

argv <- parse_args(p)

suppressPackageStartupMessages(library(tidyverse))

out_dir <- file.path(argv$out_dir)

if(is.na(argv$gtdb_tag)) {
  stop("Taxonomy tags must be specified!
       Use '-h' for help!")
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

if (!is.na(argv$num_genomes) & 
    !(argv$num_genomes) > 0 
    ) {
  stop("--num_genomes must be > 0!")
}


if (!is.na(argv$num_genomes) & 
    (!as.integer(argv$num_genomes) > 0 | !is.integer(argv$num_genomes))
) {
  stop("--num_genomes must be > 0!")
}

if (!is.na(argv$tax_level) &
    !argv$tax_level %in% c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {
  stop("'--tax_level' must be one of [Phylum|Class|Order|Family|Genus|Species]")
}

if (!is.na(argv$tax_level) &
    is.na(argv$num_genomes)) {
  stop("When using '--tax_level' you must also specify '--num_genomes'!")
}

download_format <- str_split_1(argv$F, pattern = ",")
if (!all(download_format %in% c("genbank", "fasta", "rm", "features", "gff", "protein-fasta", "genpept","wgs",
  "cds-fasta", "rna-fna", "rna-fasta", "assembly-report", "assembly-stats", "all"))) {
  stop("-F must be one of ['genbank', 'fasta', 'rm', 'features', 'gff', 'protein-fasta', 
       'genpept','wgs', 'cds-fasta', 'rna-fna', 'rna-fasta', 'assembly-report',
       'assembly-stats', 'all']!")
}

if (!is.na(argv$seed)) {
  set.seed(argv$seed)
}

to_pick <- str_split_1(argv$gtdb_tag, pattern = ",")

#prepare GTDB taxonomy database
if(am) {
  message("Reading metadata for Archaea.")
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
  message("Done.")
} else if(at) {
  message("Reading taxonomy for Archaea.")
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
  message("Done.")
} else {
  message("Skipping Archaea!")
  gtdb_arc <- NULL
}

if(bm) {
  message("Reading metadata for Bacteria")
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
  message("Done.")
} else if(bt) {
  message("Reading taxonomy for Bacteria.")
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
  message("Done.")
} else {
  message("Skipping Bacteria!")
  gtdb_bac <- NULL
}

gtdb_tax <- bind_rows(gtdb_arc, gtdb_bac) %>% 
  mutate(
    section = case_when(str_detect(accession, "^GB_") ~ "GB",
                        str_detect(accession, "^RS_") ~ "RS"),
    accession = str_remove(accession, "^(GB|RS)_")
  )

rm(gtdb_arc, gtdb_bac)

message("Selecting accessions...")
tax_selection <- gtdb_tax %>%
  {
    if (argv$representatives) {
      filter(., gtdb_representative)
    } else {
      .
    }   
  } %>%
  select(accession,
         section,
         Domain,
         Phylum,
         Class,
         Order,
         Family,
         Genus,
         Species) %>%
  pivot_longer(names_to = "taxonomy", cols = -c(accession, section)) %>%
  filter(value %in% to_pick) %>%
  distinct(accession) %>% 
  left_join(gtdb_tax, by = "accession") %>% 
  {
    if (!is.na(argv$tax_level)) {
      select(., accession, section, matches(argv$tax_level)) %>% 
      group_by(across(argv$tax_level)) %>% 
      slice_sample(n = argv$num_genomes) %>% 
      ungroup()
    } else if (!is.na(argv$num_genomes)) { 
        slice_sample(., n = argv$num_genomes)
    } else {
      .
    }
  } %>% 
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

  GB_file <- file.path(tmp_dir, "accession_to_download_GB.txt")
  
  tax_selection %>%
    filter(section == "GB") %>%
    select(accession) %>%
    write_tsv(file = GB_file, col_names = FALSE)
  
  if (!argv$just_list) {
    message("Downloading ", GB, " genomes from GenBank") 
    ncbi_genome_download_arguments <- str_c(sep = " ",
      "-A", GB_file, "-s genbank -F", argv$F, "-p 4 -o", tmp_dir, "--no-cache", NCBI_groups
    )
    system2("ncbi-genome-download", args = ncbi_genome_download_arguments)
  }
}

#download RefSeq
if (RS > 0) {
  
  RS_file <- file.path(tmp_dir, "accession_to_download_RS.txt")
  tax_selection %>%
    filter(section == "RS") %>%
    select(accession) %>%
    write_tsv(file = RS_file, col_names = FALSE)
  
  if (!argv$just_list) {
    message("Downloading ", RS, " genomes from RefSeq")
    ncbi_genome_download_arguments <- str_c(sep = " ",
      "-A", RS_file, "-s refseq -F", argv$F, "-p 4 -o", tmp_dir, "--no-cache", NCBI_groups
    )
    system2("ncbi-genome-download", args = ncbi_genome_download_arguments)
  }
}

dir.create(out_dir, recursive = TRUE)
if (argv$just_list) {
  invisible(
    file.copy(
      from = GB_file,
      to = file.path(out_dir, "accession_to_download_GB.txt")
    ))
  invisible(
    file.copy(
      from = RS_file,
      to = file.path(out_dir, "accession_to_download_RS.txt")
    )
  )
  unlink(tmp_dir, recursive = TRUE)
  message("All done!")
} else {
  #move genome files to a new directory
  
  ncbi_genomes <-
    list.files(
      tmp_dir,
      pattern = "\\.gz$",
      full.names = TRUE,
      recursive = TRUE
    )
  
  new_names <- ncbi_genomes %>% 
    basename() %>% 
    str_replace("_[:alpha:]+[:digit:]+v[:digit:]+_", "_")
  
  invisible(
    file.copy(
      from = file.path(ncbi_genomes),
      to = file.path(out_dir, new_names),
      overwrite = FALSE
    )
  )
  
  if(!argv$keep_intermediate) {
    unlink(tmp_dir, recursive = TRUE)
  }
  
  #check that everything went ok
  downloaded_files <-
    list.files(out_dir, pattern = "*.gz", full.names = FALSE)
  downloaded_accession <- str_extract(basename(downloaded_files), "(GCA|GCF)_[:digit:]+\\.[:digit:]+")
  missing_genomes <- tax_selection %>%
    filter(!accession %in% downloaded_accession) %>% 
    pull(accession)
  
  
  if (length(missing_genomes) == 0) {
    message("***!!! All genomes succesfully downloaded !!!***")
  } else {
    message("Missing these genomes:")
    gtdb_tax %>% 
      filter(accession %in% missing_genomes) %>% 
      select(accession, Domain, Phylum, Genus) %>% 
      print(n=20)
  }
} 
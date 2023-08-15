#!/usr/bin/env Rscript

library(argparser, quietly = TRUE)

p <- arg_parser(
  hide.opts = TRUE,
  "Version 2.2, by Michal Strejcek @ UCT Prague
  Downloads NCBI assemblies from RefSeq or GenBank based on GTDB taxonomy labels.

  At least one of Arc/Bac taxonomy/metadata file needs to be specified.
  Taxonomy/metadata files can be accessed at https://data.gtdb.ecogenomic.org/releases/latest/

  Depends on:
  'conda install r-tidyverse r-argparser r-furrr ncbi-datasets-cli>=15.11.0'
  'wget https://data.gtdb.ecogenomic.org/releases/latest/ar53_metadata.tar.gz'
  'wget https://data.gtdb.ecogenomic.org/releases/latest/bac120_metadata.tar.gz'
  'wget https://data.gtdb.ecogenomic.org/releases/latest/VERSION.txt'
"
)

p <-
  add_argument(p, "--metadata", help = "File location for 'metadata_rXX.tsv(.gz)'", default = "gtdb_metadata_r*.tsv.gz")
p <-
  add_argument(p, "--out_dir", help = "Output directory", default = NA)
p <-
  add_argument(p,
               "--gtdb_tag",
               help = "Comma-separated list of gtdb taxonomy, e.g. 'p__Desantisbacteria,g__Escherichia,g__Bacilli'",
               default = NA)
p <-
  add_argument(p,
               "--representatives",
               help = "[TRUE|FALSE] Download only GTDB representative or all genomes.",
               default = TRUE)
p <-
  add_argument(p,
               "--num_genomes",
               help = "Pick [integer] number of random genomes per the '--gtdb_tag' entry.",
               default = NA,
               type = "integer")
p <-
  add_argument(p,
               "--tax_level",
               help = "One of [phylum, class, order, family, genus, species]. To be used with '--num_genomes'. Move the '--num_genomes' selection to a specified taxnomic rank within the '--gtdb_tag' entries. Eg. '--gtdb_tag f__Pseudomonadaceae,c__Heimdallarchaeia --num_genomes 3 --tax_level genus' will download three random genomes of all genera within Pseudomonadaceae family and Heimdallarchaeia class.",
               default = NA)
p <-
  add_argument(p,
               "--contigs2genomes",
               help = "[Path] to a file containing contigs/scaffolds to their respective genome file.",
               default = NA)
p <-
  add_argument(p,
               "--out_tax",
               help = "[Path] to a file containing taxonomy for the accessions.",
               default = NA)
p <-
  add_argument(p,
               "--mimag",
               help = "One of [low, medium, high, all] filters MIMAG quality of MAGs",
               default = "all")
p <-
  add_argument(p,
               "--rename",
               help = "[TRUE|FALSE] Rename downloaded genome files to accession.fna and a sequence headers starts with the assembly accession.",
               default = TRUE)
p <-
  add_argument(p, "--just_acc_list", help = "Creates the accession list but no download.", flag = TRUE)
p <-
  add_argument(p,
               "--seed",
               help = "Seed number for reproducible random sampling",
               type = "integer",
               default = NA)

argv <- parse_args(p)

suppressPackageStartupMessages(library(tidyverse))

out_dir <- file.path(argv$out_dir)

if (is.na(argv$gtdb_tag)) {
  stop("Taxonomy tags must be specified!
       Use '-h' for help!")
}


if (is.na(argv$metadata)) {
  stop("GTDB metadata file must be specified!")
}

if (!is.na(argv$num_genomes) &
    !(argv$num_genomes) > 0) {
  stop("--num_genomes must be > 0!")
}

if (!is.na(argv$num_genomes) &
    (!as.integer(argv$num_genomes) > 0 |
     !is.integer(argv$num_genomes))) {
  stop("--num_genomes must be > 0!")
}

if (!is.na(argv$tax_level) &
    !argv$tax_level %in% c("phylum", "class", "order", "family", "genus", "species")) {
  stop("'--tax_level' must be one of [phylum|class|order|family|genus|species]")
}

if (!is.na(argv$tax_level) &
    is.na(argv$num_genomes)) {
  stop("When using '--tax_level' you must also specify '--num_genomes'!")
}

if (!is.na(argv$mimag) &
    !argv$mimag %in% c("low", "medium", "high", "all")) {
  stop("'--mimag' must be one of [low|medium|high|all")
}

if (!is.na(argv$seed)) {
  set.seed(argv$seed)
}

if (is.na(out_dir)) {
  stop("Output directory must be specified!")
}

to_pick <- str_split_1(argv$gtdb_tag, pattern = ",")

message("Reading GTDB metadata...")
gtdb_meta <-
  file.path(argv$metadata) %>%
  read_tsv(
    col_types = cols_only(
      accession = 'c',
      gtdb_taxonomy = 'c',
      gtdb_representative = 'l',
      mimag_low_quality = 'l',
      mimag_medium_quality = 'l',
      mimag_high_quality = 'l'
    )
  ) %>%
  separate(
    gtdb_taxonomy,
    into = c(
      "domain",
      "phylum",
      "class",
      "order",
      "family",
      "genus",
      "species"
    ),
    sep = ";"
  ) %>%
  mutate(accession = str_remove(accession, "^(GB|RS)_"))
message("Done.")

message("Selecting accessions...")
tax_selection <- gtdb_meta %>%
  {
    if (argv$representatives) {
      filter(., gtdb_representative)
    } else {
      .
    }
  } %>% {
    if (argv$mimag == "low") {
      filter(., mimag_low_quality)
    } else if (argv$mimag == "medium") {
      filter(., mimag_medium_quality)
    } else if (argv$mimag == "high") {
      filter(., mimag_high_quality)
    } else {
      .
    }
  } %>%
  select(accession,
         domain,
         phylum,
         class,
         order,
         family,
         genus,
         species) %>%
  pivot_longer(names_to = "taxonomy", cols = -c(accession)) %>%
  filter(value %in% to_pick) %>%
  distinct(accession) %>%
  left_join(gtdb_meta, by = "accession") %>%
  {
    if (!is.na(argv$tax_level)) {
      select(., accession, matches(argv$tax_level)) %>%
        group_by(across(argv$tax_level)) %>%
        slice_sample(n = argv$num_genomes) %>%
        ungroup()
    } else if (!is.na(argv$num_genomes)) {
      slice_sample(., n = argv$num_genomes)
    } else {
      .
    }
  } %>%
  distinct(accession)

num_acc <- nrow(tax_selection)
message("Found ", num_acc, " assembly accessions.")

if (num_acc == 0) {
  stop("No accession to download!")
}

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

acc_file <- file.path(out_dir, "accessions.txt")

tax_selection %>%
  select(accession) %>%
  write_tsv(file = acc_file, col_names = FALSE)

if (argv$just_acc_list) {
  message("All done!")
  quit("no", runLast = FALSE)
}

message("Downloading genomes...")
dehydrated_file <- file.path(out_dir, "dehydrated.zip")
dehydrate_arguments <- str_c(
  sep = " ",
  "download",
  "genome",
  "accession",
  "--inputfile",
  acc_file,
  "--dehydrated",
  "--filename",
  dehydrated_file
)
unzip_arguments <- str_c(sep = " ",
                         "-q",
                         dehydrated_file,
                         "-d",
                         out_dir)
rehydrate_arguments <- str_c(sep = " ",
                             "rehydrate",
                             # "--gzip",
                             "--directory", out_dir)
system2("datasets", args = dehydrate_arguments)

if (!file.exists(dehydrated_file)) {
  stop("Something went wrong with the dataset download!")
}

system2("unzip", args = unzip_arguments)
system2("datasets", args = rehydrate_arguments)

#check that everything went ok
downloaded_genomes <-
  list.files(
    out_dir,
    pattern = "*.fna",
    full.names = FALSE,
    recursive = TRUE
  )
downloaded_accession <-
  str_extract(basename(downloaded_genomes),
              "(GCA|GCF)_[:digit:]+\\.[:digit:]+")
missing_genomes <- tax_selection %>%
  filter(!accession %in% downloaded_accession) %>%
  pull(accession)


if (length(missing_genomes) == 0) {
  message("***!!! All genomes succesfully downloaded !!!***")
} else {
  message("Missing these genomes:")
  gtdb_meta %>%
    filter(accession %in% missing_genomes) %>%
    select(accession, domain, phylum, genus) %>%
    print(n = 20)
}

if (argv$rename) {
  ncbi_genomes <-
    list.files(
      out_dir,
      pattern = "\\.fna$",
      full.names = TRUE,
      recursive = TRUE
    )
  
  new_names <- ncbi_genomes %>%
    basename() %>%
    str_extract("^GC._[:digit:]+\\.[:digit:]+") %>%
    str_c(".fna")
  
  invisible(file.copy(
    from = file.path(ncbi_genomes),
    to = file.path(out_dir, new_names),
    overwrite = FALSE
  ))
  
  unlink(recursive = TRUE,
         c(
           file.path(out_dir, "README.md"),
           file.path(out_dir, "accessions.txt"),
           file.path(out_dir, "dehydrated.zip"),
           file.path(out_dir, "ncbi_dataset")
         ))
}

if (!is.na(argv$contigs2genomes)) {
  message("Generating contigs2genomes file...")
  
  if (!dir.exists(dirname(argv$contigs2genomes))) {
    dir.create(dirname(argv$contigs2genomes), recursive = TRUE)
  }
  
  fna_files <-
    list.files(
      out_dir,
      pattern = "*\\.fna$",
      recursive = TRUE,
      full.names = TRUE
    )
  
  collect_contigs <- function(x) {
    awk_cmd <- sprintf("awk '$1 ~ /^>/' %s", x)
    tibble(genomes = basename(x),
           contigs = readLines(pipe(sprintf(awk_cmd))))
  }
  library(furrr)
  plan(multisession)
  future_map(fna_files, .f = collect_contigs, .progress = TRUE) %>%
    list_rbind() %>%
    mutate(
      genomes = str_remove(genomes, "\\.fna$"),
      contigs = str_remove(contigs, " .*"),
      contigs = str_remove(contigs, "^>")
    ) %>%
    write_tsv(argv$contigs2genomes, col_names = FALSE)
  message("Done.")
}

if (!is.na(argv$out_tax)) {
  if (!dir.exists(dirname(argv$out_tax))) {
    dir.create(dirname(argv$out_tax), recursive = TRUE)
  }
  
  tax_selection %>%
    left_join(gtdb_meta, by = "accession") %>%
    select(accession, domain, phylum, class, order, family, genus, species) %>%
    write_tsv(argv$out_tax)
}

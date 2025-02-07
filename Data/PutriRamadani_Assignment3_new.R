# Weekly assignment 3
# Putri Ramadani


# Question 1: file downloaded

# Question 2
# Load the necessary library
library(biomaRt)

# Load gene list
genes <- read.delim("Congenital muscular dystrophy.tsv", sep = "\t", header = TRUE)

# check column names
colnames(genes)


# Question 3
# Connect to the Ensembl GRCh37 dataset
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", 
                   dataset = "hsapiens_gene_ensembl", 
                   host = "https://grch37.ensembl.org")

# Question 3.1
# Get gene coordinates
# Use genes data frame and the correct column name
gene_coordinates <- getBM(
  attributes = c("chromosome_name", "start_position", "end_position", "ensembl_gene_id"),
  filters = "ensembl_gene_id",
  values = genes$EnsemblId.GRch37., 
  mart = ensembl
)

# Save to BED file
write.table(gene_coordinates, "panel.bed", sep = "\t", quote = FALSE, row.names = FALSE)


# Question 3.2
# Get exon coordinates 
exon_data <- getBM(
  attributes = c("chromosome_name", "exon_chrom_start", "exon_chrom_end", "ensembl_exon_id"),
  filters = "ensembl_gene_id",
  values = genes$EnsemblId.GRch37.,
  mart = ensembl
)

# Save to BED file
write.table(exon_data, "panel_exons.bed", sep = "\t", quote = FALSE, row.names = FALSE)


# Question 3.3
# Get transcript and CDS length data
transcripts <- getBM(
  attributes = c("ensembl_gene_id", "ensembl_transcript_id", "transcript_length", "cds_length"),
  filters = "ensembl_gene_id",
  values = genes$EnsemblId.GRch37.,
  mart = ensembl
)

# Find the longest transcript for each gene
longest_transcripts <- transcripts[which.max(transcripts$cds_length), ]
longest_transcripts
# --- output --- #
# ensembl_gene_id      ensembl_transcript_id   transcript_length cds_length
# 282 ENSG00000131018       ENST00000367255             27748      26394
# -------------- #

# Get exons for the longest transcript
longest_exon_data <- getBM(
  attributes = c("chromosome_name", "exon_chrom_start", "exon_chrom_end", "ensembl_exon_id"),
  filters = "ensembl_transcript_id",
  values = longest_transcripts$ensembl_transcript_id,
  mart = ensembl
)

# Save to BED file
write.table(longest_exon_data, "panel_exons_longest_transcript.bed", sep = "\t", quote = FALSE, row.names = FALSE)



# Question 4
# Find the shortest gene with a coding sequence
shortest_gene <- transcripts[which.min(transcripts$transcript_length), ]
shortest_gene

# Get the coding sequence (CDS)
cds_sequence <- getSequence(
  id = shortest_gene$ensembl_transcript_id,
  type = "ensembl_transcript_id",
  seqType = "cdna",
  mart = ensembl
)
# -- output -- #
# ensembl_gene_id       ensembl_transcript_id    transcript_length cds_length
# 698 ENSG00000198947       ENST00000448370               123        123
# ------------ #

# Save as a FASTA file
fasta_header <- paste0(">", shortest_gene$ensembl_transcript_id, " | ", shortest_gene$ensembl_gene_id, " | Shortest CDS")
fasta_content <- paste(fasta_header, cds_sequence$cdna, sep = "\n")

writeLines(fasta_content, "shortest_gene.fasta")



# Question 5: file downloaded


# Question 6 
# Determine how many variants in the vcf PASS quality filters
# Load necessary libraries
library(dplyr)
library(readr)
library(VariantAnnotation)

# Load the VCF file
library(VariantAnnotation)

# Read VCF file
vcf_file <- "PGPC_0001_S1.flt.subset.vcf"
vcf <- readVcf(vcf_file, "hg19")

# View the structure of the VCF
vcf

# Filter for variants that pass the quality filter
variants_pass <- rowRanges(vcf)[rowRanges(vcf)$FILTER == "PASS"]

# Get the number of PASS variants
num_pass_variants <- length(variants_pass)
cat("Number of variants that pass the quality filter:", num_pass_variants, "\n")
# Number of variants that pass the quality filter: 42535 


# Question 7
# Load necessary libraries
library(VariantAnnotation)
library(GenomicRanges)

# Read the VCF file
vcf_file <- "PGPC_0001_S1.flt.subset.vcf"
vcf <- readVcf(vcf_file, "hg19")

# Load and create the panel region (GRanges object from BED file)
panel <- read.table("panel.bed", sep = "\t", header = TRUE)
panel_region <- GRanges(
  seqnames = panel$chromosome_name,
  ranges = IRanges(start = panel$start_position, end = panel$end_position),
  names = panel$ensembl_gene_id
)

# Check chromosome names in the panel
seqlevels(panel_region)

# Check chromosome names in the VCF
seqlevels(rowRanges(vcf))

# Update chromosome names to UCSC style if necessary
seqlevelsStyle(panel_region) <- "UCSC"

# Verify that the chromosome names are now compatible
seqlevels(panel_region)
seqlevels(rowRanges(vcf))


# Subset VCF based on overlaps with panel regions
overlaps <- findOverlaps(rowRanges(vcf), panel_region)
filtered_variants <- vcf[queryHits(overlaps)]

# Expand variants with multiple ALT alleles 
expanded_variants <- expand(filtered_variants)

# Filter for PASS variants
pass_variants <- expanded_variants[rowRanges(expanded_variants)$FILTER == "PASS"]
pass_variants

# Extract the necessary data from the pass_variants object
CHROM <- as.character(seqnames(pass_variants))
POS <- start(ranges(pass_variants))

# Create the Data Frame
variant_data <- data.frame(
  Coordinate = paste0(CHROM, ":", POS),
  REF = as.character(ref(pass_variants)),
  ALT = as.character(unlist(alt(pass_variants))),
  QUAL = qual(pass_variants),
  FILTER = rowRanges(pass_variants)$FILTER
)

# Save the TSV File 
write.table(variant_data, "variants.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


# Question 8 
# Load Required Libraries 
library(VariantAnnotation)
library(GenomicRanges)
library(MafDb.TOPMed.freeze5.hg19)

# Load the Filtered Variants 
# Read the TSV file generated in Question 7
variant_data <- read.table("variants.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Extract CHROM and POS from the Coordinate column
variant_data$CHROM <- sub(":.*", "", variant_data$Coordinate)
variant_data$POS <- as.integer(sub(".*:", "", variant_data$Coordinate))

# Create a GRanges Object for TopMed Annotation 
variant_ranges <- GRanges(
  seqnames = variant_data$CHROM,
  ranges = IRanges(start = variant_data$POS, end = variant_data$POS), # check again here
  REF = variant_data$REF,
  ALT = variant_data$ALT
)

# Ensure chromosome naming style matches between the variant ranges and the MAF database
seqlevelsStyle(variant_ranges) <- "UCSC"

# Annotate Variants with MAF from TopMed
mafdb <- MafDb.TOPMed.freeze5.hg19
snp_pos_maf <- gscores(mafdb, variant_ranges)

# Add AF Values to the Variant Data 
variant_data$AF <- snp_pos_maf$AF

# Filter for Rare Variants (MAF < 0.01) 
rare_variants <- variant_data[variant_data$AF < 0.01 & !is.na(variant_data$AF), ]
dim(rare_variants)
# output: 57 rare variants

# Save the Annotated Variants to a TSV File
write.table(rare_variants, "variants_rare.tsv", sep = "\t", quote = FALSE, row.names = FALSE)



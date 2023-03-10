# DESCRIPTION: Checking the data, exploring, rough analyses
# AUTHOR: Julia Romanowska
# DATE CREATED: 2022-10-18
# DATE MODIFIED:

# SETUP ----
library(karyoploteR)
library(here)
library(ensemblfetchR)
library(biomaRt)
library(tidyverse)

regul_regs_file <- here("DATA", "all_regul_regs.rds")
genes_regs_file <- here("DATA", "all_genes_regs.rds")

# READ DATA ----
cpgs_stable <- read_csv2(here("DATA", "3_stabsel_cpgs_anno.csv"))

cpgs_stable_tidy <- cpgs_stable %>%
	rename(min_clock = '6_cpg_clock') %>%
	select(cpg:min_clock, CHR, MAPINFO) %>%
	mutate(
		CHR = as.factor(CHR)
	)
cpgs_stable_tidy

cpgs_stable_tidy_GRanges <- makeGRangesFromDataFrame(
		as.data.frame(cpgs_stable_tidy) %>%
			select(seqnames = CHR, start = MAPINFO, everything()) %>%
			mutate(
				end = start + 1,
				seqnames = paste0("chr", seqnames)
			),
		keep.extra.columns = TRUE
	)
cpgs_stable_tidy_GRanges

# EXPLORE ----
skimr::skim(cpgs_stable_tidy)

# statistics
cpgs_stable_tidy %>% janitor::tabyl(CHR) %>%
	arrange(percent)

# FETCH ANNOTATION ----
check_cpgs_in_region <- function(chromosome_name, chromosome_start, chromosome_end){
	cur_cpgs_stable <- cpgs_stable_tidy %>%
		filter(CHR == chromosome_name) %>%
		filter(MAPINFO >= chromosome_start &
					 	MAPINFO <= chromosome_end)
	if(nrow(cur_cpgs_stable) == 0){
		return(NA)
	}
	return(list(cur_cpgs_stable$cpg))
}


## regulatory regions ----
if(!file.exists(regul_regs_file)){
	ensembl_regul <- useEnsembl("regulation",
	    dataset = "hsapiens_regulatory_feature",
	    GRCh = 37)
	my_attribs <- c(
	  "regulatory_stable_id",
	  "chromosome_name",
	  "chromosome_start",
	  "chromosome_end",
	  "feature_type_name",
	  "feature_type_description",
	  "so_accession",
	  "activity",
	  "epigenome_name",
	  "epigenome_description"
	)
	my_filters <- c(
	  "chromosome_name",
	  "start", "end"
	)
	
	data_in <- cpgs_stable_tidy %>%
		select(CHR, MAPINFO) %>%
		mutate(
			end = MAPINFO + 1,
			CHR = as.integer(levels(CHR)[CHR])
		)
	names(data_in) <- my_filters
	data_in
	
	regul_regs <- grabRegulRegions(
		positions = data_in,
		mart = ensembl_regul,
		attribs = my_attribs,
		filters = my_filters,
		asGRanges = FALSE
	)
	regul_regs <- as_tibble(regul_regs)
	
	# match regions with CpGs
	regul_regs <- regul_regs %>%
		rowwise() %>%
		mutate(
			cpg_id = check_cpgs_in_region(chromosome_name, chromosome_start, chromosome_end)
		) %>%
		unnest(cols = cpg_id)
	
	saveRDS(regul_regs, regul_regs_file)
} else {
	regul_regs <- readRDS(regul_regs_file)
}

regul_regs %>% glimpse()

# let's narrow down to the 'epigenome' that is within cord blood
regul_regs_all <- regul_regs
regul_regs <- regul_regs_all %>%
	filter(str_detect(epigenome_name, "CB"))
regul_regs %>% glimpse()

knitr::kable(
	regul_regs %>%
		distinct(epigenome_name, epigenome_description)
)

regul_regs %>%
	distinct(regulatory_stable_id, cpg_id) %>%
	count(cpg_id, sort = TRUE)
regul_regs %>%
	distinct(regulatory_stable_id, cpg_id, feature_type_name) %>%
	count(feature_type_name, sort = TRUE)

# what about activity in various cells?
regul_regs %>%
	distinct(regulatory_stable_id, activity, epigenome_name,
					 epigenome_description, cpg_id) %>%
	janitor::tabyl(activity, epigenome_name)

# which are active?
regul_regs %>%
	filter(activity == "ACTIVE") %>%
	distinct(regulatory_stable_id, epigenome_name, epigenome_description, cpg_id,
					 feature_type_name, activity)

# only the 'minimal clock' probes
regul_regs_min_clock <- regul_regs %>%
		right_join(cpgs_stable_tidy %>% filter(min_clock), by = c("cpg_id" = "cpg"))
DT::datatable(regul_regs_min_clock)

regul_regs_min_clock %>%
	distinct(regulatory_stable_id, feature_type_name, cpg_id)

## genes ----
if(!file.exists(genes_regs_file)){
	my_attribs <- c(
	  "ensembl_gene_id",
	  "chromosome_name",
	  "start_position",
	  "end_position",
	  "description",
	  "external_gene_name"
	)
	
	genes_regs <- grabGenesEnsembl(
		positions = data_in,
		attribs = my_attribs,
		filters = my_filters,
		asGRanges = FALSE,
		genome_ver = 37
	)
	genes_regs <- as_tibble(genes_regs)

	genes_regs <- genes_regs %>%
		rowwise() %>%
		mutate(
			cpg_id = check_cpgs_in_region(chromosome_name, start_position, end_position)
		) %>%
		unnest(cols = cpg_id)

	saveRDS(genes_regs, genes_regs_file)
} else {
	genes_regs <- readRDS(genes_regs_file)
}

genes_regs
genes_regs %>% count(cpg_id, sort = TRUE)
genes_regs %>%
	count(ensembl_gene_id, sort = TRUE) %>%
	left_join(genes_regs %>% distinct(ensembl_gene_id, external_gene_name))

# cat(genes_regs %>% distinct(external_gene_name) %>% pull(), sep = "\n")

## plot positions ----
regul_regs_GRanges <- makeGRangesFromDataFrame(
		as.data.frame(
			regul_regs %>%
				select(regulatory_stable_id, starts_with("chrom"),
								 starts_with("feature"), so_accession, cpg_id) %>%
				distinct()
			) %>%
			select(
				seqnames = chromosome_name,
				start = chromosome_start,
				end = chromosome_end,
				everything()) %>%
			mutate(
				seqnames = paste0("chr", seqnames)
			),
		keep.extra.columns = TRUE
	)

genes_regs_GRanges <- makeGRangesFromDataFrame(
		as.data.frame(genes_regs) %>%
			select(
				seqnames = chromosome_name,
				start = start_position,
				end = end_position,
				everything()) %>%
			mutate(
				seqnames = paste0("chr", seqnames)
			),
		keep.extra.columns = TRUE
	)

png(
	filename = here("FIGURES", "cpgs_genomic_annotation.png"),
	width = 8,
	height = 12,
	units = "in",
	res = 150
)

pp <- getDefaultPlotParams(plot.type = 2)
pp$data1height <- 300
pp$data2height <- 600
kp <- plotKaryotype(
	genome = "hg19",
	chromosomes = "autosomal",
	plot.type = 2,
	plot.params = pp
)
kpPlotMarkers(
	kp,
	data = cpgs_stable_tidy_GRanges,
	labels = cpgs_stable_tidy$cpg,
	text.orientation = "horizontal",
	r1 = 0.5,
	cex = 0.5,
	# ignore.chromosome.ends = TRUE,
	label.color = "orange",
	line.color = "orange",
	marker.parts = c(0.4, 0.5, 0.1),
	label.dist = 0.00001,
	max.iter = 10000
)
# mark those CpGs that are in the minimal set
kpPlotMarkers(
	kp,
	data = cpgs_stable_tidy_GRanges[cpgs_stable_tidy$min_clock,],
	labels = cpgs_stable_tidy %>%
		filter(min_clock) %>% pull(cpg),
	text.orientation = "horizontal",
	r1 = 0.5,
	cex = 0.5,
	# ignore.chromosome.ends = TRUE,
	label.color = "red",
	line.color = "red",
	marker.parts = c(0.4, 0.5, 0.1),
	label.dist = 0.00001,
	max.iter = 10000
)
# plot background for regulatory regions and genes
kpDataBackground(
	kp,
	r0 = 0.5,
	r1 = 1.1,
	data.panel = 2,
	color = "gray80"
)
kpDataBackground(
	kp,
	r0 = 0,
	r1 = 0.6,
	data.panel = 2,
	color = "#DEFFD4"
)
# plot regulatory regions
kpPlotMarkers(
	kp,
	data = regul_regs_GRanges,
	labels = regul_regs_GRanges$feature_type_name,
	text.orientation = "horizontal",
	data.panel = 2,
	r1 = 0.3,
	cex = 0.5,
	# ignore.chromosome.ends = TRUE,
	label.color = "darkgreen",
	line.color = "darkgreen",
	marker.parts = c(0.4, 0.5, 0.1),
	label.dist = 0.00001,
	max.iter = 10000
)
# plot regulatory regions
kpPlotMarkers(
	kp,
	data = genes_regs_GRanges,
	labels = genes_regs_GRanges$external_gene_name,
	text.orientation = "horizontal",
	data.panel = 2,
	r1 = 1,
	cex = 0.5,
	# ignore.chromosome.ends = TRUE,
	label.color = "black",
	line.color = "black",
	marker.parts = c(0.4, 0.5, 0.1),
	label.dist = 0.00001,
	max.iter = 10000
)

dev.off()

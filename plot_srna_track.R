#!/usr/bin/env Rscript --vanilla
options(scipen = 999)

# require platform "unix". Sorry windows !
if(.Platform$OS.type != "unix") {
  stop('.Platform$OS.type of your device must be "unix"!')
}

# Load required libraries, quietly
if(!suppressPackageStartupMessages(require(data.table))) {
  stop("R package data.table is required, but could not be loaded. Make sure it is installed!")
}
if(!suppressPackageStartupMessages(require(docopt))) {
  stop("R package docopt is required, but could not be loaded. Make sure it is installed!")
}
if(!suppressPackageStartupMessages(require(tidyverse))) {
  stop("R package tidyverse is required, but could not be loaded. Make sure it is installed!")
}
if(!suppressPackageStartupMessages(require(cowplot))) {
  stop("R package cowplot is required, but could not be loaded. Make sure it is installed!")
}
if(!suppressPackageStartupMessages(require(IRanges))) {
  stop("R package IRanges is required, but could not be loaded. Make sure it is installed!")
}


# parse coords
pdf_file <- "test.pdf"
coords <- 'RNA1:260-310'
parsed_coords <- str_match(coords, '([^:]+):(\\d+)-(\\d+)')
chrom <- parsed_coords[1,2]
userStart <- as.numeric(parsed_coords[1,3])
userEnd <- as.numeric(parsed_coords[1,4])

# coords size limit is 100,000
if((userEnd - userStart) > 100000) {
  stop("Maximum interval size is 100,000 nts. Please revise your coordinates.")
}

filelist <- list.files("5_map2cmv/",".*-1.*.csv")

sRNAData <- data.frame()

for (file in filelist) {
  tmp <- fread(paste0("5_map2cmv/", file))
  tmp <- tmp %>% filter(Chromosome == chrom, Position %in% c(userStart:userEnd))
  tmp$RPM[tmp$Strand == "-"] <- tmp$RPM[tmp$Strand == "-"] * -1
  sRNAData <- rbind(sRNAData, tmp)
}

# Transformation to tidy format suitable for plotting
# sRNAData <- pivot_longer(sRNAData, cols=3:8, names_to = "RNA Size", values_to = "RPM")

# define color palette
sRNAcols = c("lightgray", # <21nts
             "blue",      # 21nts
             "mediumseagreen", # 22nts
             "orange", # 23 nts
             "tomato",    # 24nts
             "darkgray")  # >24nts

# define the order to list the Length categories
Lorder = c("<21","21","22","23", "24",">24")

# Get plus and minus, and remove useless 0 rows. That will save time during plotting
# Must keep first and last positions.

PosMin = min(sRNAData$Position)
PosMax = max(sRNAData$Position)

cov_plus <- filter(sRNAData, Strand == '+', (RPM != 0 | Position %in% c(PosMin, PosMax)))
cov_minus <- filter(sRNAData, Strand == '-', (RPM != 0 | Position %in% c(PosMin, PosMax)))


psRNA <- ggplot() +
  geom_col(data = cov_plus, 
           aes(x = Position, y = RPM, 
               fill = factor(`RNA Size`, levels=Lorder)),
           width = 1) +
  geom_col(data = cov_minus, 
           aes(x = Position, y = RPM, 
               fill = factor(`RNA Size`, levels=Lorder)),
           width = 1) +
  scale_fill_manual(values = sRNAcols, name = "RNA length") +
  facet_grid(bamName ~ .) +
  geom_hline(yintercept = 0) +
  theme_classic() +
  labs(x = coords) +
  coord_cartesian(xlim = c(userStart, userEnd))

# Save the plot
# Sizing of output plot depends on number of bams OR readgroups and depth of mRNA packing
# vInches = nrow(bams)

RGlist <- NULL
if(is.null(RGlist)) {
  # seqRows = nrow(bams)
  seqRows = length(filelist)
} else {
  seqRows = nrow(rgs)
}


vInches <- seqRows
finalP <- psRNA

ggsave2(pdf_file, plot = finalP, width = 7, height = vInches, units = "in")

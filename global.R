library(stringr)
library(reshape2)
library(RColorBrewer)
library(shiny)
library(shinyWidgets)
library(DT)
library(plotly)
library(omicsViewer)
library(heatmaply)
library(randomcoloR)

source("00_util.R")
dff <- readRDS("finalTable_filter_unipep1_ef06.Rds")
dff$bp.domain.PID[dff$bp.gene.name == "DOK2"] <- "+"
ii <- nchar(dff$pY.sequence) != "12" | substr(dff$pY.sequence, 6, 7) != "pY"
unique(dff$pY.sequence[ii])
smap <- c(
  "NNYQFpY" = "NNYQFpY_____",
  "LQPNNpYQFC" = "LQPNNpYQFC__",
  "RRRpYIVRKR" = "__RRRpYIVRKR",
  "TLNLCpYSA" = "TLNLCpYSA___",
  "SNDCSpYGG" = "SNDCSpYGG___",
  "GPAPQpY" = "GPAPQpY_____",
  "PYSPCpYPDPR" = "PYSPCpYPDPR_",
  "GSKKIpYDI" = "GSKKIpYDI___",
  "FTDNSpY" = "FTDNSpY_____",
  "NpYVSCCK" = "____NpYVSCCK",
  "LWNPTpYRS" = "LWNPTpYRS___",
  "KLMDTpYDS" = "KLMDTpYDS___",
  "AALGApYV" = "AALGApYV____",
  "SLFYNpYSML" = "SLFYNpYSML__")
for (i in names(smap))
  dff$pY.sequence[dff$pY.sequence == i] <- smap[i]
dff$pY.sequence <- sub("p", "", dff$pY.sequence)

dat <- dff
rtk <- unique(dat$rtk)
pySites <- lapply(split(dat$pY.id2, dat$rtk), unique)

ii <- dat$bp.domain.PID == "+" | dat$bp.domain.SH2 == "+"
dat_domain <- dat[ii, ]
pl_domain <- unique(dat_domain$bp.gene.name)
bp <- unique(c(pl_domain, dat$bp.gene.name))

pi <- readRDS("shiny_RTKuniprot_featureExtract_expSum.RDS")
names(pi) <- sapply(pi, function(x) as.character(x$summary$proteinName[1]))

motif_bg <- omicsViewer:::aaFreq(unique(dat$pY.sequence))
tt <- "The phosphotyrosine interactome landscape of human receptor tyrosine kinases"

scc <- c(
  "Interacting gene name" = "bp.gene.name", 
  "Enrichment Score" = "af.filter.batch",
  "Phosphotyrosine site" = "pYs",
  "Phosphotyrosine sequence" = "pY.sequence",
  "Type pY" = "pY.type",
  "PTB domain" = "bp.domain.PID",
  "SH2 domain" = "bp.domain.SH2"
)

scc2 <- c(
  "Interacting gene name" = "bp.gene.name", 
  "Enrichment Score" = "af.filter.batch",
  "Phosphotyrosine site" = "pY.gene.name",
  "Phosphotyrosine sequence" = "pY.sequence",
  "Type pY" = "pY.type",
  "PTB domain" = "bp.domain.PID",
  "SH2 domain" = "bp.domain.SH2"
)

source("functions.R")
source("L1_module_pY.R")
source("L1_module_bindingProtein.R")
source("L1_module_heatmap.R")
source("L2_module_motif.R")



source("Functions.R")
library(pacman)
print(Sys.time())

pacman::p_load(eulerr, ggalt, ggbeeswarm, OrgMassSpecR, robustbase, plyr,
               stringi, stringr, tidyr, patchwork, ggplot2, reshape2,
               Cairo, gridGraphics, data.table, devtools, diann, gridExtra, ggpubr,
               ggpointdensity, viridis, scales, matrixStats, dplyr, targets,
               MASS, flux, psych, sva, ggExtra, gprofiler2, NCmisc, ggrepel, ggridges, rstatix, dtplyr,
               lawstat,lme4, lsmeans, boot,msEmpiRe,rsvd,enviPat,qs,PupillometryR,gtools,stats, plotly,
               htmlwidgets,htmltools,gghalves,bayestestR,distributions)

############ SINGLE
SC_fpath <- "/work/slavov/JD/DATA_SN"
SC_fpath_Ext <- "/work/slavov/JD/DATA_SN_extracted"

############ BULK
bulk_nuc_fpaths_tims <- "/work/slavov/JD/DATA_bulk_nuc" #BR 1, 2, 3, 4, 5, 6
bulk_WC_fpaths_tims <- "/work/slavov/JD/DATA_bulk_WC" #BR 1, 2, 3, 4, 5, 6

############ BULK PURITY
bulk_purity_fpath_timsTOFSCP <- "/work/slavov/JD/DATA_bulk_purity" 

############ siRNA
bulk_nuc_siRNA_all_fpath <- "/work/slavov/JD/DATA_siRNA_nuc_new"
bulk_WC_siRNA_all_fpath <- "/work/slavov/JD/DATA_siRNA_WC_new"

############ META
meta_path <- "/work/slavov/JD/META/Meta_SC_SN_04122024_Prep1_2_3_6_final.tsv"
meta_bulk_path_timsTOFSCP <- "/work/slavov/JD/META/Meta_THP1Bulk_01172024_timsTOFSCP.tsv"
meta_bulk_path_timsTOFSCP_siRNA <- "/work/slavov/JD/META/Meta_THP1Bulk_siRNA_01222024_v3.tsv"

############ Human yeast
HY_fpath <- "/work/slavov/JD/DATA_HY"


############ MISC:
masses_fpath <- "/work/slavov/JD/Misc_files/uniprot_masses.tsv" # Protein mass info
SP <- "/work/slavov/JD/Misc_files/swissprot_genes_10072023.tsv" # Swissprot info
HL_path <- "/work/slavov/JD/Misc_files/SILAC_protein_turnOver/protein half lives high qual-Table 1.tsv" # Half life data
markers_fpath <- "/work/slavov/JD/Misc_files/Subcellular_protein_07102022.tsv" # Subcellular protein markers
NPC_regions <- "/work/slavov/JD/Misc_files/NUP_regions.tsv"
NPC_scaffoldProts <- "/work/slavov/JD/Misc_files/NPC_scaffoldProts_Mathieson2018.tsv"
all_TFs <- "/work/slavov/JD/Misc_files/TF_names_v_1.01.txt"
all_TFs <- read.delim(all_TFs, header=F)

DA_TFs_to_label <- c("JUN","JUNB","FOS","NFKB1","NFKB2","REL","RELA","ATF3")
NUPs_to_label <- c("TPR")

all_TFs <- all_TFs[all_TFs$V1%!in%DA_TFs_to_label,]


nucleoporins <- c("NUP160","NUP85","SEH1L","NUP96","SEC13","NUP107","NUP133","NUP43","NUP37","ELYS","AHCTF1",
                  "NUP205","NUP188","NUP93","NUP155","NUP53","NUP54","NUP58","NUP62","NUP98","NDC1","NUP210","POM121C","POM121","RAE1","NUP42",
                  "NUP88","NUP214","DDX19","GLE1","NUP358","NUP153","NUP50","TPR","AAAS")


for_NPC_annot <- c("NUP205","NUP62","NUP210","NUP50","TPR","NUP160")

enhancers <- c("NOP56","MYBBP1A","BCLAF1","DHX15","RRS1","NUP205","BOP1","RBM6","MDN1")
dampeners <- c("PLEC","SND1","IMMT","HSD17B4","RPSA","VIM","TPR")
keep_KDs <- c("DHX15_kd","BOP1_kd","RBM6_kd","MDN1_kd","PLEC_kd","IMMT_kd","VIM_kd",
              "TPR_kd","SND1_kd","HSD17B4_kd","RPSA_kd","NOP56_kd","MYBBP1A_kd",
              "BCLAF1_kd","RRS1_kd","NUP205_kd")

keep_PSEA_new <- c("GO:0034063", "GO:0072673", "GO:0006810", "GO:0016973",
                   "GO:0140223", "GO:0042254", "GO:0002183", "GO:0000054", "GO:0044391", "GO:0045824", 
                   "GO:0003729", "GO:0017056", "GO:0006352","GO:0032496", 
                   "GO:0005643","GO:0015031","GO:0006338","GO:0000502","GO:0006367","GO:0071014","GO:0006954",
                   "GO:0051020","GO:0061608","GO:0090685","GO:0006412","GO:0033673","GO:0006914","GO:0008303",
                   "GO:0035091","GO:1902807","GO:0034440", "GO:0030695", "GO:0006952","GO:0008379", "GO:0051123",
                   "GO:0042254", "GO:0022626", "GO:0043039","GO:0003746","GO:0003743","GO:0001682","GO:0050685")

keep_PSEA_WC <- c("GO:0006913","GO:0006412", "GO:0005759", "GO:0044391", "GO:0042254","GO:0005681", "GO:0006091","GO:0030127",
                  "GO:0032543","GO:0006163","GO:0099503","GO:0006950","GO:0006913","GO:0009060","GO:0046390","GO:0051236","GO:0005832",
                  "GO:0051168","GO:0008152","GO:0019646", "GO:0016477","GO:0070847","GO:0046488","GO:00324960","GO:0106310","GO:0017056")

PSEA_GO_terms_siRNAWC <- c("GO:0008305","GO:0004197","GO:0005967","GO:0051208","GO:0050801","GO:0005764","GO:0008610","GO:0016049","GO:0051059","GO:0009617",
                           "GO:0000785","GO:0046483","GO:0000278","GO:0006281","GO:0006911","GO:0033554","GO:0006611","GO:0042254","GO:0006913",
                           "GO:0046324","GO:1902107","GO:0045087","GO:0043122","GO:0030225","GO:0051170")

DNA_repair_complex <- c("ASCC1", "ASCC2", "ASCC3", "ATM", "FANCD2", "FANCI", "KDM1A", "LIG4", "PRKDC", "RCOR1", "TP53BP1", "XRCC4", "XRCC5", "XRCC6", "PAXX")


translation_cyt_action <- c("GO:0022626","GO:0003743","GO:0017101")
translation_nuc_action <- c("GO:0042254", "GO:0006397")

poorly_characterized <- c("CCDC12", "TRIQK", "HNRNPUL2", "URB1")

SN_GOterms_PSEA <- c("GO:0042254","GO:0006364","GO:0002181","GO:0006913","GO:0051292",
                     "GO:0045023","GO:0022625","GO:0022627","GO:0005643","GO:0005681",
                     "GO:0140142","GO:0017116","GO:0043603")

SN_GOterms_PSEA_06032024 <- c("GO:0016072", "GO:0005681", "GO:0002181", "GO:0005925", 
                          "GO:0051301", "GO:0045088", "GO:0006886", "GO:0046931", "GO:0005682", "GO:00081870", 
                          "GO:0043039", "GO:0042555")

transport_dyn <- c("GO:0051169", "GO:0051168","GO:0005643", "GO:0051292", "GO:0044614",
                   "GO:0006913", "GO:0046907", "GO:0046822","GO:0032386","GO:0006886",
                   "GO:0031074","KEGG:03013","GO:0044613","GO:0065002","REAC:R-HSA-3301854",
                   "GO:0005635","GO:0005643","GO:0006607","GO:0006999","GO:0031080","GO:0031965",
                   "GO:0042306","GO:0042307","GO:0042564","GO:0044613","GO:0044614","GO:0044615",
                   "GO:0046824","GO:0051664","GO:0090435","GO:0140142","REAC:R-HSA-3301854")

GOs_siRNA <-  c("GO:1990745","GO:0045842","GO:0071159","GO:0071363","GO:0032770","GO:0032722","GO:0032620",
                "GO:0002269","GO:2000192", "GO:0030225","GO:0034142","GO:0050729","GO:0060907","GO:0050999","GO:0030687","GO:0000028","GO:0030688",
                "GO:0051168","GO:0005049","GO:0000054","GO:0000445","GO:0051170","GO:0061676","GO:0006607","GO:0042564","GO:0061608",
                "GO:0042307","GO:0006606","GO:0032496","GO:0090382","GO:0043162")

GOs_dynamics <- c("GO:0006607","GO:0051020", "GO:0000054", "GO:0031063", "GO:0045089", 
                  "GO:0006886", "GO:0051301", "GO:0036402", "GO:0051123", "GO:0005643", "GO:0038061", "GO:0034599",
                  "GO:0006606","GO:0051168","GO:0008134","GO:0000278")

# Vectors for each category
nuclear_transport_and_ribosome_biogenesis = c(
  "GO:1990745", # EARP complex
  "GO:0030687", # preribosome, large subunit precursor
  "GO:0000028", # ribosomal small subunit assembly
  "GO:0030688", # preribosome, small subunit precursor
  "GO:0051168", # nuclear export
  "GO:0005049", # nuclear export signal receptor activity
  "GO:0000054", # ribosomal subunit export from nucleus
  "GO:0000445", # THO complex part of transcription export complex
  "GO:0061676", # importin-alpha family protein binding
  "GO:0006607", # NLS-bearing protein import into nucleus
  "GO:0042564", # NLS-dependent protein nuclear import complex
  "GO:0061608", # nuclear import signal receptor activity
  "GO:0042307", # positive regulation of protein import into nucleus
  "GO:0006606"  # protein import into nucleus
)

inflammatory_response_and_immune_activation = c(
  "GO:0034142", # toll-like receptor 4 signaling pathway
  "GO:0050729", # positive regulation of inflammatory response
  "GO:0060907", # positive regulation of macrophage cytokine production
  "GO:0002269", # leukocyte activation involved in inflammatory response
  "GO:0032620", # interleukin-17 production
  "GO:0032722", # positive regulation of chemokine production
  "GO:0032496",  # response to lipopolysaccharide
  "GO:0071159" # NF-kappaB complex
  
)

macrophage_function = c(
  "GO:0032770", # positive regulation of monooxygenase activity
  "GO:0050999", # regulation of nitric-oxide synthase activity
  "GO:0030225", # macrophage differentiation
  "GO:0090382", # phagosome maturation
  "GO:2000192" # negative regulation of fatty acid transport
)

protein_modification_and_stability = c(
  "GO:0051438", # regulation of ubiquitin-protein transferase activity
  "GO:0043162", # ubiquitin-dependent protein catabolic process via the multivesicular body sorting pathway
  "GO:0035173", # histone kinase activity
  "GO:0008168", # methyltransferase activity
  "GO:0031648"  # protein destabilization
)

cell_cycle_regulation_and_growth = c(
  "GO:0045842", # positive regulation of mitotic metaphase/anaphase transition
  "GO:0051301",  # cell division
  "GO:0030307",  # positive regulation of cell growth
  "GO:0071363" # cellular response to vascular endothelial growth factor stimulus
)



source_targets <- function(path) {
  source(path)$value
}

targets_list <- c(
  source_targets("BulkSamples_targets.R"),
  source_targets("SCoPEDIA_benchmarking_targets.R"),
  source_targets("SingleNuclei_targets.R"),
  source_targets("WriteTables_targets.R"),
  source_targets("KnockDowns_targets.R"),
  source_targets("SavePlots_targets.R")
)

targets_list <- unlist(targets_list, recursive = FALSE)

targets_list


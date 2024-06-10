
# Here we analyze the knock-down data

list(
  #######################################################
  ############# Bulk siRNA KD analyses:  ################
  #######################################################
  
  # siRNA NUC and WCs, get data
  tar_target(data_nuc_bulk_siRNA, pD_get_data(bulk_nuc_siRNA_all_fpath, meta_bulk_path_timsTOFSCP_siRNA, filter=F, proteotypic=T, ppm=5, mono_only=T, correct_IsoEnv=T),format = "qs"),
  tar_target(data_wc_bulk_siRNA, pD_get_data(bulk_WC_siRNA_all_fpath, meta_bulk_path_timsTOFSCP_siRNA, filter=F, proteotypic=T, ppm=5, mono_only=T, correct_IsoEnv=T),format = "qs"),

  # siRNA NUC: Investigate changes in protein transport
  tar_target(data_nuc_bulk_siRNA_fil, rm_data(data_nuc_bulk_siRNA, rm=c("gJD4645_B-B12_1_6214","gJD4646_B-B13_1_6215","gJD4644_B-B11_1_6210")),format = "qs"),  #these files are mixed in with WC-sets and produced a large batch effect, just remove them.
  tar_target(nuc_maxLFQ_siRNA, diann_maxlfq(data_nuc_bulk_siRNA_fil, group.header="Genes", id.header = "seqcharge", quantity.header = "Ms1.Area_iso"),format = "qs"), 
  tar_target(nuc_maxLFQ_siRNA_m, pD_melt_MaxLFQ(nuc_maxLFQ_siRNA, meta_bulk_path_timsTOFSCP_siRNA),format = "qs"), 
  tar_target(nuc_imp_siRNA, pD_impute_bulk(nuc_maxLFQ_siRNA_m, missing.prot.frac=0.95, missing.cell.frac=0.95, k=3, Label=T, BioRep=T, LC_batch=T),format = "qs"), 
  tar_target(nuc_BC_bulk_siRNA, pD_bulk_BC(nuc_imp_siRNA, meta_bulk_path_timsTOFSCP_siRNA, Label=T, BioRep=F, LC_batch=F, Cell_Batch=F,use_mod=F),format = "qs"),
  tar_target(Effect_on_transport_siRNA_weighted_25, pD_effect_on_transport(nuc_BC_bulk_siRNA, meta_bulk_path_timsTOFSCP_siRNA, MsEmpire_bulk_BR123456_merged, Transport_Score_cor, titlename="top25", mostdif=25, useall=F),format = "qs"), #use top 25 most differential proteins
  tar_target(Effect_on_transport_siRNA_weighted_100, pD_effect_on_transport(nuc_BC_bulk_siRNA, meta_bulk_path_timsTOFSCP_siRNA, MsEmpire_bulk_BR123456_merged, Transport_Score_cor, titlename="top100", mostdif=100, useall=F),format = "qs"), #use top 100 most differential proteins
  tar_target(Effect_on_transport_siRNA_weighted_qval01, pD_effect_on_transport(nuc_BC_bulk_siRNA, meta_bulk_path_timsTOFSCP_siRNA, MsEmpire_bulk_BR123456_merged, Transport_Score_cor, titlename="qval01", mostdif=50, useall=T),format = "qs"), #use 1% FDR differential proteins
  tar_target(cors_DampEnhance, pD_plot_cors_DE(Transport_Score_cor$cors, Enhancers=enhancers,Dampeners=dampeners),format = "qs"),
  
  # checking that the proteins were actually KD'ed in nucleus
  tar_target(bulk_DA_siRNA_NUC, pD_bulk_DA_all(nuc_BC_bulk_siRNA, collapse_conds=F),format = "qs"),
  tar_target(siRNA_NUC_all, pD_plot_DA_v2(bulk_DA_siRNA_NUC, xlab="Log2, KD / Neg KD", titleLab = "Nuclei"),format = "qs"),
  
  
  # siRNA WC:
  tar_target(WC_maxLFQ_siRNA, diann_maxlfq(data_wc_bulk_siRNA, group.header="Genes", id.header = "seqcharge", quantity.header = "Ms1.Area_iso"),format = "qs"),
  tar_target(WC_maxLFQ_siRNA_m, pD_melt_MaxLFQ(WC_maxLFQ_siRNA, meta_bulk_path_timsTOFSCP_siRNA),format = "qs"),
  tar_target(WC_imp_siRNA, pD_impute_bulk(WC_maxLFQ_siRNA_m, missing.prot.frac=0.95, missing.cell.frac=0.95, k=3, Label=T, BioRep=T, LC_batch=T),format = "qs"),
  tar_target(WC_BC_bulk_siRNA, pD_bulk_BC(WC_imp_siRNA, meta_bulk_path_timsTOFSCP_siRNA, Label=T, BioRep=F, LC_batch=F, Cell_Batch=F,use_mod=F),format = "qs"),
  
  # checking that the proteins were actually KD'ed in whole-cell abundance
  tar_target(bulk_DA_siRNA_WC, pD_bulk_DA_all(WC_BC_bulk_siRNA, collapse_conds=F),format = "qs"),
  tar_target(siRNA_WC_all, pD_plot_DA_v2(bulk_DA_siRNA_WC, xlab="Log2, KD / Neg KD", titleLab = "Whole-cells"),format = "qs"),
  
  # are the whole-cell proteomes indicative of cell-states that are primed to mediate LPS-induced nucleocytoplasmic transport?
  tar_target(PSEA_KDs_WC, pD_PSEA_bulk_siRNA(bulk_DA_siRNA_WC,keep_KDs,min_prot=3,min_fraction = 0.05,ontology_sources =c("GO")),format = "qs"),
  tar_target(PSEA_siRNA_KDs_WC_corOrder, pD_plot_PSEA_siRNA_bulk_cors_order(PSEA_KDs_WC, bulk_DA_siRNA_WC, Effect_on_transport_siRNA_weighted_100, PSEA_GO_terms_siRNAWC, Transport_Score_cor, KD_genes = c(dampeners,enhancers)),format = "qs")
  
)

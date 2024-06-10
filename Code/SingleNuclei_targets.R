
# Here we analyze the single nuclei data
# Note: "SN" = "Single Nuclei"

list(
  
  #######################################################
  #####  Single nucleus quant and quality control:  #####
  #######################################################
  
  # for filtering out poorly prepared/acquired single nuclei
  tar_target(data, pD_get_data(SC_fpath,meta_path,filter=F, proteotypic=F, correct_IsoEnv=F),format = "qs"), 
  tar_target(keep_cells, pD_filter_count(data,TQV_threshold=0.1,CQV_threshold=0.1,min_count=50),format = "qs"), 
  tar_target(data_proteotypic, pD_filter_proteotypic(data),format = "qs"), 
  tar_target(med_car, pD_median_car(data_proteotypic, car="mTRAQ8"),format = "qs"), 

  tar_target(data_good_cells, pD_filterCells(data_proteotypic,keep_cells),format = "qs"), 

  # for filtering out precursors with systematically poor quant
  tar_target(norm, pD_car_norm_data(data_good_cells, Run ="Run", seqcharge ="seqcharge", Celltype="Celltype", Quant="Ms1.Area"),format = "qs"), 
  tar_target(prec_correct, pD_rm_poorQuant(norm, min_proportion=0.01, compression_factor=6),format = "qs"), 
  
  # MaxLFQ
  tar_target(SN_maxLFQprep, pD_SN_maxLFQPrep(norm, med_car),format = "qs"), 
  tar_target(SN_maxLFQprep_select, pD_intersect_seqcharge(SN_maxLFQprep, prec_correct),format = "qs"), 
  tar_target(SN_correct_forLC, pD_spline_correct_forLC(SN_maxLFQprep_select),format = "qs"),
  tar_target(SN_maxLFQ, diann_maxlfq(SN_correct_forLC, group.header="Genes", id.header = "seqcharge", quantity.header = "adj_val"),format = "qs"), 
  tar_target(dat_SN_maxLFQ_m, pD_melt_MaxLFQ(SN_maxLFQ, meta_path),format = "qs"), 
  
  # assess purity of single nuclei and remove ones with poor nuclear enrichment
  tar_target(Align_nucs_MLFQ, pD_nuc_align(dat_SN_maxLFQ_m, dat_bulk_purity),format = "qs"), 
  tar_target(AUC_purity_MLFQ, pD_AUC_nucPurity(Align_nucs_MLFQ, maxAUC=5),format = "qs"), 
  tar_target(nucs_MLFQ, pD_filter_purity(dat_SN_maxLFQ_m, AUC_purity_MLFQ, maxAUC = 5),format = "qs"), 
  
  # normalize, impute, batch correct
  tar_target(nucs_norm_MLFQ, pD_prot_norm(nucs_MLFQ, Quant="norm_prot"),format = "qs"), 
  tar_target(nucs_imputed_MLFQ_log2, pD_impute_log2(nucs_norm_MLFQ, missing.prot.frac = 0.95, missing.cell.frac = 0.95, k=5),format = "qs"), 
  tar_target(nucs_BC_MLFQ, pD_batchCorrect_labs_LC(nucs_imputed_MLFQ_log2, meta_path),format = "qs"), 
  tar_target(pre_post_filter, pD_prepost(data_proteotypic, nucs_BC_MLFQ$unimputed.BC),format = "qs"), 
  tar_target(Overlap_prot_ab, pD_overlap(nucs_BC_MLFQ, prot=c("NUP205","VIM")),format = "qs"), 
  
  
  #######################################################
  ######### Getting weights from bulk data:  ############
  #######################################################
  
  # Same bulk data as before, but filtering for higher quality and bulk biological replicates which correspond.... 
  # to the single nuclei (our single nuclei correspond to bio-reps 1,2,3,6 from bulk data)
  tar_target(output_MSEmpPrep_ind_0.05filt, pD_MsEmpire_prep_specificBR(data_nuc_bulk,BR=c(1,2,3,6),conds=c("NT","min_10","min_30","min_60"),TQVfilt=0.05,CQfilt=0.05),format = "qs"), 
  tar_target(MsEmpire_bulk_ind_0.05filt,  pD_MsEmpire_run_specific(output_MSEmpPrep_ind_0.05filt),format = "qs"), 
  tar_target(MsEmpire_bulk_ind_0.05filt_merged,  pD_merge_MSEmpire(MsEmpire_bulk_ind_0.05filt,quant="p.val"),format = "qs"), 
  tar_target(MsEmpire_bulk_ind_0.05filt_bartel, pD_BartelTest(MsEmpire_bulk_ind_0.05filt_merged),format = "qs"), 
  
  
  #######################################################
  ############## Biological analyses:  ##################
  #######################################################
  
  # PCA
  tar_target(PCA, pD_weightedBulk_PCA(nucs_BC_MLFQ, MsEmpire_bulk_ind_0.05filt_merged, MsEmpire_bulk_ind_0.05filt_bartel, quant="p.val"),format = "qs"), 

  # Transport score
  tar_target(Transport_Score, pD_TransportScore(nucs_BC_MLFQ, MsEmpire_bulk_ind_0.05filt_merged, frac=0.75, quant="p.val",FC="log2FC"),format = "qs"), 
 
  # correlation to transport score
  tar_target(Transport_Score_cor, pD_TransportScore_cor(Transport_Score, keep=c("NUC_10","NUC_30","NUC_60"),keepGenes=c("NUP205","VIM")),format = "qs"), 
  tar_target(PSEA_TS_cor, pD_PSEA_simple_v2(Transport_Score_cor$cors, min_prot=3, min_fraction = 0.05, ontology_sources =c("GO")),format = "qs"), 
  tar_target(PSEA_TS_cor_plot, pD_SN_PSEA_plot(PSEA_TS_cor, GOs=SN_GOterms_PSEA_06032024),format = "qs"), 

  # Nuclear pore complex analysis
  tar_target(Nucleoporins_TS_cor, pD_TS_cor_plot(Transport_Score_cor, NPC_regions, for_NPC_annot, NPC_scaffoldProts, genes=nucleoporins),format = "qs"), 
  tar_target(NPC_TS_cor, pD_NPC_Zscore_cor(Transport_Score, nucleoporins),format = "qs"), 
  tar_target(NPC_halfLife_cor, pD_halfLives(Transport_Score_cor, HL_path, annot = c("RANBP2", "NUP153", "NUP62", "NUP214", "NUP93", "NUP98", "NUP88", "NUP50")),format = "qs"), 
  tar_target(NPC_passiveDiffusion, pD_passive_NPC(Transport_Score, nucleoporins, SN=c("gJD4449_B-L24_1_5732mTRAQ4", "gJD2547_B-K8_1_2546mTRAQ4")),format = "qs") 
  
)



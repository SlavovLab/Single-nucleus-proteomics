
# Here we analyze the bulk nuclei and whole-cell data

list(
  #######################################################
  ################ Nuclei bulk analysis:  ###############
  #######################################################
  tar_target(data_nuc_bulk, pD_get_data(bulk_nuc_fpaths_tims, meta_bulk_path_timsTOFSCP, filter=F, proteotypic=T, ppm=5, mono_only=T,correct_IsoEnv=T),format = "qs"), 

  # BULK NUC MS-Empire
  tar_target(output_MSEmpPrep_BR123456, pD_MsEmpire_prep_specificBR(data_nuc_bulk,BR=c(1,2,3,4,5,6),conds=c("NT","min_10","min_30","min_60"),TQVfilt=1.1,CQfilt=1.1),format = "qs"), 
  tar_target(MsEmpire_bulk_BR123456,  pD_MsEmpire_run_specific(output_MSEmpPrep_BR123456),format = "qs"), 
  tar_target(MsEmpire_bulk_BR123456_merged,  pD_merge_MSEmpire(MsEmpire_bulk_BR123456,quant="p.val"),format = "qs"), 
  tar_target(MSEmpire_bulk_BR123456_merged_bartel, pD_BartelTest(MsEmpire_bulk_BR123456_merged),format = "qs"), 
  tar_target(timeSeries_TFs_nuc_DNARepair, pD_timeseries_MSEmp(MsEmpire_bulk_BR123456_merged, tofilt = DNA_repair_complex, c("FANCI","ASCC3","PAXX","PRKDC","XRCC6","XRCC5"), 
                                                               p.adj_filter=1.1,filename="DNA_repair", spec_width=3.3,spec_height=4.1, multiple=F),format = "qs"), 
  tar_target(timeSeries_TFs_nuc, pD_timeseries_MSEmp(MsEmpire_bulk_BR123456_merged, tofilt = NULL, DA_TFs_to_label, p.adj_filter=0.05,filename="TFs",
                                                     spec_width=3.3,spec_height=4.1, multiple=F),format = "qs"), 
  tar_target(timeSeries_allTFs_nuc, pD_timeseries_MSEmp(MsEmpire_bulk_BR123456_merged, tofilt = NULL, all_TFs, p.adj_filter=0.05,filename="allTFs",
                                                        spec_width=3.3,spec_height=4.1, multiple=F),format = "qs"), 

  #all_TFs
  tar_target(timeSeries_TFs_WC, pD_timeseries_MSEmp(MsEmpire_bulk_BR123456_nofilt_WC_merged, tofilt = NULL, DA_TFs_to_label, p.adj_filter=0.05,filename="TFs_WC",
                                                    spec_width=3.3,spec_height=4.1, multiple=F),format = "qs"), 
  tar_target(combine_timeSeries, pD_combine_data(dat = list(MsEmpire_bulk_BR123456_nofilt_WC_merged,MsEmpire_bulk_BR123456_merged),labels=c("WC","NUC")),format = "qs"), 
  tar_target(nuc_filteredpval, pD_pval_filt(MsEmpire_bulk_BR123456_merged,genefilter=nucleoporins,threshold=0.05),format = "qs"), 
  tar_target(timeSeries_NUPs_nucWC, pD_timeseries_MSEmp(combine_timeSeries, tofilt = nuc_filteredpval$uniprots, NUPs_to_label, p.adj_filter=1.1,filename="NUPS_nuc",
                                                        spec_width=5,spec_height=2.3, multiple=T),format = "qs"), 
  
  # BULK NUC passive diffusion
  tar_target(Nuc_massDependence, pD_passive_diffusion(MsEmpire_bulk_BR123456_merged, masses_fpath),format = "qs"), 
  tar_target(Nuc_massFC_dynamics, pD_mass_FC_time(MsEmpire_bulk_BR123456_merged, masses_fpath, cond=c("min_60")),format = "qs"), 

  # BULK NUC time series
  tar_target(PSEA_bulk_nuc_10, pD_PSEA_simple(MsEmpire_bulk_BR123456_merged, cond="min_10", min_prot=3, min_fraction = 0.05, ontology_sources =c("GO")),format = "qs"), 
  tar_target(PSEA_bulk_nuc_30, pD_PSEA_simple(MsEmpire_bulk_BR123456_merged, cond="min_30", min_prot=3, min_fraction = 0.05, ontology_sources =c("GO")),format = "qs"), 
  tar_target(PSEA_bulk_nuc_60, pD_PSEA_simple(MsEmpire_bulk_BR123456_merged, cond="min_60", min_prot=3, min_fraction = 0.05, ontology_sources =c("GO")),format = "qs"), 
  tar_target(bulk_PSEA_joined, pD_PSEA_join_bulk(PSEA_bulk_nuc_10, PSEA_bulk_nuc_30, PSEA_bulk_nuc_60, MsEmpire_bulk_BR123456_merged, Cond_col ="Cond", quant="log2FC"),format = "qs"), 
  tar_target(PSEA_bulk_nuc_plot, pD_plot_PSEA_bulk(bulk_PSEA_joined, MsEmpire_bulk_BR123456_merged, keep=keep_PSEA_new,
                                                      conds=c("min_10","min_30","min_60"),title_plot="bulk_nuc", spec_width=6.6,spec_height=8.5),format = "qs"), 
  
  #######################################################
  ################ WC bulk analysis:  ###################
  #######################################################
  tar_target(data_wc_bulk, pD_get_data(bulk_WC_fpaths_tims, meta_bulk_path_timsTOFSCP, filter=F, proteotypic=T, ppm=5, mono_only=T, correct_IsoEnv=T),format = "qs"), 

  # BULK WC MS-Empire
  tar_target(output_MSEmpPrep_BR123456_nofilt_WC, pD_MsEmpire_prep_specificBR(data_wc_bulk, BR=c(1,2,3,4,5,6),conds=c("NT","min_10","min_30","min_60"),TQVfilt=1.1,CQfilt=1.1),format = "qs"), 
  tar_target(MsEmpire_bulk_BR123456_nofilt_WC,  pD_MsEmpire_run_specific(output_MSEmpPrep_BR123456_nofilt_WC),format = "qs"), 
  tar_target(MsEmpire_bulk_BR123456_nofilt_WC_merged,  pD_merge_MSEmpire(MsEmpire_bulk_BR123456_nofilt_WC,quant="p.val"),format = "qs"), 
  
  # BULK WC and NUC
  tar_target(PSEA_bulk_nuc_10_WC, pD_PSEA_simple(MsEmpire_bulk_BR123456_nofilt_WC_merged, cond="min_10", min_prot=3, min_fraction = 0.05, ontology_sources =c("GO")),format = "qs"), 
  tar_target(PSEA_bulk_nuc_30_WC, pD_PSEA_simple(MsEmpire_bulk_BR123456_nofilt_WC_merged, cond="min_30", min_prot=3, min_fraction = 0.05, ontology_sources =c("GO")),format = "qs"), 
  tar_target(PSEA_bulk_nuc_60_WC, pD_PSEA_simple(MsEmpire_bulk_BR123456_nofilt_WC_merged, cond="min_60", min_prot=3, min_fraction = 0.05, ontology_sources =c("GO")),format = "qs"), 
  tar_target(bulk_PSEA_joined_WC, pD_PSEA_join_bulk(PSEA_bulk_nuc_10_WC, PSEA_bulk_nuc_30_WC, PSEA_bulk_nuc_60_WC, MsEmpire_bulk_BR123456_nofilt_WC_merged, Cond_col ="Cond", quant="log2FC"),format = "qs"),
  tar_target(PSEA_bulk_WC_plot, pD_plot_PSEA_bulk(bulk_PSEA_joined_WC, MsEmpire_bulk_BR123456_nofilt_WC_merged, keep=keep_PSEA_WC, 
                                                     conds=c("min_10","min_30","min_60"),title_plot="bulk_WC", spec_width=6.6,spec_height=5.5),format = "qs"), 
  
  # VOLCANO PLOT BULK WC and NUC
  tar_target(DA_plot_bulk_NUC_WC, pD_plot_DA_WC_NUC(MsEmpire_bulk_BR123456_merged, MsEmpire_bulk_BR123456_nofilt_WC_merged, SelectCond="min_60"),format = "qs"),
  
  #######################################################
  ############## Bulk purity analysis:  #################
  #######################################################
  tar_target(data_nuc_bulk_purity, pD_get_data(bulk_purity_fpath_timsTOFSCP, meta_bulk_path_timsTOFSCP, filter=F, TQV_threshold=1.1, CQV_threshold=0.3, proteotypic=T, ppm=5, mono_only=T, correct_IsoEnv=T),format = "qs"),
  tar_target(dat_bulk_maxLFQ_purity, diann_maxlfq(data_nuc_bulk_purity, group.header="Genes", id.header = "seqcharge", quantity.header = "Ms1.Area_iso"),format = "qs"),
  tar_target(dat_bulk_maxLFQ_m_purity, pD_melt_MaxLFQ(dat_bulk_maxLFQ_purity, meta_bulk_path_timsTOFSCP),format = "qs"),
  tar_target(dat_bulk_purity, pD_bulk_purity(dat_bulk_maxLFQ_m_purity, markers_fpath),format = "qs")

)

  
  
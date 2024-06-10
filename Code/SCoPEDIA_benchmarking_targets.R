
# Here we analyze the SCoPE-DIA benchamrking experiments for human nuclei mixed with yeast


list(
  
  #######################################################
  ##########  Human yeast quant evaluations:  ###########
  #######################################################

  tar_target(HY_dat, pD_get_data(HY_fpath, meta_bulk_path_timsTOFSCP, filter=F, proteotypic=F, ppm=5, mono_only=T, correct_IsoEnv=T),format = "qs"), 
  
  tar_target(HY_keepPrecs_6, pD_HY_keepPrecs(HY_dat, compression_factor=6),format="qs"), 
  tar_target(HY_filtered_6, pD_HY_filt(HY_dat, HY_keepPrecs_6$keep_precs),format="qs"), 
   
  tar_target(HY_keepPrecs_4, pD_HY_keepPrecs(HY_dat, compression_factor=4),format="qs"), 
  tar_target(HY_filtered_4, pD_HY_filt(HY_dat, HY_keepPrecs_4$keep_precs),format="qs"), 
   
  tar_target(HY_keepPrecs_2, pD_HY_keepPrecs(HY_dat, compression_factor=2),format="qs"), 
  tar_target(HY_filtered_2, pD_HY_filt(HY_dat, HY_keepPrecs_2$keep_precs),format="qs"), 
  
  tar_target(HY_keepPrecs_all, pD_HY_keepPrecs(HY_dat, compression_factor=Inf),format="qs"), 
  tar_target(HY_filtered_all, pD_HY_filt(HY_dat, HY_keepPrecs_all$keep_precs),format="qs"), 
  
  tar_target(HY_comparisons, pD_compare_HY_filt(HY_filtered_all, HY_filtered_2, HY_filtered_4, HY_filtered_6),format="qs") 

)

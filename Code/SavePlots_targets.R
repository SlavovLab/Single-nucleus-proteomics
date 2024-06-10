

list(
  
  ####################### save plots
  
  ### bulk
  tar_target(bulk_purity, save_plots(dat_bulk_purity$plot_Conds_density, myfilename="Figures/Bulk_purity.pdf",mywidth=4.7,myheight=4.7)),

  #### siRNA KDs
  tar_target(siRNA_WC_plot, save_plots_list(siRNA_WC_all$plot_list, myfilename="Figures/Supp/siRNA_KDs_WC",mywidth=3,myheight=3)),
  tar_target(siRNA_NUC_plot, save_plots_list(siRNA_NUC_all$plot_list, myfilename="Figures/Supp/siRNA_KDs_NUC",mywidth=3,myheight=3)),
  
  ### HY quant
  tar_target(HY_compression_hist, save_plots(HY_keepPrecs_6$compression_plot, myfilename="Figures/Supp/HY_compression_6x_hist.pdf",mywidth=4.7,myheight=4.7)),
  tar_target(HY_compression_coverage, save_plots(HY_comparisons$prot_count, myfilename="Figures/Supp/HY_compression_coverage.pdf",mywidth=4.7,myheight=4.7)),
  tar_target(HY_compression_accuracy, save_plots(HY_comparisons$prot_int, myfilename="Figures/Supp/HY_compression_accuracy.pdf",mywidth=4.7,myheight=4.7)),
  
  #### single nuc tehcnical
  tar_target(SN_keepnuclei_plot, save_plots(keep_cells$coverage_plot, myfilename="Figures/Supp/SN_keepNuclei.pdf",mywidth=6,myheight=4.7)),
  tar_target(SN_keepprecs_plot, save_plots(prec_correct$prec_plot, myfilename="Figures/Supp/SN_keepPrecs.pdf",mywidth=4.7,myheight=4.7)),
  
  tar_target(SN_AUC_hist, save_plots(AUC_purity_MLFQ$plot_AUC, myfilename="Figures/Supp/SN_AUC_hist.pdf",mywidth=4.7,myheight=4.7)),
  tar_target(SN_AUC, save_plots(AUC_purity_MLFQ$AUC_plot, myfilename="Figures/Supp/SN_AUC.pdf",mywidth=4.7,myheight=4.7)),
  tar_target(SN_AUC_filt, save_plots(AUC_purity_MLFQ$AUC_plot_filt, myfilename="Figures/Supp/SN_AUC_filt.pdf",mywidth=4.7,myheight=4.7))
  
)


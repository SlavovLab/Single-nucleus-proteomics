
list(
  
  # Bulk
  tar_target(writeBulkNuc, write.table(MsEmpire_bulk_BR123456_merged,
                                       file="Output_results/Bulk_Nuclei_DA.tsv", row.names = FALSE,sep="\t", quote=T),format = "qs"),
  
  tar_target(writeBulkWC, write.table(MsEmpire_bulk_BR123456_nofilt_WC_merged,
                                       file="Output_results/Bulk_WC_DA.tsv", row.names = FALSE,sep="\t", quote=T),format = "qs"),
  
  tar_target(Write_Bulk_PSEA_WC, write.table(bulk_PSEA_joined_WC %>% dplyr::arrange(desc(med_rel)),
                                                     file="Output_results/Bulk_WC_PSEA.tsv", row.names = FALSE,sep="\t", quote=T),format = "qs"),
  
  tar_target(Write_Bulk_PSEA_Nuc, write.table(bulk_PSEA_joined %>% dplyr::arrange(desc(med_rel)),
                                          file="Output_results/Bulk_Nuclei_PSEA.tsv", row.names = FALSE,sep="\t", quote=T),format = "qs"),
  
  #Single nuclei
  tar_target(Write_Cor_to_TS, write.table(Transport_Score_cor$cors,
                                      file="Output_results/Single_Nuclei_Cor_to_TS.tsv", row.names = FALSE,sep="\t", quote=T),format = "qs"),
  
  tar_target(Write_Cor_to_TS_timespecific_BR, write.table(Transport_Score_cor$cors_timeSpecific_BR,
                                                   file="Output_results/Single_Nuclei_Cor_to_TS_timespecific_BR.tsv", row.names = FALSE,sep="\t", quote=T),format = "qs"),
  
  tar_target(Write_SingleNuc, write.table(nucs_BC_MLFQ$imputed.BC_updated,
                                      file="Output_results/Single_Nuclei_Imputed.tsv", row.names = FALSE,sep="\t", quote=T),format = "qs"),
  tar_target(Write_SingleNuc_noimp, write.table(nucs_BC_MLFQ$mat.sc.imp_NA,
                                          file="Output_results/Single_Nuclei_NoImputation.tsv", row.names = FALSE,sep="\t", quote=T),format = "qs"),
  
  tar_target(Write_SingleNuc_CorTS_PSEA, write.table(PSEA_TS_cor %>% dplyr::arrange(desc(med_ab)),
              file="Output_results/Singel_Nuclei_PSEA_cor_to_TS.tsv", row.names = FALSE,sep="\t", quote=T),format = "qs"),
  
  tar_target(Write_SingleNuc_Zscore_cor_TS, write.table(Effect_on_transport_siRNA_weighted_100$cors_BR,
                                            file="Output_results/Single_Nuclei_Zscore_correlation_to_transport.tsv", row.names = FALSE,sep="\t", quote=T),format = "qs"),
  
  #Knock-downs
  tar_target(Write_KDs_Zs, write.table(Effect_on_transport_siRNA_weighted_100$Zs,
                                                file="Output_results/Knockdowns_Zscores.tsv", row.names = FALSE,sep="\t"),format = "qs"),
  
  tar_target(Write_KDs_PSEA_WC, write.table(PSEA_siRNA_KDs_WC_corOrder,
                                              file="Output_results/Knockdowns_WC_PSEA.tsv", row.names = FALSE,sep="\t", quote=T),format = "qs")
  
)

'%!in%' <- function(x,y)!('%in%'(x,y))

DIANN_read_folder <- function(folder_path){
  file_list <- list.files(folder_path)
  data_list <- list()
  for (file_name in file_list) {
    file_path <- file.path(folder_path, file_name)
    data <- data.frame(fread(file_path, 
                             select = c('Run','Lib.PG.Q.Value','Precursor.Id','Stripped.Sequence',
                                        'Precursor.Charge','Ms1.Area','Precursor.Translated','Protein.Group','Translated.Q.Value',
                                        'Proteotypic','Genes',"Protein.Names",'Modified.Sequence','Channel.Q.Value', 'First.Protein.Description')))
    data_list[[file_name]] <- data
  }
  df <- do.call("rbind", data_list)
  return(df)
}

pD_isotopicCO_modified <- function(df,ppm=5,monoisotopic_only=T){
  
  cnames <- colnames(df)
  df <- pD_seqcharge(df)
  df$seqcharge_run <- paste0(df$seqcharge, df$Run)
  
  mTRAQ0 <- list(C=7, H=12, N=2, O=1)
  mTRAQ4 <- list(C=4, H=12, N=1, O=1)
  mTRAQ8 <- list(C=1, H=12, N=0, O=1)
  
  # count number of each mTRAQ label type per precursor:
  df$d0 <- str_count(df$Modified.Sequence, "mTRAQ0|mTRAQ-K-0|mTRAQ-n-0")
  df$d4 <- str_count(df$Modified.Sequence, "mTRAQ4|mTRAQ-K-4|mTRAQ-n-4")
  df$d8 <- str_count(df$Modified.Sequence, "mTRAQ8|mTRAQ-K-8|mTRAQ-n-8")
  
  uni_prec <- df %>% dplyr::distinct(Modified.Sequence, .keep_all=T)
  
  deamid <- list(N=-1, H=-1, O=1)
  ox <- list(O=1)
  carba <- list(C=2, H=3, N=1, O=1)
  
  precs <- vector(mode = "list", length = nrow(uni_prec))
  # get chemical formula for each peptide
  for(i in 1:nrow(uni_prec)){
    tempseq <- paste0(uni_prec$Stripped.Sequence[i])
    modseq <- paste0(uni_prec$Modified.Sequence[i])
    numLab_d0 <- uni_prec$d0[i]
    numLab_d4 <- uni_prec$d4[i]
    numLab_d8 <- uni_prec$d8[i]
    
    ## other modifications
    d <- str_count(uni_prec$Precursor.Id[i], "UniMod:7")
    o <- str_count(uni_prec$Precursor.Id[i], "UniMod:35")
    c <- str_count(uni_prec$Precursor.Id[i], "UniMod:4")
    
    ch <- as.numeric(paste0(uni_prec$Precursor.Charge[i]))
    el_temp <- OrgMassSpecR::ConvertPeptide(tempseq)
    el_temp$C <- el_temp$C+(mTRAQ0$C*numLab_d0)+(mTRAQ4$C*numLab_d4)+(mTRAQ8$C*numLab_d8)+(carba$C*c)
    el_temp$H <- el_temp$H+(mTRAQ0$H*numLab_d0)+(mTRAQ4$H*numLab_d4)+(mTRAQ8$H*numLab_d8)+(carba$H*c)
    el_temp$N <- el_temp$N+(mTRAQ0$N*numLab_d0)+(mTRAQ4$N*numLab_d4)+(mTRAQ8$N*numLab_d8)+(carba$N*c)
    el_temp$O <- el_temp$O+(mTRAQ0$O*numLab_d0)+(mTRAQ4$O*numLab_d4)+(mTRAQ8$O*numLab_d8)+(carba$O*c)
    el_temp_df <- data.frame(el_temp)
    el_temp_df$charge <- paste0(ch)
    rownames(el_temp_df) <- modseq
    precs[[i]] <- el_temp_df
  }
  
  precs_df <- do.call("rbind", precs)
  precs_df <- precs_df[!grepl("K",row.names(precs_df)),] #remove K peptides because they have an extra label which will shift them far enough away to not be problematic (isotopically)
  
  precs_df$charge <- as.numeric(precs_df$charge)
  precs_df <- precs_df %>% dplyr::mutate(form = ifelse(S==0, paste0("C",precs_df$C,"H",precs_df$H,"N",precs_df$N,"O",precs_df$O),
                                                       paste0("C",precs_df$C,"H",precs_df$H,"N",precs_df$N,"O",precs_df$O,"S",precs_df$S)))
  data(isotopes)
  iso <- enviPat::isopattern(isotopes, precs_df$form, threshold=0.001, plotit=FALSE,
                    algo=1,rel_to=0,verbose=TRUE, return_iso_calc_amount=FALSE, charge=precs_df$charge)

  iso_env_list <- list()
  for(i in 1:length(iso)){
    expected_mz <- as.numeric(iso[[i]][1,1])
    iso_temp <- data.frame(iso[[i]])
    charge <- as.numeric(precs_df$charge[i])
    expected_mzs <- c(expected_mz, expected_mz+(1.0034)/charge, expected_mz+(2*1.0034)/charge, expected_mz+(3*1.0034)/charge, expected_mz+(4*1.0034)/charge, expected_mz+(5*1.0034)/charge)
    keep <-data.frame("X1" = sum(iso_temp$abundance[which(iso_temp$'m.z'>(expected_mzs[1]-expected_mzs[1]*ppm*0.000001)&iso_temp$'m.z'<(expected_mzs[1]+expected_mzs[1]*ppm*0.000001))]),
                      "X2" = sum(iso_temp$abundance[which(iso_temp$'m.z'>(expected_mzs[2]-expected_mzs[2]*ppm*0.000001)&iso_temp$'m.z'<(expected_mzs[2]+expected_mzs[2]*ppm*0.000001))]),
                      "X3" = sum(iso_temp$abundance[which(iso_temp$'m.z'>(expected_mzs[3]-expected_mzs[3]*ppm*0.000001)&iso_temp$'m.z'<(expected_mzs[3]+expected_mzs[3]*ppm*0.000001))]),
                      "X4" = sum(iso_temp$abundance[which(iso_temp$'m.z'>(expected_mzs[4]-expected_mzs[4]*ppm*0.000001)&iso_temp$'m.z'<(expected_mzs[4]+expected_mzs[4]*ppm*0.000001))]),
                      "X5" = sum(iso_temp$abundance[which(iso_temp$'m.z'>(expected_mzs[5]-expected_mzs[5]*ppm*0.000001)&iso_temp$'m.z'<(expected_mzs[5]+expected_mzs[5]*ppm*0.000001))]),
                      "X6" = sum(iso_temp$abundance[which(iso_temp$'m.z'>(expected_mzs[6]-expected_mzs[6]*ppm*0.000001)&iso_temp$'m.z'<(expected_mzs[6]+expected_mzs[6]*ppm*0.000001))]))
    row.names(keep) <- row.names(precs_df)[i]
    iso_env_list[[i]] <- keep
  }
  
  iso_env_df <- do.call("rbind", iso_env_list)
  iso_env_df$modseq <- row.names(iso_env_df)
  # join with the quantitative data, and correct for isotopic carryover
  dat <- df %>% left_join(iso_env_df, by =c("Modified.Sequence" = "modseq"))
  dat$num_labs <- dat$d0+dat$d4+dat$d8
  dat <- pD_channel(dat)
  if(monoisotopic_only==T){
    dat0 <- dat %>% dplyr::filter(channel_name=="mTRAQ0") %>% dplyr::mutate("y0" = (X5)/(X1+X5))
    dat4 <- dat %>% dplyr::filter(channel_name=="mTRAQ4") %>% dplyr::mutate("y4" = (X5)/(X1+X5))
    dat8 <- dat %>% dplyr::filter(channel_name=="mTRAQ8") %>% dplyr::mutate("y8" = (X5)/(X1+X5))
    dat_fin <- rbind.fill(dat0,dat4,dat8)
  } else{
    dat0 <- dat %>% dplyr::filter(channel_name=="mTRAQ0") %>% dplyr::mutate("y0" = (X5+X6)/(X1+X2+X5+X6))
    dat4 <- dat %>% dplyr::filter(channel_name=="mTRAQ4") %>% dplyr::mutate("y4" = (X5+X6)/(X1+X2+X5+X6))
    dat8 <- dat %>% dplyr::filter(channel_name=="mTRAQ8") %>% dplyr::mutate("y8" = (X5+X6)/(X1+X2+X5+X6))
    dat_fin <- rbind.fill(dat0,dat4,dat8)
  }
  
  if(monoisotopic_only==T){
    dat <- dat %>% dplyr::mutate("y" = (X5)/(X1+X5))
  } else{
    dat <- dat %>% dplyr::mutate("y" = (X5+X6)/(X1+X2+X5+X6))
  }
  dat_d <- dcast(dat_fin, seqcharge_run~channel_name,value.var="Ms1.Area")
  dat_d[is.na(dat_d)] <- 0
  
  dat_d_y <- dcast(dat, seqcharge_run~channel_name,value.var="y")
  dat_d_y <- dat_d_y %>% dplyr::rename("y0"="mTRAQ0","y4"="mTRAQ4","y8"="mTRAQ8")
  all_dat <- dat_d %>% full_join(dat_d_y, by =c("seqcharge_run"="seqcharge_run"))
  all_dat[is.na(all_dat)]<-0
  all_dat$d0_iso <- all_dat$mTRAQ0+all_dat$mTRAQ0*all_dat$y0
  all_dat$d4_iso <- all_dat$mTRAQ4 - all_dat$mTRAQ0*all_dat$y0 + all_dat$mTRAQ4*all_dat$y4 #remove what was from d0, and add what bled into d8
  all_dat$d8_iso <- all_dat$mTRAQ8 - all_dat$mTRAQ4*all_dat$y4 + all_dat$mTRAQ8*all_dat$y8 #remove what was from d4, and add what bled into d12
  #only do the isotopic correction for precursors with 1 label.. otherwise, the difference in mass is enough to avoid impactful isotopic carryover (such is the case with K-containing peptides).
  all_dat$d0_iso <- ifelse(is.na(all_dat$d0_iso), all_dat$mTRAQ0, all_dat$d0_iso)
  all_dat$d4_iso <- ifelse(is.na(all_dat$d4_iso), all_dat$mTRAQ4, all_dat$d4_iso)
  all_dat$d8_iso <- ifelse(is.na(all_dat$d8_iso), all_dat$mTRAQ8, all_dat$d8_iso)
  all_dat[all_dat<0]<-0 #replace negatives with 0
  all_dat_d0 <- all_dat %>% dplyr::select("seqcharge_run", "d0_iso") %>% dplyr::rename("Ms1.Area_iso" = "d0_iso") %>% na.omit()
  all_dat_d4 <- all_dat %>% dplyr::select("seqcharge_run", "d4_iso")%>% dplyr::rename("Ms1.Area_iso" = "d4_iso") %>% na.omit()
  all_dat_d8 <- all_dat %>% dplyr::select("seqcharge_run", "d8_iso")%>% dplyr::rename("Ms1.Area_iso" = "d8_iso") %>% na.omit()
  d0_dat <- dat[grepl("mTRAQ0|mTRAQ-K-0|mTRAQ-n-0", dat$Modified.Sequence),]
  d4_dat <- dat[grepl("mTRAQ4|mTRAQ-K-4|mTRAQ-n-4", dat$Modified.Sequence),]
  d8_dat <- dat[grepl("mTRAQ8|mTRAQ-K-8|mTRAQ-n-8", dat$Modified.Sequence),]
  
  d0_dat <- d0_dat %>% left_join(all_dat_d0, by =c("seqcharge_run" = "seqcharge_run"))
  d4_dat <- d4_dat %>% left_join(all_dat_d4, by =c("seqcharge_run" = "seqcharge_run"))
  d8_dat <- d8_dat %>% left_join(all_dat_d8, by =c("seqcharge_run" = "seqcharge_run"))
  dat <- rbind(d0_dat, d4_dat, d8_dat)
  dat$Ms1.Area_iso[is.na(dat$Ms1.Area_iso)] <-0
  dat <- dat %>% dplyr::select(all_of(cnames), "Ms1.Area_iso")
  rm(df)
  return(dat)
}

pD_seqcharge <- function(df){
  df$seqcharge <- paste0(df$Modified.Sequence, df$Precursor.Charge)
  #remove the channel-specific info
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ0\\)")
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ\\)")
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ-K-0\\)")
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ-n-0\\)")
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ4\\)")
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ-K-4\\)")
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ-n-4\\)")
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ8\\)")
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ-K-8\\)")
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ-n-8\\)")
  return(df)
}

pD_channel <- function(df){
  df <- df %>% dplyr::mutate("channel_name" = ifelse(grepl("-0|mTRAQ0", Modified.Sequence), "mTRAQ0",
                                                     ifelse(grepl("-4|mTRAQ4", Modified.Sequence), "mTRAQ4", "mTRAQ8")))
  return(df)
}

pD_filter <- function(dat,TQV_threshold=0.1,CQV_threshold=0.1){
  dat <- dat[dat$Translated.Q.Value<TQV_threshold,]
  dat <- dat[dat$Channel.Q.Value<CQV_threshold,]
  return(dat)
}

pD_get_data <- function(fpaths, meta_path, filter=T,TQV_threshold=0.5,CQV_threshold=0.5,proteotypic=T, ppm=5, mono_only=T,correct_IsoEnv=T){
  dat <- DIANN_read_folder(fpaths)
  if(proteotypic==T){
    dat <- dat[dat$Proteotypic==T,]
  }
  if(correct_IsoEnv==T){
    dat <- pD_isotopicCO_modified(dat, ppm=ppm, monoisotopic_only=mono_only) #isopopic envelope correction
  }
  m <- fread(meta_path)
  dat <- pD_seqcharge(dat)
  dat <- pD_channel(dat)
  dat <- dat[dat$Lib.PG.Q.Value<0.01,]
  if(filter==T){
    dat <- pD_filter(dat, TQV_threshold, CQV_threshold)
  }
  dat <- dat %>% dplyr::rename("prot" = "Protein.Group")
  dat <- dat %>% dplyr::select("Run","channel_name","Ms1.Area", any_of(c("Ms1.Area_iso")), "Proteotypic","Precursor.Translated","prot", "Protein.Names","Genes","First.Protein.Description","seqcharge","Translated.Q.Value","Channel.Q.Value")
  dat$run_chan <- paste0(dat$Run, dat$channel_name)
  m$run_chan <- paste0(m$Raw, "mTRAQ",m$Label)
  m$batch <- paste0(m$Label,"_",m$BioRep,"_",m$LC_batch)
  dat<-dat %>% inner_join(m, by =c("run_chan"="run_chan"))
  dat$File.Name <- dat$run_chan
  return(dat)
}

pD_MsEmpire_prep_specificBR <- function(Nuc_dat,BR=c(1,2,3,6),conds=c("NT","min_10","min_30","min_60"),TQVfilt=0.05,CQfilt=1.1){
  dat_fpaths_fin <- as.vector(NULL)
  mappings_fpaths_fin <- as.vector(NULL)
  BRs <- as.vector(NULL)
  Nuc_dat$id <- paste0(Nuc_dat$Genes,".",Nuc_dat$seqcharge)
  Nuc_dat <- Nuc_dat[Nuc_dat$Channel.Q.Value<CQfilt,]
  Nuc_dat <- Nuc_dat[Nuc_dat$Translated.Q.Value<TQVfilt,]
  for(i in 1:length(BR)){
    Nuc_temp <- Nuc_dat[Nuc_dat$BioRep==BR[i],]
    for(k in 2:length(conds)){
      conds_temp <- conds[c(1,k)]
      Nuc <- Nuc_temp[Nuc_temp$Celltype%in%conds_temp,]
      #### join with meta data:
      Nuc <- Nuc %>% dplyr::mutate("cond" = ifelse(grepl("NT", Celltype), "c1",
                                                   ifelse(grepl("10", Celltype), "c2", 
                                                          ifelse(grepl("30", Celltype),"c3", "c4"))))
      Nuc$sample <- paste0(Nuc$cond,".rep",Nuc$Rep,Nuc$BioRep)
      Nuc <- Nuc %>% dplyr::select("id","Ms1.Area_iso","sample")
      Nuc_d <- reshape2::dcast(Nuc,id~sample, value.var = "Ms1.Area_iso")
      Nuc_d[is.na(Nuc_d)] <- 0
      mappings <- Nuc  %>% dplyr::distinct(sample) %>% dplyr::mutate("condition" = ifelse(grepl("c1", sample), "0",
                                                                                          ifelse(grepl("c2", sample), "1", 
                                                                                                 ifelse(grepl("c3", sample), "2", "3"))))
      time <- Sys.time()
      #make unique filepaths, and store them so they can be loaded back in future
      dat_fpaths <- paste0("nuc_bulk_BR",BR[i],"_cond",conds[k],time,".tsv")
      mappings_fpath <- paste0("nuc_mapping_bulk_BR",BR[i],"_cond",conds[k],time,".tsv")
      write.table(Nuc_d,file=dat_fpaths, row.names = FALSE,sep="\t")
      write.table(mappings,file=mappings_fpath, row.names = FALSE,sep="\t")
      dat_fpaths_fin <- c(dat_fpaths,dat_fpaths_fin)
      mappings_fpaths_fin <- c(mappings_fpath,mappings_fpaths_fin)
      BRs <- c(paste0(BR[i]),BRs)
    }
  }
  fin <- list(dat_fpaths_fin=dat_fpaths_fin, mappings_fpaths_fin=mappings_fpaths_fin,BRs=BRs)
  return(fin)
}

pD_MsEmpire_run_specific <- function(dat_fpaths_msEmp){
  dat <- dat_fpaths_msEmp$dat_fpaths_fin
  mappings <- dat_fpaths_msEmp$mappings_fpaths_fin
  BRs <- dat_fpaths_msEmp$BRs
  empty_list <- list()
  for(i in 1:length(mappings)){
    data <- msEmpiRe::read.standard(dat[i], mappings[i],
                                    prot.id.generator=function(pep) unlist(strsplit(pep, "\\."))[1],
                                    signal_pattern="c.*rep.*")
    
    #extract the first two conditions
    conditions <- msEmpiRe::extract_conditions(data)
    conditions_temp <- conditions[, c(1,2)]
    
    #removing peptides that are detected in less than 2 samples per condition
    data_temp <- msEmpiRe::filter_detection_rate(data, condition=conditions_temp, rate=2)
    system.time(data_temp <- msEmpiRe::normalize(data_temp))
    set.seed(1234)
    
    #main analysis, system.time is just to measure the time of processing
    system.time(result <- msEmpiRe::de.ana(data_temp))
    result$BioRep <- paste0(BRs[i])
    result$Cond <- paste0(dat[i])
    empty_list[[i]] <- data.frame(result)
  }
  df <- do.call("rbind", empty_list)
  df$Gene <- df$prot.id

  return(df)
}

pD_merge_MSEmpire <- function(bulk_dat,quant="p.val"){
  bulk_dat$zscore <- p.to.Z(bulk_dat[[quant]]) #convert pval to zscore
  bulk_dat$zscore[is.na(bulk_dat$zscore)] <- max(bulk_dat$zscore,na.rm=T) 
  bulk_dat$zscore[is.infinite(bulk_dat$zscore)] <- max(bulk_dat$zscore,na.rm=T)
  #get zscore from each indiivudal biorep at bulk-level, take the mean across those bioreps and square it
  bulk_dat <- bulk_dat %>% dplyr::mutate("Cond"=ifelse(grepl("min_60",Cond),"min_60",
                                                       ifelse(grepl("min_30",Cond),"min_30","min_10")))
  bulk_dat <- bulk_dat %>% group_by(Gene,Cond) %>% add_count() %>%
    dplyr::mutate("z_joined" = (sum(sign(log2FC)*zscore,na.rm=T))/sqrt(n)) %>% 
    dplyr::mutate("log2FC"=mean(log2FC,na.rm=T)) %>% ungroup() %>%
    dplyr::distinct(Gene,Cond, .keep_all=T) %>% 
    dplyr::mutate(p.val = Z.to.p(z_joined)) %>% group_by(Cond) %>%
    dplyr::mutate(p.adj = p.adjust(p.val, method = "BH")) %>% ungroup() %>%
    dplyr::select(Gene,p.val,p.adj,Cond,log2FC)
  return(bulk_dat)
}

pD_BartelTest <- function(bulk_DA_all_rm_BR1_30min, p.adj_filter=1.1){
  bulk_DA <- bulk_DA_all_rm_BR1_30min %>% dplyr::mutate("time"=Cond) %>% 
    add_count(Gene) %>% filter(n==3) %>% group_by(Gene) %>% filter(median(p.adj)<p.adj_filter) %>%
    ungroup() %>% 
    dplyr::select("Gene","time","log2FC")
  pvals <- bulk_DA_all_rm_BR1_30min %>% group_by(Gene) %>%
    summarize(p.adj = min(p.adj)) %>% ungroup()
  bulk_DA_zero <- bulk_DA %>% distinct(Gene,.keep_all = T)
  bulk_DA_zero$time <- "min_0"
  bulk_DA_zero$log2FC <- 0
  both <- rbind(bulk_DA_zero,bulk_DA)
  both$time <- factor(both$time, levels=c("min_0","min_10","min_30","min_60"))
  vonNeumannBartel_result <- both  %>% arrange(time) %>% dplyr::group_by(Gene) %>% 
    dplyr::summarize(VNR = lawstat::bartels.test(log2FC)$parameter[1]) %>% ungroup() %>% 
    left_join(pvals, by =c("Gene"= "Gene")) %>% 
    dplyr::mutate("sig" = ifelse(p.adj<0.01, "pval<0.01", ifelse(p.adj>0.1, "pval>0.1","remove"))) %>%
    filter(!grepl("remove",sig))
  Bartel_plot <- ggplot(vonNeumannBartel_result, aes(x=VNR,fill=sig,y=sig, height = after_stat(density))) + 
    geom_density_ridges(stat = "binline", bins = 20, scale = 0.95, draw_baseline = FALSE)

  df <- vonNeumannBartel_result
  fin <- list(Bartel_plot=Bartel_plot, df=df)
  return(fin)
}

pD_timeseries_MSEmp <- function(dat, tofilt = NULL, DA_prots_to_label, p.adj_filter=1.1,
                                filename="TimeSeries", spec_width=6.5,spec_height=4.3, multiple=F){
  highFC <- dat[abs(dat$log2FC)>2,]
  highFC <- unique(highFC$Gene) 
  highFC <- highFC[highFC%!in%DA_prots_to_label]
  if(!is.null(tofilt)){
    dat <- dat %>% dplyr::filter(Gene%in%tofilt)  #filter for select proteins if a vector of them is provided
  }
  
  if(multiple==T){ #in case your input has cellular multiple fractions (e.g. whole-cells and nuclei)
    dat <- dat %>% dplyr::mutate("frac"= Fraction) %>% dplyr::mutate("GeneFrac" = paste0(Gene,"_",frac))
    keepProts <- dat %>% add_count(Gene,Fraction) %>% filter(n==3) %>% filter(Cond=="min_60") %>%
      distinct(Gene,Fraction,.keep_all=T) %>% filter(p.adj<p.adj_filter) %>% dplyr::select("GeneFrac","Cond","p.adj")
    bulk_DA <- dat %>% dplyr::mutate("time"=Cond) %>% filter(GeneFrac%in%keepProts$GeneFrac) %>%
      ungroup() %>% dplyr::select("Gene","time","log2FC","GeneFrac","Fraction") 
  } else{
    dat <- dat %>% dplyr::mutate("frac"="placeholder") %>% dplyr::mutate("GeneFrac" = paste0(Gene,"_",frac))
    keepProts <- dat %>% add_count(Gene) %>% filter(n==3) %>% filter(Cond=="min_60") %>% 
      distinct(Gene,.keep_all=T) %>% filter(p.adj<p.adj_filter) %>% dplyr::select("GeneFrac","Cond","p.adj") 
    bulk_DA <- dat %>% dplyr::mutate("time"=Cond) %>% filter(GeneFrac%in%keepProts$GeneFrac) %>%
      ungroup() %>% dplyr::select("Gene","time","log2FC","GeneFrac") 
  }
  bulk_DA_zero <- bulk_DA %>% distinct(GeneFrac,.keep_all = T)
  bulk_DA_zero$time <- "min_0"
  bulk_DA_zero$log2FC <- 0
  both <- rbind(bulk_DA_zero,bulk_DA)
  both$type <- "P-value < 0.05"
  # simulate a null distribution
  sim <- both
  set.seed(123)
  sim$log2FC <- rnorm(sim$log2FC, mean=0, sd=1)
  sim$type <- "Simulated null"
  all <- rbind(both,sim) 
  
  all$time <- factor(all$time, levels=c("min_0","min_10","min_30","min_60"))
  if(multiple==T){
    vonNeumannBartel_result <- all  %>% arrange(time) %>% dplyr::group_by(GeneFrac,type,Fraction) %>% dplyr::summarize(VNR = lawstat::bartels.test(log2FC)$parameter[1]) %>% ungroup()
  } else{
    vonNeumannBartel_result <- all  %>% arrange(time) %>% dplyr::group_by(GeneFrac,type) %>% dplyr::summarize(VNR = lawstat::bartels.test(log2FC)$parameter[1]) %>% ungroup()
  }
  df <- data.frame(vonNeumannBartel_result)
  df$type <- as.character(df$type)
  df$VNR <- as.numeric(df$VNR)
  simres <- df %>% filter(df$type=="Simulated null")
  realres <- df %>% filter(df$type=="P-value < 0.05")
  
  result <- t.test(simres$VNR,realres$VNR)
  
  Bartel_plot <- ggplot(df, aes(x=VNR, fill=type)) + 
    geom_histogram(alpha=0.9,color="black",binwidth = 0.25) + 
    facet_grid(~type) + theme_bw() + theme(legend.position = "none") +
    labs(x="Bartel's rank VNR", fill = "Dif abundance", y="Proteins", subtitle = paste0("pval: ",result$p.value)) +
    scale_fill_manual(values =c("#3d405b","gray45"))
  if("Fraction" %in% names(bulk_DA)) {
    Bartel_plot <- ggplot(df, aes(x=VNR, fill=type)) + 
      geom_histogram(alpha=0.9,color="black",binwidth = 0.25) + 
      facet_grid(~Fraction+type) + theme_bw() + theme(legend.position = "none") +
      labs(x="Bartel's rank VNR", fill = "Dif abundance", y="Proteins", subtitle = paste0("pval: ",result$p.value)) +
      scale_fill_manual(values =c("#3d405b","gray45"))
  }
  ggsave(paste0("Figures/Bartel_",filename,".pdf"),Bartel_plot,width=3.8,height=3.3)
  
  dats <- both %>% dplyr::mutate("time_actual" = ifelse(grepl("10",time),10,
                                                        ifelse(grepl("30",time),30,
                                                               ifelse(grepl("60",time),60,0))))
  timeseries_plot <- ggplot(dats,aes(x=time_actual,y=log2FC,group=GeneFrac))+geom_line(alpha=0.1) + theme_classic() +
    geom_line(data=dats %>% filter(Gene%in%highFC), alpha=0.9, color="blue") + xlim(0,63)+
    geom_label_repel(data=dats %>% filter(Gene%in%highFC) %>% filter(time=="min_60"), aes(label=Gene),max.overlaps = Inf,size=3)+
    geom_line(data=dats %>% filter(Gene%in%DA_prots_to_label), alpha=0.9, color="red") + xlim(0,63)+
    geom_label_repel(data=dats %>% filter(Gene%in%DA_prots_to_label) %>% filter(time=="min_60"), aes(label=Gene),max.overlaps = Inf,size=3)
  
  if("Fraction" %in% names(dats)) {
    timeseries_plot <- ggplot(dats,aes(x=time_actual,y=log2FC,group=GeneFrac))+geom_line(alpha=0.4) + theme_classic() +
      geom_line(data=dats %>% filter(Gene%in%highFC), alpha=0.9, color="blue") + xlim(0,63)+
      geom_label_repel(data=dats %>% filter(Gene%in%highFC) %>% filter(time=="min_60"), aes(label=Gene),max.overlaps = Inf,size=3) +
      geom_line(data=dats %>% filter(Gene%in%DA_prots_to_label), alpha=1, color="red") +
      geom_label_repel(data=dats %>% filter(Gene%in%DA_prots_to_label) %>% filter(time=="min_60"), aes(label=Gene),max.overlaps = Inf,size=3) +
      facet_wrap(~Fraction) 

  }
  ggsave(paste0("Figures/TimeSeries_",filename,".pdf"),timeseries_plot,width=spec_width,height=spec_height)
  
  
  df_nofilt <- both %>% arrange(time) %>% dplyr::group_by(Gene) %>% dplyr::summarize(VNR = lawstat::bartels.test(log2FC)$parameter[1]) %>% ungroup()
  df_nofilt <- data.frame(df_nofilt)
  fin <- list(Bartel_plot=Bartel_plot, timeseries_plot=timeseries_plot, df=df_nofilt)
  return(fin)
}

pD_combine_data <- function(dat = list(MsEmpire_bulk_BR123456_nofilt_WC_merged,MsEmpire_bulk_BR123456_merged),labels=c("WC","NUC")){
  fin_list <- list()
  for(i in 1:length(dat)){
    temp <- dat[[i]]
    temp$Fraction <- paste0(labels[i])
    fin_list[[i]] <- temp
  }
  fin <- do.call("rbind", fin_list)
  return(fin)
}

pD_pval_filt <- function(dat,genefilter=nucleoporins,threshold=0.05){
  dat <- dat %>% filter(Gene%in%genefilter) %>% filter(p.adj<threshold)
  uniprots <- unique(dat$Gene)
  fin <- list(uniprots=uniprots,dat=dat)
  return(fin)
}

pD_passive_diffusion <- function(bulk_DA, masses_fpath, cond=c("min_60")){
  bulk_DA_lim <- bulk_DA
  bulk_DA_lim <- bulk_DA_lim[bulk_DA_lim$p.adj<0.05,] %>% distinct(Gene,Cond,.keep_all=T)
  uniprot_masses <- read.delim(masses_fpath)
  dat4 <- bulk_DA_lim %>% inner_join(uniprot_masses, by =c("Gene" = "Gene.Names..primary."))
  dat4 <- na.omit(dat4)
  dat4_cors <-dat4 %>% group_by(Cond) %>%dplyr::mutate(sp_cor = paste0("sp r = ",round(cor(abs(log2FC),log10(Mass), method = "spearman"),2))) %>%
    dplyr::mutate(pr_cor = paste0("pr r = ",round(cor(abs(log2FC),log10(Mass), method = "pearson"),2))) %>% 
    dplyr::mutate(sp_cor_pval = paste0("sp pval = ",cor.test(abs(log2FC),log10(Mass),method="spearman")$p.value)) %>% ungroup() %>%
    distinct(Cond,.keep_all=T)
  
  a <- ggplot(dat4, aes(x=log10(Mass), y =abs(log2FC))) + 
    geom_point(shape=21, color="black",fill="grey70",alpha=0.7) + 
    facet_wrap(~Cond)+
    geom_smooth(method = "lm") +    
    stat_cor(method="spearman")+
    #geom_vline(xintercept = log10(40000), color="red", linetype="dashed",linewidth=0.75) +
    labs(title = "Protein mass and transport",
         x = "Log10, Protein mass (Da)", y= "|Log2, fold change|") + theme_bw()
  ggsave("Figures/Bulk_passive_dif.pdf",a,width=6.5,height=4)
  return(a)
}

pD_mass_FC_time <- function(bulk_DA, masses_fpath, cond=c("min_60")){
  bulk_DA_lim <- bulk_DA
  keep_prots <- bulk_DA_lim %>% add_count(Gene) %>% filter(n==3 & p.adj<0.05) %>% distinct(Gene,.keep_all=T)
  bulk_DA_lim <- bulk_DA_lim[bulk_DA_lim$Gene%in%keep_prots$Gene,]
  uniprot_masses <- read.delim(masses_fpath)
  dat4 <- bulk_DA_lim %>% inner_join(uniprot_masses, by =c("Gene" = "Gene.Names..primary."),relationship = "many-to-many")
  dat_zeros <- dat4 %>% distinct(Gene,.keep_all=T) %>% dplyr::mutate("Cond"="NT") %>% dplyr::mutate("log2FC"=0)
  all <- rbind(dat_zeros,dat4)
  all <- all %>% dplyr::mutate("time"=ifelse(grepl("NT",Cond),0,
                                             ifelse(grepl("10",Cond),10,
                                                    ifelse(grepl("30",Cond),30,60))))
  all <- na.omit(all)
  all$log10_Mass <- log10(all$Mass)
  set.seed(123)
  predictions <- all %>%
    group_by(Gene) %>%
    do({
      #interpolate data acros timepoints with 3rd degree poly
      model <- lm(log2FC ~ poly(time, 3, raw = TRUE), data = .)
      time_new <- seq(min(.$time), max(.$time), length.out = 300)
      predicted_log2FC <- predict(model, newdata = data.frame(time = time_new))
      data.frame(Gene = unique(.$Gene), time = time_new, log2FC_predicted = predicted_log2FC)
    }) %>%
    ungroup() 
  
  predictions_mass <- predictions %>% inner_join(uniprot_masses, by =c("Gene" = "Gene.Names..primary."),relationship = "many-to-many")
  predictions_mass$log10_Mass <- log10(predictions_mass$Mass)
  
  time_bins <- 300
  mass_bins <- 10
  mass_breaks <- seq(min(predictions_mass$log10_Mass), max(predictions_mass$log10_Mass), length.out = mass_bins + 1)
  mass_bin_labels <- sprintf('[%.2f, %.2f)', head(mass_breaks, -1), tail(mass_breaks, -1))
  time_breaks <- seq(min(predictions_mass$time, na.rm = TRUE), max(predictions_mass$time, na.rm = TRUE), length.out = time_bins + 1)
  predictions_mass$bin_time <- base::cut(predictions_mass$time, breaks = time_breaks, include.lowest = TRUE, labels = FALSE)

  ###
  setDT(predictions_mass)
  
  #sort proteins by mass
  predictions_mass <- predictions_mass[order(Mass)]
  n <- nrow(predictions_mass)
  progress_interval <- n
  predictions_mass[, Avg_log2FC := sapply(.I, function(idx) {
    if (idx %% progress_interval == 0 || idx == n) {
      cat(sprintf("Processing %d of %d rows...\n", idx, n))
      flush.console() 
    }
    calculate_moving_average_dt(predictions_mass, idx)
  }), by = bin_time]
  
  predictions_mass$Mass_fac <- as.factor(log10(predictions_mass$Mass))
  weighted_avg_time <- predictions_mass %>%
    group_by(Mass_fac,Mass,Gene,Protein.names) %>%
    summarise(weighted_avg_time = sum(bin_time * (Avg_log2FC)) / sum(Avg_log2FC)) %>%
    ungroup()
  
  averaged_data <- predictions_mass %>% group_by(Mass) %>% dplyr::mutate("adj_val" = Avg_log2FC/max(Avg_log2FC)) %>% ungroup()
  averaged_data <- averaged_data %>% arrange(Mass)
  
  dyn_mass <- ggplot(averaged_data, aes(y = bin_time/5, x = Mass_fac)) +
    geom_tile(aes(fill = adj_val)) +
    scale_fill_viridis_c(option="mako") +
    labs(y = "Time (min)", x = "Protein (smallest to largest)", fill = "relative |Log2, Treated / NT|") +
    theme_classic() + theme(axis.text.x = element_blank(),
                            axis.ticks.x = element_blank(),
                            legend.position="top") +
    geom_point(data=weighted_avg_time, aes(y=weighted_avg_time/5,x=Mass_fac),color="red",alpha=0.4,size=1) 
  ggsave("Figures/bulk_mass_dynamics.pdf",dyn_mass,width=4.8,height=6.7)
  sp_cor <- cor((log10(weighted_avg_time$Mass)),(weighted_avg_time$weighted_avg_time),method="spearman")
  pr_cor <- cor((log10(weighted_avg_time$Mass)),(weighted_avg_time$weighted_avg_time),method="pearson")
  
  center_time <- ggplot(weighted_avg_time, aes(x=log10(Mass), y=weighted_avg_time/5)) + 
    geom_point(shape=21,fill="grey70",color="black",alpha=0.7,size=2) + 
    theme_classic() + geom_vline(xintercept=log10(40000),color="black",linetype="dashed",linewidth=0.75) +
    labs(y="Log2(FC) weighted average time",x="Log10, protein mass (Da)",
         subtitle=paste0("prots=",nrow(weighted_avg_time))) +
    stat_cor(method = "spearman") +
    geom_smooth(data = weighted_avg_time %>% filter(Mass<40000), method = "lm",linewidth=0.75) +
    geom_smooth(data = weighted_avg_time %>% filter(Mass>40000), method = "lm",linewidth=0.75)
  ggsave("Figures/bulk_mass_dynamics_time.pdf",center_time,width=3.1,height=2.8)
  
  return(list(all=all,dyn_mass=dyn_mass,center_time=center_time))
}

calculate_moving_average_dt <- function(df, index, n_neighbors = 40) {
  current_mass <- df$Mass[index]
  
  #find neighbors within bin_time
  df_time_bin <- df[bin_time == df$bin_time[index]]
  greater_indices <- which(df_time_bin$Mass > current_mass)
  lesser_indices <- which(df_time_bin$Mass < current_mass)
  
  if (length(greater_indices) > n_neighbors) {
    greater_indices <- head(greater_indices, n_neighbors)
  }
  if (length(lesser_indices) > n_neighbors) {
    lesser_indices <- tail(lesser_indices, n_neighbors)
  }
  
  #calculating the median across all neighbors
  neighbors_log2FC <- abs(df_time_bin$log2FC_predicted[c(lesser_indices, greater_indices)])
  return(median(neighbors_log2FC, na.rm = TRUE))
}

pD_PSEA_simple <- function(dat, cond="min_10", min_prot=3, min_fraction=0.05,ontology_sources =c("GO")){
  dat <- dat %>% dplyr::select("log2FC","Gene","Cond") %>% 
    dplyr::arrange(desc(log2FC))
  dat <- dat[dat$Cond==paste0(cond),]
  res <- gprofiler2::gost(dat$Gene, organism="hsapiens", ordered_query=T, user_threshold=0.01, sources=ontology_sources,correction_method='fdr',evcodes=T)
  res2 <- res$result
  res2 <- res2[res2$intersection_size >= min_prot,] #require certain number of proteins per GO term
  res2 <- res2[(1-res2$precision) >= min_fraction,] #require certain fraction of completeness
  res_meta <- res$meta
  res2$med_ab <- 0 #place holder
  res2$mean_ab <- 0 #place holder
  colnum <- ncol(res2)
  for(k in 1:nrow(res2)){
    res2_temp <- res2[k,]
    res2_temp <- unlist(strsplit(res2_temp$intersection,split=","))
    dat_temp_abs <- dat[dat$Gene%in%res2_temp,]
    dat_temp_abs_med <- median(dat_temp_abs$log2FC,na.rm=T) #get median value
    dat_temp_abs_mean <- mean(dat_temp_abs$log2FC,na.rm=T) #get mean value
    
    res2[k,colnum-1] <- dat_temp_abs_med
    res2[k,colnum] <- dat_temp_abs_mean
  }
  return(res2)
}

pD_PSEA_simple_v2 <- function(dat, min_prot=3, min_fraction=0.05,ontology_sources =c("GO")){
  dat <- dat %>% dplyr::select("pr_r","Genes") %>% 
    dplyr::arrange(desc(pr_r))
  res <- gprofiler2::gost(dat$Genes, organism="hsapiens", ordered_query=T, user_threshold=0.01, sources=ontology_sources,correction_method='fdr',evcodes=T)
  res2 <- res$result
  res2 <- res2[res2$intersection_size >= min_prot,] #require certain number of proteins per GO term
  res2 <- res2[(1-res2$precision) >= min_fraction,] #require certain fraction of completeness
  res_meta <- res$meta
  res2$med_ab <- 0 #place holder
  res2$mean_ab <- 0 #place holder
  colnum <- ncol(res2)
  for(k in 1:nrow(res2)){
    res2_temp <- res2[k,]
    res2_temp <- unlist(strsplit(res2_temp$intersection,split=","))
    dat_temp_abs <- dat[dat$Genes%in%res2_temp,]
    dat_temp_abs_med <- median(dat_temp_abs$pr_r,na.rm=T) #get median value
    dat_temp_abs_mean <- mean(dat_temp_abs$pr_r,na.rm=T) #get mean value
    
    res2[k,colnum-1] <- dat_temp_abs_med
    res2[k,colnum] <- dat_temp_abs_mean
  }
  
  # columns to character vector if list
  res2 <- convert_list_columns(res2)
  return(res2)
}

convert_list_columns <- function(df) {
  for (col in colnames(df)) {
    if (is.list(df[[col]])) {
      df[[col]] <- sapply(df[[col]], toString)
    }
  }
  return(df)
}

pD_PSEA_join_bulk <- function(PSEA_bulk_nuc_10, PSEA_bulk_nuc_30, PSEA_bulk_nuc_60, MsEmpire_bulk_BR123456_merged, Cond_col ="Cond", quant="log2FC"){
  PSEA_bulk_nuc_10$Cond <- "min_10"
  PSEA_bulk_nuc_30$Cond <- "min_30"
  PSEA_bulk_nuc_60$Cond <- "min_60"
  GOresults <- rbind(PSEA_bulk_nuc_10,PSEA_bulk_nuc_30,PSEA_bulk_nuc_60)
  unicond <- unique(MsEmpire_bulk_BR123456_merged[[Cond_col]])
  uni_term_id <- unique(GOresults$term_id)
  dat_list <- list()
  for(i in 1:length(uni_term_id)){
    res2_temp <- GOresults[GOresults$term_id==paste0(uni_term_id[i]),] 
    res2_temp_genes <- data.frame("Gene" = unlist(strsplit(res2_temp$intersection,split=","))) %>% distinct(Gene)
    keep_genes <- MsEmpire_bulk_BR123456_merged[MsEmpire_bulk_BR123456_merged$Gene%in%res2_temp_genes$Gene,] %>% 
      distinct(Gene, !!sym(Cond_col)) %>% add_count(Gene) %>% 
      filter(n==length(unicond)) %>% distinct(Gene) #require each gene to be quantified in all conditions
    for(k in 1:length(unicond)){
      temp_cond <- MsEmpire_bulk_BR123456_merged[MsEmpire_bulk_BR123456_merged[[Cond_col]]==paste0(unicond[k]),] %>% 
        filter(Gene%in%keep_genes$Gene)
      dat_temp_abs_med <- median(temp_cond[[quant]],na.rm=T) #get median value
      dat_temp_abs_mean <- mean(temp_cond[[quant]],na.rm=T) #get mean value
      res2_temp2 <- res2_temp[res2_temp$Cond==paste0(unicond[k]),] %>% distinct(term_id,.keep_all=T) %>%
        dplyr::select("source","term_id","term_name","p_value")
      if(nrow(res2_temp2)==0){
        res2_temp2 <- res2_temp %>% distinct(term_id,.keep_all=T) %>%
          dplyr::select("source","term_id","term_name","p_value") %>% dplyr::mutate("p_value"=0.05) #placeholder in case the GO term is not present
      }
      res2_temp2$med_ab <- as.numeric(paste0(dat_temp_abs_med))
      res2_temp2$mean_ab <- as.numeric(paste0(dat_temp_abs_mean))
      res2_temp2$Cond <- paste0(unicond[k])
      keep_genes_temp <- paste(keep_genes, collapse = ";")
      res2_temp2$intersected_prots <- as.vector(keep_genes_temp)
      res2_temp2$num_prots <- nrow(keep_genes)
      res2_temp2$term_name_numprot <- paste0(res2_temp2$term_name," [n=",nrow(keep_genes),"]")
      dat_list[[i*(length(unicond))-length(unicond)+k]] <- res2_temp2
    }
  }
  fin <- do.call("rbind",dat_list)
  #all other values are relative to the NT, for that reason I can create the time zero by setting it all values to zero,
  # then computing relative abundances together.
  fin_zero <- fin %>% distinct(term_id, .keep_all=T)
  fin_zero$Cond <- "NT"
  fin_zero$med_ab <- 0 
  fin_zero$mean_ab <- 0
  fin2 <- rbind(fin,fin_zero)
  fin2 <- fin2 %>% group_by(term_id) %>% dplyr::mutate("med_rel" = med_ab-mean(med_ab,na.rm=T)) %>%
    dplyr::mutate("mean_rel" = mean_ab-mean(mean_ab,na.rm=T)) %>% ungroup()
  return(fin2)
}

pD_plot_PSEA_bulk <- function(PSEA_bulk, bulk, keep=c("GO:00000001"), conds=c("min_10","min_30","min_60"),title_plot="PSEA", spec_width=6.6,spec_height=8.5){
  all <- PSEA_bulk[PSEA_bulk$term_id%in%keep,]
  all$Cond <- factor(all$Cond, levels=c("NT","min_10","min_30","min_60"))
  med_data <- all %>%
    filter(Cond %in% c("min_30", "min_60")) %>%
    group_by(term_name_numprot) %>%
    summarize(mean_value = mean(med_rel))
  all$term_name_numprot <- factor(all$term_name_numprot, levels = med_data$term_name_numprot[order(med_data$mean_value, decreasing = TRUE)])
  all <- all %>% dplyr::mutate("med_rel_adj"=ifelse(med_rel< -0.5,-0.5,
                                                    ifelse(med_rel>0.5,0.5,med_rel)))
  PSEA_plot_adj_med <- ggplot(all, aes(x=Cond,y=term_name_numprot,fill=med_rel_adj)) + geom_tile(color="black") + theme_classic() + 
    scale_fill_gradientn(colors = c("#394d71", "#5f97ce", "white", "#f89b62", "#f05056"),
                         values = scales::rescale(c(min(all$med_rel_adj), 
                                                    mean(range(all$med_rel_adj)), 
                                                    max(all$med_rel_adj))),
                         na.value = "grey50",
                         guide = "colourbar") + labs(fill="Log2, Median ab.") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position="top")
  ggsave(paste0("Figures/PSEA_",title_plot,".pdf"),PSEA_plot_adj_med,width=spec_width,height=spec_height)
  
  return(PSEA_plot_adj_med)
}

pD_plot_prot_bulkHeatmap <- function(dat, bulk, keep=bulk_prots_toPlot, conds=c("min_10","min_30","min_60"),title_plot="PSEA", spec_width=6.6,spec_height=49.5){
  keepprots <- read.delim(keep) %>% filter(Gene!="PAXX") #FC far too large to be compared with everything else
  dat <- dat %>% inner_join(keepprots, by=c("Gene"="Gene"),relationship = "many-to-many")
  dat_zero <- dat %>% distinct(Gene, .keep_all=T)
  dat_zero$log2FC <- 0
  dat_zero$Cond <- "NT"
  dat <- rbind(dat,dat_zero)
  dat$Cond <- factor(dat$Cond, levels=c("NT","min_10","min_30","min_60"))
  dat <- dat %>% dplyr::group_by(Gene) %>% dplyr::mutate("rel_ab"=log2FC-mean(log2FC)) %>% ungroup()
  dat <- dat %>% left_join(uniprot_masses, by =c("Gene"="Gene.Names..primary."))
  dat$Protein.names <- sub(" \\(.*", "", dat$Protein.names)
  dat$prot <- paste0(dat$Protein.names," [",dat$Gene,"]")
  med_data <- dat %>%
    filter(Cond %in% c("min_30", "min_60")) %>%
    group_by(prot) %>%
    summarize(mean_value = mean(rel_ab))
  dat$prot <- factor(dat$prot, levels = med_data$prot[order(med_data$mean_value, decreasing = TRUE)])
  dat <- dat %>% dplyr::mutate("rel_ab_adj"=ifelse(rel_ab< -0.5,-0.5,
                                                   ifelse(rel_ab>0.5,0.5,rel_ab)))
  
  prot_plot <- ggplot(dat, aes(x=Cond,y=prot,fill=rel_ab_adj)) + geom_tile(color="black") + theme_classic() + 
    facet_wrap(~group_type,ncol=1)+
    scale_fill_gradientn(colors = c("#394d71", "#5f97ce", "white", "#f89b62", "#f05056"),
                         limits = c(-0.5, 0.5), 
                         na.value = "grey50",
                         guide = "colourbar") + labs(fill="Log2, Median ab.") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position="top")
  
  ggsave(paste0("Figures/Protein_bulk",title_plot,".pdf"),prot_plot,width=spec_width,height=spec_height)
  
  return(PSEA_plot_adj_med)
}

pD_melt_MaxLFQ <- function(dat,meta_bulk_path){
  pD_1 <- data.frame(dat)
  pD_3_m <- melt(as.matrix(pD_1))
  m <- fread(meta_bulk_path)
  m$run_chan <- paste0(m$Raw,"mTRAQ",m$Label)
  m$run_chan <- gsub('\\-', '.', m$run_chan)
  pD_3_m <- pD_3_m %>% left_join(m, by=c("Var2"="run_chan"))
  pD_3_m <- pD_3_m %>% dplyr::rename("Genes"="Var1", "run_chan"="Var2")
  pD_3_m$Gene <- pD_3_m$Genes
  pD_3_m$norm_prot <- log2(as.numeric(pD_3_m$value))
  return(pD_3_m)
}

pD_bulk_purity <- function(dat, markers_fpath){
  ##### NOTE: the histones were selected to not map to any other location other than chromosome (e.g. not even nucleoplasm)
  markers <- read.delim(markers_fpath)
  dat$Celltype_rep <- paste0(dat$Celltype,dat$BioRep)
  dat <- na.omit(dat)
  uni_cond <- unique(dat$Celltype_rep)
  dat <- dat %>% left_join(markers, by =c("Genes" = "Gene_name")) 
  dat <- dat %>% 
    dplyr::group_by(Celltype_rep, Genes, Type) %>%
    dplyr::mutate("value" = mean(value,na.rm=T)) %>% 
    ungroup() %>% distinct(Celltype_rep, Genes, Type, .keep_all=T)
  dat$value[is.na(dat$value)] <- 0
  dat <- dat[!is.na(dat$value),]
  
  uni_celltype <- unique(dat$Type)
  list_temp <- list()
  for(i in 1:length(uni_cond)){
    temp_df <- dat[dat$Celltype_rep==paste0(uni_cond[i]),] %>% dplyr::select("Type", "Genes", "value", "Subcellular")  
    temp_df <- reshape2::dcast(temp_df, Subcellular+Genes~Type, value.var = "value")
    
    temp_df$WC_log2 <- log2(temp_df$WC)
    temp_df$NUC_log2 <- log2(temp_df$NUC)
    temp_df_His <- temp_df[grepl("His|Nuc", temp_df$Subcellular),]
    temp_df_His <- temp_df_His[!is.infinite(temp_df_His$WC_log2),]
    temp_df_His <- temp_df_His[!is.infinite(temp_df_His$NUC_log2),]
    X <- temp_df_His %>% dplyr::select("WC_log2", "NUC_log2")
    X_rat <- median(temp_df_His$WC_log2-temp_df_His$NUC_log2)
    temp_df$WC_nuc <- temp_df$WC_log2 - temp_df$NUC_log2  #######
    temp_df$WC_nuc <- temp_df$WC_nuc-X_rat #normalize based on histones.  #########
    temp_df$Celltype_rep <- paste0(uni_cond[i])
    temp_df$WC_nuc[is.na(temp_df$WC_log2)] <- -8
    temp_df$WC_nuc[is.na(temp_df$NUC_log2)] <- 8
    temp_df$WC_nuc[is.na(temp_df$NUC_log2)&is.na(temp_df$WC_log2)] <- NA
    temp_df <- temp_df[!is.na(temp_df$WC_nuc),]
    
    list_temp[[i]] <- temp_df
  }
  fin <- do.call("rbind", list_temp)
  fin1 <- fin[fin$WC_nuc!=8.000&fin$WC_nuc!=-8.000,] #remove values that were filled-in
  fin_med_ratio <- fin1 %>% dplyr::group_by(Genes) %>% 
    dplyr::mutate(med_NUC_log2 = median(NUC_log2,na.rm=T)) %>% 
    dplyr::mutate(med_WC_log2 = median(WC_log2,na.rm=T)) %>% ungroup() %>%
    dplyr::mutate(med_WC_NUC_log2_rat = med_WC_log2-med_NUC_log2) %>%
    dplyr::distinct(Genes,.keep_all=T) %>%
    dplyr::select("Genes","med_WC_NUC_log2_rat","med_NUC_log2","med_WC_log2")
  fin_all <- fin[!is.na(fin$Subcellular),] #remove other prots
  plot_sepConds <- ggplot(fin_all, aes(x=(WC_nuc))) +
    facet_grid(
      rows = vars(reorder(Subcellular,-WC_nuc,na.rm=T)), scales = "free", space = "free") + 
    geom_density_ridges2(
      aes(point_color = Celltype_rep, point_fill = Celltype_rep, y=Celltype_rep, fill=Celltype_rep),
      alpha = 0.3, point_alpha = 0.5, jittered_points = TRUE, quantile_lines = TRUE, quantiles = 2) + 
    labs(x="Log2, Whole cell / Nuclei", y = "Subcellular localization", 
         title="Relative abundance of subcellular marker proteins",point_color="", point_fill="") +
    theme_bw() + theme(legend.position = "none")
  fin_med <- fin %>% dplyr::group_by(Genes) %>% dplyr::mutate(med_rat = median(WC_nuc,na.rm=T)) %>%
    add_count() %>% dplyr::filter(n/length(uni_cond)>0.75) %>% 
    ungroup() %>% #must be present in >75% conditions
    distinct(Genes, .keep_all=T)  %>% na.omit()
  plot_CombinedConds <- ggplot(fin_med, aes(x=(med_rat))) +
    facet_grid(
      rows = vars(reorder(Subcellular,-med_rat,na.rm=T)), scales = "free", space = "free") + 
    geom_density_ridges2(
      aes(point_color = Subcellular, point_fill = Subcellular, y=Subcellular, fill=Subcellular),
      alpha = 0.3, point_alpha = 0.5, jittered_points = TRUE, quantile_lines = TRUE, quantiles = 2) + 
    labs(x="Log2, Whole cell / Nuclei", y = "Subcellular localization", 
         title="Relative abundance of subcellular marker proteins",point_color="", point_fill="") +
    theme_bw() + theme(legend.position = "none")
  
  plot_CombinedConds_box <- ggplot(fin_med, aes(x=(med_rat), y=Subcellular, fill=Subcellular)) +
    facet_grid(
      rows = vars(reorder(Subcellular,-med_rat,na.rm=T)), scales = "free", space = "free") + 
    geom_boxplot(alpha=0.8,outlier.shape = NULL) + geom_point(alpha=0.2)+
    labs(x="Log2, Whole cell / Nuclei", y = "Subcellular localization", 
         title="Relative abundance of subcellular marker proteins",point_color="", point_fill="") +
    theme_bw() + theme(legend.position = "none",
                       strip.background = element_blank(),
                       strip.text.y = element_blank())
  
  mean_values <- fin_med %>%
    group_by(Subcellular) %>%
    summarize(mean_value = mean(med_rat))
  fin_med$Subcellular <- factor(fin_med$Subcellular, levels = mean_values$Subcellular[order(mean_values$mean_value)])
  fin_med <- fin_med %>% dplyr::mutate("Subcell3" = ifelse(grepl("Nuc|His", Subcellular), "Nucleus","Other"))
  plot_Conds_density <- ggplot(fin_med, aes(x=(med_rat))) +
    facet_grid(
      rows = vars(reorder(Subcell3,-med_rat,na.rm=T)), scales = "free", space = "free") +
    geom_density_ridges2(
      aes(point_color = Subcellular, point_fill = Subcellular, y=Subcellular, fill=Subcellular),
      alpha = 0.3, point_alpha = 0.5, jittered_points = TRUE, quantile_lines = TRUE, quantiles = 2) + 
    labs(x="Log2, Whole cell / Nuclei", y = "Subcellular localization", 
         title="Relative abundance of subcellular marker proteins",point_color="", point_fill="") +
    theme_bw() + theme(legend.position = "none")
  
  fin_med_ratio <- fin_med_ratio %>% dplyr::rename("Gene"="Genes")
  fin <- fin %>% dplyr::rename("Gene"="Genes")
  
  dat_fin <- list(fin_med_ratio = fin_med_ratio, fin = fin, plot_sepConds=plot_sepConds, plot_CombinedConds_box=plot_CombinedConds_box,
                  plot_Conds_density=plot_Conds_density)
  return(dat_fin)
}

pD_HY_keepPrecs <- function(dat,compression_factor=6){
  dat2 <- dat %>% dplyr::mutate("Species"=ifelse(grepl("HUMAN",Protein.Names),"H. sapiens",
                                                 "S. cerevisiae"))
  dat2 <- pD_rmMixSpec(dat2)
  dat2 <- dat2 %>% dplyr::mutate("log2MS1" = log2(Ms1.Area)) %>% filter(log2MS1>0) 
  dat2 <- dat2 %>% dplyr::mutate("Celltype_seqcharge_prot" = paste0(Celltype,"_",seqcharge,"_",Protein.Names))
  dat04 <- dat2 %>% filter(grepl("mTRAQ0|mTRAQ4",channel_name)) %>% dplyr::group_by(Celltype,seqcharge) %>% dplyr::mutate("norm_med_single"=log2(median(2^log2MS1,na.rm=T))) %>% dplyr::ungroup() %>%
    ungroup() %>% dplyr::distinct(Celltype,seqcharge,.keep_all=T) %>% dplyr::select("Celltype_seqcharge_prot","norm_med_single") 
  dat8 <- dat2 %>% filter(grepl("mTRAQ8",channel_name)) %>% dplyr::group_by(Celltype,seqcharge) %>% dplyr::mutate("norm_med"=log2(median(2^log2MS1,na.rm=T))) %>% dplyr::ungroup() %>%
    ungroup() %>% dplyr::distinct(Celltype,seqcharge,.keep_all=T)
  dat_d048 <- dat04 %>% inner_join(dat8, by =c("Celltype_seqcharge_prot"="Celltype_seqcharge_prot"))
  dat_d048 <- dat_d048 %>% dplyr::mutate("SC_car_FC" = norm_med_single-norm_med) 
  dat_d048$car <- as.numeric(str_sub(dat_d048$Celltype, end = -2)) 
  dat_d048_precs <- dat_d048 %>% filter(Celltype=="50x") %>% dplyr::mutate(dif = (2^SC_car_FC)*(car))
  
  compression_plot <- ggplot(dat_d048_precs, aes(x=(SC_car_FC))) + geom_histogram() +
    annotate("rect", xmin=-Inf, xmax=log2(compression_factor*1/50), ymin=-Inf, ymax=Inf, color="gray10", fill="green",alpha=0.2) +
    geom_vline(xintercept=log2(1/50),color="red",linewidth=1) +
    labs(x="Log2, median(Single nuc / 50x Carrier) precusors", y= "Precursors quantified") +
    theme_classic()
  dat_d048_precsKeep <- dat_d048_precs %>% filter(log2(dif)<log2(compression_factor))
  keep_precs <- unique(dat_d048_precsKeep$seqcharge)
  
  fin <- list(keep_precs=keep_precs,compression_plot=compression_plot)
  return(fin)
}

pD_rmMixSpec <- function(df){
  df$HY<-F
  df$HY[grepl("HUMAN",df$Protein.Names) & grepl("YEAST",df$Protein.Names) & grepl("ECOLI",df$Protein.Names)] <-T
  df$HY[grepl("HUMAN",df$Protein.Names) & grepl("ECOLI",df$Protein.Names)] <-T
  df$HY[grepl("HUMAN",df$Protein.Names) & grepl("YEAST",df$Protein.Names)] <-T
  df$HY[grepl("ECOLI",df$Protein.Names) & grepl("YEAST",df$Protein.Names)] <-T
  table(df$HY)
  df <- df[grepl("FALSE", df$HY),] %>% dplyr::select(-"HY")
  
  return(df)
}

pD_HY_filt <- function(dat, keep_precs){
  dat <- dat[dat$seqcharge%in%keep_precs,]
  dat2 <- dat %>% dplyr::mutate("Species"=ifelse(grepl("HUMAN",Protein.Names),"H. sapiens",
                                                 "S. cerevisiae"))
  dat2 <- pD_rmMixSpec(dat2)
  #count
  dat2count <- dat2 %>% group_by(run_chan) %>% add_count() %>% 
    ungroup() %>% distinct(run_chan, .keep_all=T) %>% dplyr::select("run_chan","Celltype", "channel_name", "n")
  dat2count <- na.omit(dat2count)
  dat2count$Celltype <- factor(dat2count$Celltype, levels=c("0x", "1x", "5x","10x","25x","50x"))
  prec_count <- ggplot(dat2count %>% filter(!grepl("mTRAQ8",channel_name)), 
                       aes(x=Celltype,y=n,color=channel_name)) + geom_beeswarm(size=2,alpha=0.8) + theme_bw() +
    ylim(0,25000) + 
    labs(y="Precursors quantified in single-nuclei channels",x="Carrier amount",title="")
  
  means <- dat2 %>% dplyr::group_by(channel_name,Celltype,seqcharge) %>% filter(Ms1.Area_iso>0) %>% 
    dplyr::summarize("mean_prec_ab" = mean(Ms1.Area_iso,na.rm=T)) %>% ungroup() %>% group_by(seqcharge) %>%
    add_count() %>% filter(n==6*3) %>% ungroup() 
  
  means$Celltype <- factor(means$Celltype, levels=c("0x", "1x", "5x","10x","25x","50x"))
  abundances_plot <- ggplot(means, aes(x=log2(mean_prec_ab),y=Celltype,fill=channel_name)) + geom_density_ridges2(alpha=0.5) +
    theme_bw() + labs(y="Carrier amount", x = "Log2, precursor abundance")
  
  means_sumary <- means %>% filter(channel_name=="mTRAQ8") %>% group_by(Celltype) %>%
    dplyr::summarize(median(mean_prec_ab))
  
  dat2$log2MS1 <- log2(dat2$Ms1.Area_iso)
  dat2 <- dat2[dat2$log2MS1>0,]
  dat2 <- dat2 %>% dplyr::group_by(run_chan) %>% dplyr::mutate("norm"=log2MS1-median(log2MS1,na.rm=T)) %>% dplyr::ungroup()
  dat2 <- dat2 %>% dplyr::group_by(channel_name,Celltype,seqcharge) %>% dplyr::mutate("norm_med"=median(norm,na.rm=T)) %>% dplyr::ungroup() %>%
    dplyr::distinct(channel_name,Celltype,seqcharge,.keep_all=T)
  
  dat2$Celltype_seqcharge_prot <- paste0(dat2$Celltype,"_",dat2$seqcharge,"_",dat2$Protein.Names)
  dat2_d0 <- dat2[grepl("mTRAQ0",dat2$channel_name),] %>% dplyr::select("Celltype_seqcharge_prot","norm_med") %>% dplyr::rename("norm_med_d0"="norm_med")
  dat2_d4 <- dat2[grepl("mTRAQ4",dat2$channel_name),]
  
  dat2_d04 <- dat2_d0 %>% inner_join(dat2_d4, by =c("Celltype_seqcharge_prot"="Celltype_seqcharge_prot"))
  
  
  dat2_d04$log2FC <- dat2_d04$norm_med_d0-dat2_d04$norm_med
  dat2_d04$Celltype <- factor(dat2_d04$Celltype, levels=c("0x", "1x", "5x","10x","25x","50x"))
  dat2_d04_H <- dat2_d04 %>% filter(Species=="H. sapiens") %>% group_by(Celltype) %>%
    dplyr::summarize("med_H" = median(log2FC,na.rm=T)) %>% ungroup()
  dat2_d04 <- dat2_d04 %>% left_join(dat2_d04_H, by =c("Celltype" = "Celltype"))
  dat2_d04$log2FC <- dat2_d04$log2FC-dat2_d04$med_H
  prec_all <-  ggplot(dat2_d04, aes(x=Celltype,y=log2FC,fill=Species)) + geom_boxplot() + theme_bw() +
    geom_hline(yintercept=2, linetype="dashed",color="blue")+
    labs(y="Log2, d0/d4",x="Carrier amount",title="Precursor-level, all")
  
  dat2_d04_int <- dat2_d04 %>% dplyr::select(-c("med_H")) %>% add_count(seqcharge) %>% filter(n==6)
  dat2_d04_intH <- dat2_d04_int %>% filter(Species=="H. sapiens") %>% group_by(Celltype) %>%
    dplyr::summarize("med_H" = median(log2FC,na.rm=T)) %>% ungroup()
  dat2_d04_int <- dat2_d04_int %>% left_join(dat2_d04_intH, by =c("Celltype" = "Celltype"))
  dat2_d04_int$log2FC <- dat2_d04_int$log2FC-dat2_d04_int$med_H
  prec_int <-  ggplot(dat2_d04_int, aes(x=Celltype,y=log2FC,fill=Species)) + geom_boxplot(width=0.5)+ theme_bw()+
    geom_hline(yintercept=2, linetype="dashed",color="blue") + 
    labs(y="Log2, d0/d4",x="Carrier amount",title="Precursor-level, intersected")
  
  meta <- dat %>% distinct(run_chan,.keep_all=T) %>% dplyr::select("run_chan","Celltype","channel_name")
  
  
  ############ maxLFQ
  dat_maxLFQ <- diann_maxlfq(dat, group.header="Protein.Names", id.header = "seqcharge", quantity.header = "Ms1.Area_iso")
  dat_maxLFQ <- data.frame(dat_maxLFQ)
  dat_maxLFQ$Protein.Names <- row.names(dat_maxLFQ)
  dat_maxLFQ_m <- melt(dat_maxLFQ)
  dat_maxLFQ_m <- dat_maxLFQ_m %>% dplyr::mutate("Species"=ifelse(grepl("HUMAN",Protein.Names),"H. sapiens",
                                                                  "S. cerevisiae"))
  dat_maxLFQ_m <- pD_rmMixSpec(dat_maxLFQ_m)
  
  dat_maxLFQ_m$log2MS1 <- log2(dat_maxLFQ_m$value)
  dat_maxLFQ_m <- dat_maxLFQ_m[dat_maxLFQ_m$log2MS1>0,]
  dat_maxLFQ_m <- dat_maxLFQ_m %>% dplyr::group_by(variable) %>% dplyr::mutate("norm"=log2MS1-median(log2MS1,na.rm=T)) %>% dplyr::ungroup()
  dat_maxLFQ_m$variable <- gsub(".", "-", dat_maxLFQ_m$variable, fixed = TRUE) #in case column names changed hyphen to dot
  dat_maxLFQ_m <- dat_maxLFQ_m %>% left_join(meta, by=c("variable"="run_chan"))
  #count
  dat_maxLFQ_m_count <- dat_maxLFQ_m %>% group_by(variable) %>% add_count() %>% 
    ungroup() %>% distinct(variable, .keep_all=T) %>% dplyr::select("variable","Celltype", "channel_name", "n")
  dat_maxLFQ_m_count <- na.omit(dat_maxLFQ_m_count)
  dat_maxLFQ_m_count$Celltype <- factor(dat_maxLFQ_m_count$Celltype, levels=c("0x", "1x", "5x","10x","25x","50x"))
  prot_count <- ggplot(dat_maxLFQ_m_count %>% filter(!grepl("mTRAQ8",channel_name)), 
                       aes(x=Celltype,y=n,color=channel_name)) + geom_beeswarm(size=2,alpha=0.8) + theme_bw() +
    ylim(0,3000) + 
    labs(y="Proteins quantified in single-nuclei channels",x="Carrier amount",title="")
  
  
  dat_maxLFQ_m <- dat_maxLFQ_m %>% dplyr::group_by(channel_name,Celltype,Protein.Names) %>% 
    dplyr::mutate("norm_med"=median(norm,na.rm=T)) %>% dplyr::ungroup() %>%
    dplyr::distinct(channel_name,Celltype,Protein.Names,.keep_all=T)
  dat_maxLFQ_m$car_Prot <- paste0(dat_maxLFQ_m$Protein.Names,"_",dat_maxLFQ_m$Celltype)

  ###
  dat2_prot_d0 <- dat_maxLFQ_m[grepl("mTRAQ0",dat_maxLFQ_m$channel_name),] %>% dplyr::select("car_Prot","norm_med") %>% dplyr::rename("norm_med_d0"="norm_med")
  dat2_prot_d4 <- dat_maxLFQ_m[grepl("mTRAQ4",dat_maxLFQ_m$channel_name),]
  
  dat2_prot_d04 <- dat2_prot_d0 %>% inner_join(dat2_prot_d4, by =c("car_Prot"="car_Prot"))
  
  dat2_prot_d04$log2FC <- dat2_prot_d04$norm_med_d0-dat2_prot_d04$norm_med
  dat2_prot_d04$Celltype <- factor(dat2_prot_d04$Celltype, levels=c("0x", "1x", "5x","10x","25x","50x"))
  
  dat2_prot_d04_H <- dat2_prot_d04 %>% filter(Species=="H. sapiens") %>% group_by(Celltype) %>%
    dplyr::summarize("med_H" = median(log2FC,na.rm=T)) %>% ungroup()
  dat2_prot_d04 <- dat2_prot_d04 %>% left_join(dat2_prot_d04_H, by =c("Celltype" = "Celltype"))
  dat2_prot_d04$log2FC <- dat2_prot_d04$log2FC-dat2_prot_d04$med_H
  prot_all <-  ggplot(dat2_prot_d04, aes(x=Celltype,y=log2FC,fill=Species)) + geom_boxplot(width=0.5) + theme_bw() +
    geom_hline(yintercept=2, linetype="dashed",color="blue")+
    labs(y="Log2, d0/d4",x="Carrier amount",title="Protein-level, all")
  
  dat2_prot_d04int <- dat2_prot_d04 %>% dplyr::select(-c("med_H")) %>% add_count(Protein.Names) %>% filter(n==6)
  dat2_prot_d04int_H <- dat2_prot_d04int %>% filter(Species=="H. sapiens") %>% group_by(Celltype) %>%
    dplyr::summarize("med_H" = median(log2FC,na.rm=T)) %>% ungroup()
  dat2_prot_d04int <- dat2_prot_d04int %>% left_join(dat2_prot_d04int_H, by =c("Celltype" = "Celltype"))
  dat2_prot_d04int$log2FC <- dat2_prot_d04int$log2FC-dat2_prot_d04int$med_H
  prot_int <-  ggplot(dat2_prot_d04int, aes(x=Celltype,y=log2FC,fill=Species)) + geom_boxplot(width=0.5)+ theme_bw()+
    geom_hline(yintercept=2, linetype="dashed",color="blue") + 
    labs(y="Log2, d0/d4",x="Carrier amount",title="Protein-level, intersected") 
  
  fin <- list(dat2_prot_d04=dat2_prot_d04,abundances_plot=abundances_plot, prec_count=prec_count, prec_all=prec_all, prec_int=prec_int, prot_count=prot_count, prot_all=prot_all, prot_int=prot_int)
  return(fin)
}

pD_compare_HY_filt <- function(HY_all, HY_filtered, HY_filtered2, HY_filtered3){
  all <- HY_all$dat2_prot_d04
  all$type <- "None"
  filt <- HY_filtered$dat2_prot_d04
  filt$type <- "2x (heavy filt)"
  filt1 <- HY_filtered2$dat2_prot_d04
  filt1$type <- "4x (mod. filt)"
  filt2 <- HY_filtered3$dat2_prot_d04
  filt2$type <- "6x (light filt)"
  both <- rbind(all, filt,filt1,filt2) %>% filter(Celltype=="50x")
  both_filt <- both %>% add_count(car_Prot) %>% filter(n==4)
  upper_quant <- quantile(both_filt$log2FC,0.999)
  lower_quant <- quantile(both_filt$log2FC,0.001)
  
  prot_int <-  ggplot(both_filt, aes(x=type,y=log2FC,fill=Species)) + geom_boxplot(width=0.5)+ theme_classic()+
    labs(y="Log2, d0/d4",x="Compression filtering",title="Protein-level, intersected") +
    geom_hline(yintercept=2, linetype="dashed",color="#e56730",linewidth=0.75)+
    geom_hline(yintercept=0, linetype="dashed",color="#109abf",linewidth=0.75)+
    scale_fill_manual(values=c("#68aabc","#e08a6c"))

  upper_quant <- quantile(both$log2FC,0.999)
  lower_quant <- quantile(both$log2FC,0.001)
  prot_all <-  ggplot(both, aes(x=type,y=log2FC,fill=Species)) + geom_boxplot(width=0.5)+ theme_classic()+
    labs(y="Log2, d0/d4",x="Compression filtering",title="Protein-level, all")+
    geom_hline(yintercept=2, linetype="dashed",color="#e56730",linewidth=0.75)+
    geom_hline(yintercept=0, linetype="dashed",color="#109abf",linewidth=0.75)+
    scale_fill_manual(values=c("#68aabc","#e08a6c"))

  prot_counts <- both %>% add_count(type) %>% distinct(type,.keep_all=T) %>% dplyr::select("type","n")
  prot_count <-  ggplot(prot_counts, aes(x=type,y=n)) + geom_bar(stat="identity")+ theme_classic()+
    labs(y="Proteins quantified",x="Compression filtering",title="")
  fin <- list(prot_int=prot_int,prot_all=prot_all,prot_count=prot_count)
  return(fin)
}

pD_filter_proteotypic <- function(dat){
  dat <- dat[dat$Proteotypic==T,]
  return(dat)
}

pD_filter_count <- function(data,TQV_threshold=0.1,CQV_threshold=0.1,min_count=50){
  all_groups <- data %>% 
    distinct(run_chan,.keep_all=T) %>% 
    dplyr::mutate(count = 0) #initialize counts to 0
  
  filtered_counts <- data %>% 
    filter(Translated.Q.Value < TQV_threshold & Channel.Q.Value < CQV_threshold & Ms1.Area>0) %>%
    group_by(run_chan) %>% 
    dplyr::summarise(count = n(), .groups = 'drop') 
  
  final_counts <- all_groups %>%
    left_join(filtered_counts, by = "run_chan", suffix = c(".all_groups", ".filtered")) %>%
    dplyr::mutate(count = coalesce(count.filtered, count.all_groups)) %>% ungroup() %>%
    select(run_chan, count,Celltype) %>% dplyr::mutate("real"=ifelse(grepl("neg",Celltype),"Negative control","Nucleus"))
  retained <- final_counts[final_counts$count>=min_count,] %>% filter(!grepl("mTRAQ8",run_chan)) %>% filter(Celltype!="neg") %>% nrow(.)
  total <- final_counts %>% filter(!grepl("mTRAQ8",run_chan)) %>% filter(Celltype!="neg") %>% nrow(.)
  
  coverage_plot <- ggplot(final_counts %>% filter(!grepl("mTRAQ8",run_chan)), aes(x = count, fill = real)) +
    labs(x = "Precursors with high quality", color = "", y = "Nuclei") +
    geom_rect(inherit.aes = FALSE, 
              mapping = aes(xmin = 50, xmax = Inf, ymin = 0, ymax = Inf, fill = NULL), 
              color = "gray10", fill = "green") + 
    geom_histogram() +
    theme_bw() + scale_fill_manual(values=c("red2","grey30"))+
    geom_text(x=900, y=330, label=paste0(retained," retained, ",total," total")) +
    theme(legend.position = c(0.8, 0.85), legend.text = element_text(size = 12),
          legend.background = element_rect(fill = "white", size = 0.2, linetype = "solid", colour = "gray40"),
          legend.title = element_blank(),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"))
  
  final_counts_filt <- final_counts %>% filter(count>=min_count) %>% filter(!grepl("mTRAQ8",run_chan))
  return(list(final_counts_filt=final_counts_filt, coverage_plot=coverage_plot))
}

pD_median_car <- function(dat, car="mTRAQ8"){
  dat1 <- dat %>% dplyr::group_by(seqcharge) %>% dplyr::filter(channel_name==car) %>%
    dplyr::mutate(mean_car= mean(Ms1.Area)) %>% distinct(seqcharge,.keep_all=T) %>%
    dplyr::select(seqcharge, mean_car)
  return(dat1)
}

pD_filterCells <- function(data_noFilt,keep_cells){
  dat <- data_noFilt[data_noFilt$run_chan%in%keep_cells$final_counts_filt$run_chan,]
  dat <- dat[dat$Celltype!="neg",] #remove neg controls if any
  datcar <- data_noFilt[data_noFilt$Celltype=="Carrier",]
  both <- rbind(dat,datcar)
  return(both)
}

pD_car_norm_data <- function(dat, Run ="Run", seqcharge ="seqcharge", Celltype="Celltype", Quant="Ms1.Area"){
  dat$pep_run <- paste0(dat$Run,dat$seqcharge)
  datcar <- dat[grepl("Carrier", dat$Celltype),] %>% dplyr::select("pep_run",all_of(Quant)) %>%
    dplyr::rename("Carrier_MS1"=all_of(Quant)) %>% distinct(pep_run,.keep_all=T)
  dat <- dat[!grepl("Carrier", dat$Celltype),] %>% distinct(seqcharge,run_chan,.keep_all=T)
  dat_fin <- dat %>% left_join(datcar, by =c("pep_run"="pep_run"))
  dat_fin$pep_rat <- dat_fin[[Quant]]/dat_fin$Carrier_MS1
  dat_fin <- dat_fin[which(dat_fin$pep_rat>0 & dat_fin$pep_rat<Inf),]
  return(dat_fin)
}

pD_rm_poorQuant <- function(dat, min_proportion=0.1, compression_factor=6){
  n_uniruns <- length(unique(dat$run_chan)) #number of single nuc 
  precs1 <- dat %>% dplyr::group_by(seqcharge) %>% add_count() %>%
    dplyr::mutate("med_rat" = median(pep_rat*Carrier_amount)) %>% #account for carrier amount
    dplyr::mutate(n_proportion = n/n_uniruns)%>% distinct(seqcharge, .keep_all = T) %>%
    dplyr::select(seqcharge,med_rat,n_proportion)
  
  prec_plot <- ggplot(precs1, aes(x=log10(med_rat),y=n_proportion)) +
    annotate("rect", xmin=-Inf, xmax=log10(compression_factor), ymin=min_proportion, ymax=Inf, color="gray10", fill="green") +
    geom_point(alpha=0.07) +
    geom_vline(xintercept=0,color="red",linewidth=1) +
    labs(x="Log10, median(Single nuc / Carrier) precusors", y= "Fraction of nuclei quantified") +
    theme_bw()
  
  dat <- dat %>% dplyr::group_by(seqcharge) %>% add_count() %>%
    dplyr::mutate("med_compress" = log10(median(pep_rat*Carrier_amount))) %>% dplyr::ungroup() %>%
    dplyr::filter((med_compress)<log10(compression_factor)) %>% dplyr::mutate(n_proportion = n/n_uniruns) %>%
    dplyr::filter(n_proportion>min_proportion) %>% dplyr::select(-c("n"))
  
  dat_list <- list(dat = dat, prec_plot = prec_plot)
  
  return(dat_list)
}

pD_SN_maxLFQPrep <- function(dat, car){
  dat1 <- dat %>% left_join(car, by =c("seqcharge"="seqcharge")) %>% 
    dplyr::mutate("val" = pep_rat*mean_car)
  return(dat1)
}

pD_intersect_seqcharge <- function(dat, keep){
  keep_temp <- keep$dat
  dat <- dat[dat$seqcharge%in%keep_temp$seqcharge,]
  dat$File.Name <- dat$run_chan
  return(dat)
}

pD_nuc_align <- function(df,dat_bulk_purity){
  #by having a known proportion of WC:nucleus for each protein, and taking a ratio of single nucleus:reference nucleus, we can quantify enrichment as a 3rd degree polynomial
  NT_WCNUC <- dat_bulk_purity$fin_med_ratio
  dis_prot_all<- df %>% inner_join(NT_WCNUC, by =c("Genes" = "Gene"))
  dis_prot_all <- dis_prot_all[!is.na(dis_prot_all$norm_prot),]
  dis_prot_all <- dis_prot_all[!is.na(dis_prot_all$med_WC_NUC_log2_rat),]
  unirunchans <- unique(dis_prot_all$run_chan)
  dis_prot_all <- dis_prot_all %>% dplyr::mutate(norm_bulk_NUC_log2 = med_NUC_log2-median(med_NUC_log2,na.rm=T))
  dis_prot_all <- dis_prot_all %>% dplyr::group_by(run_chan) %>% dplyr::mutate(med_norm_prot = norm_prot-median(norm_prot,na.rm=T)) %>% ungroup()
  
  dis_prot_all$SN_vs_BN <- dis_prot_all$med_norm_prot - dis_prot_all$norm_bulk_NUC_log2
  list_dfs <- list()
  #fit the data
  for(i in 1:length(unirunchans)){
    uni1 <- dis_prot_all[dis_prot_all$run_chan==unirunchans[i],]
    temp <- uni1 %>% arrange(med_WC_NUC_log2_rat)
    model <- rlm((temp$SN_vs_BN) ~ poly(temp$med_WC_NUC_log2_rat,3)) #robust 3rd degree polynomial
    temp$predlm = predict(model,data.frame((temp$SN_vs_BN)))
    temp$prot_norm_adj <- temp$SN_vs_BN-temp$predlm #if we wanted to correct for the differences in purity (probably quite noisy, I don't recommend using this)
    list_dfs[[i]] <- temp
  }
  df <- do.call("rbind", list_dfs)
  df <- df %>% dplyr::rename("WC_nuc" = "med_WC_NUC_log2_rat")
  return(df)
}

pD_AUC_nucPurity <- function(df,maxAUC=5){
  #calculate area under the curve as a way of quantifying lack of nuclear enrichment
  ints <- df %>% dplyr::group_by(run_chan) %>% dplyr::filter(WC_nuc>-2&WC_nuc<7.5) %>%
    dplyr::mutate("abs_AUC" = flux::auc(WC_nuc,abs(predlm))) %>% ungroup() 
  
  ints_dis <- ints %>% distinct(run_chan, .keep_all=T)
  maxpredlm <- max(ints$predlm)
  minpredlm <- min(ints$predlm)
  
  AUC_plot <- ggplot(ints, aes(x=WC_nuc, y=predlm, group=run_chan, color=abs_AUC)) + geom_line(alpha=0.5) +
    scale_color_viridis(limits=c(0,max(ints$abs_AUC))) + theme_bw() + geom_abline(slope=1,color="red",linetype="dashed")+
    labs(title=paste0(length(unique(ints$run_chan))," nuclei"),
         x="Log2, Reference whole cell / Reference nuclei", 
         y="Log2, Single nuclei / Reference nuclei") + ylim(minpredlm,maxpredlm)
  
  ints_filt <- ints[ints$abs_AUC<maxAUC,]
  AUC_plot_filt <- ggplot(ints_filt, aes(x=WC_nuc, y=predlm, group=run_chan, color=abs_AUC)) + geom_line(alpha=0.5) +
    scale_color_viridis(limits=c(0,max(ints$abs_AUC))) + theme_bw() + geom_abline(slope=1,color="red",linetype="dashed")+
    labs(title=paste0(length(unique(ints_filt$run_chan))," nuclei"),
         x="Log2, Reference whole cell / Reference nuclei", 
         y="Log2, Single nuclei / Reference nuclei") + ylim(minpredlm,maxpredlm)
  
  plot_AUC <-  ggplot(ints_dis, aes(x=abs_AUC)) + 
    annotate("rect", xmin=-Inf, xmax=maxAUC, ymin=-Inf, ymax=Inf, color="gray10", fill="green") +
    geom_histogram()+
    labs(x="AUC, nuclear impurity", y="Count of nuclei") +
    theme_bw()

  fin <- list(ints=ints, AUC_plot=AUC_plot,AUC_plot_filt=AUC_plot_filt,plot_AUC=plot_AUC)
  return(fin)
}

pD_filter_purity <- function(df, purity, maxAUC=5){
  keeps <- purity$ints[purity$ints$abs_AUC<maxAUC,]
  df <- df[df$run_chan%in%keeps$run_chan,]
  return(df)
}

pD_prot_norm <- function(df, Quant = "norm_prot"){
  df$Quantnew <- df[[Quant]]
  dat_fin <- df %>% dplyr::group_by(run_chan) %>%
    dplyr::mutate("norm_prot" = Quantnew - median(Quantnew, na.rm=T)) %>%
    ungroup() %>% dplyr::group_by(Genes) %>%
    dplyr::mutate("norm_prot" = norm_prot - mean(norm_prot, na.rm=T)) %>% ungroup()
  return(dat_fin)
}

pD_impute_log2 <- function(df, missing.prot.frac=0.7, missing.cell.frac=0.95, k=10){
  dat_fin1 <- df %>% dplyr::select("Genes", "run_chan", "norm_prot")
  dat_fin1$Genes <- as.factor(dat_fin1$Genes)
  dat_fin1$norm_prot <- as.numeric(dat_fin1$norm_prot)
  dat_fin1$run_chan <- gsub('\\.', '-', dat_fin1$run_chan) 
  dat_fin1$run_chan <- as.factor(dat_fin1$run_chan)
  dat_fin2<-dcast(dat_fin1, Genes~run_chan, value.var = "norm_prot", fill=NA)
  dat3<-as.matrix(dat_fin2[,-1]); row.names(dat3)<-dat_fin2[,1]
  
  # impute single celldata
  pre_impute <-(filt.mat.rc(dat3,missing.prot.frac,missing.cell.frac))
  sc.imp<-cr_norm_log(pre_impute)
  sc.imp <- hknn_weighted_EucNorm(sc.imp, k)
  sc.imp<-cr_norm_log((sc.imp))
  
  dat <- list(pre_impute = pre_impute, sc.imp = sc.imp)
  
  return(dat)
}

filt.mat.rc <- function(mat, pct.r,pct.c){
  #from H. Specht SCoPE2
  kr<-c()
  for(k in 1:nrow(mat)){
    pct.na<-length(which(is.na(mat[k,]))) / length(mat[k,])
    if(pct.na <= pct.r){ kr<-c(kr,k)}
  }
  mat<-mat[kr,]
  kc<-c()
  for(k in 1:ncol(mat)){
    pct.na<-length(which(is.na(mat[,k]))) / length(mat[,k])
    if(pct.na <= pct.c){ kc<-c(kc,k)}
  }
  mat<-mat[,kc]
  return(mat)
}

cr_norm_log<-function(dat){
  #normalize proteins in each sample by median protein abundance in each sample
  for(k in 1:ncol(dat)){
    dat[,k]<-dat[,k]-median(dat[,k], na.rm = T)
  }
  
  #normalize each protein by mean across samples
  for(k in 1:nrow(dat)){
    dat[k,]<-dat[k,]-mean(dat[k,], na.rm = T)
  }
  
  return(dat)
}

hknn_weighted_EucNorm <-function(dat, k){
  
  # create a copy of the data, NA values to be filled in later
  dat.imp<-dat
  dist.mat<-as.matrix( dist(t(dat)) )
  counts <- pairwiseCount(dat)
  dist.mat <- dist.mat/counts
  cnames<-colnames(dist.mat)
  for(X in cnames){
    
    # find the distances of all other columns to that column 
    distances<-dist.mat[, X]
    distances.ordered<-distances[order(distances, decreasing = F)]
    dat.reordered<-dat[ , names(distances.ordered ) ]
    vec<-dat[, X]
    na.index<-which( is.na(vec) )
    
    for(i in na.index){
  
      # find most similar columns that have a non-NA value in this row
      closest.columns<-names( which( !is.na(dat.reordered[i, ])  ) )
      distances.ordered_mat <- as.matrix(distances.ordered)
      
      #from the samples that are ordered in terms of their similarity, select the ones which have a non-NA value
      distances.ordered_mat <- distances.ordered_mat[rownames(distances.ordered_mat)%in%closest.columns,] 
      closest_vec <- distances.ordered_mat[1:(k)] #select the first 'k' samples with their euclidian distances
      closest_vec_inv <- 1/closest_vec #inverted euclidian distances (small distances = more similar = should have greater weight)
      total <- sum(closest_vec_inv) #sum the inverted euclidian distances
      
      if( length(closest.columns)>k ){
        # replace NA in column X with the weighted average of the same row in k of the most similar columns
        vec[i]<-sum(closest_vec_inv*as.matrix(dat[ i, closest.columns[1:k] ] ))/total
      }
      
      if( length(closest.columns)<=k ){
        j <- length(closest.columns)
        closest_vec <- distances.ordered_mat[1:(j)]
        closest_vec_inv <- 1/closest_vec
        total <- sum(closest_vec_inv)
        # replace NA in column X with the weighted avg of the same row in all of the most similar columns
        vec[i]<-sum(closest_vec_inv*as.matrix(dat[ i, closest.columns[1:j] ] ))/total
      }
    }
    # replace each column with the new imputed values column
    dat.imp[,X]<-vec
  }
  return(dat.imp)
  
}

pD_batchCorrect_labs_LC <- function(dat, meta_fpath){
  sc.imp <- dat$sc.imp
  m <- fread(meta_fpath)
  m$run_chan <- paste0(m$Raw, "mTRAQ",m$Label)
  no_neg <- m[!grepl("neg|Carrier",m$Celltype),]
  no_neg$lab_set <- paste0(no_neg$Label,no_neg$LC_batch)
  uni_BR <- unique(no_neg$BioRep)
  temp_list <- list()
  
  # correct for labels and LC-batch before biological replicate (otherwise potentially too many batches to fit the data well)
  for(i in 1:length(uni_BR)){
    no_neg_temp <- no_neg[no_neg$BioRep==uni_BR[i],]
    temp <- sc.imp[,colnames(sc.imp)%in%no_neg_temp$run_chan]
    batch.covs <-no_neg_temp$lab_set[match(colnames(temp), no_neg_temp$run_chan)]
    lab_set_counts <- table(batch.covs)  # count the number of cells in each lab_set
    valid_lab_sets <- names(lab_set_counts[lab_set_counts >= 10]) # only keep lab_set groups that have 10 or more cells
    batch.covs <- batch.covs[batch.covs%in%valid_lab_sets]
    no_neg_temp <- no_neg_temp[no_neg_temp$lab_set%in%valid_lab_sets]
    temp <- temp[,colnames(temp)%in%no_neg_temp$run_chan] # only keep those cell batches with >10 
    temp<-cr_norm_log((temp))
    batch.covs <-no_neg_temp$lab_set[match(colnames(temp), no_neg_temp$run_chan)]
    mod<-data.frame(no_neg_temp$Celltype[match(colnames(temp), no_neg_temp$run_chan)]); colnames(mod)<-"celltype"
    mod<-model.matrix(~as.factor(celltype), data=mod)
    matrix.sc.batch_temp <- ComBat(temp, batch=batch.covs, mod=mod)
    matrix.sc.batch_temp<-cr_norm_log((matrix.sc.batch_temp))
    temp_list[[i]] <- matrix.sc.batch_temp
  }
  
  # correct for biological replicates if there's >1
  if(length(uni_BR>1)){
    matrix.sc.batch <- do.call("cbind",temp_list)
    imputed.BC<-cr_norm_log((matrix.sc.batch))
    batch.covs <-no_neg$BioRep[match(colnames(imputed.BC), no_neg$run_chan)]
    mod<-data.frame(no_neg$Celltype[match(colnames(imputed.BC), no_neg$run_chan)]); colnames(mod)<-"celltype"
    mod<-model.matrix(~as.factor(celltype), data=mod)
    imputed.BC_updated <- ComBat(imputed.BC, batch=batch.covs, mod=mod)
    imputed.BC_updated<-cr_norm_log((imputed.BC_updated))
  } else {
    matrix.sc.batch <- do.call("cbind",temp_list)
    imputed.BC_updated<-cr_norm_log((matrix.sc.batch))
  }
  
  # data to use:
  mat.sc.imp_NA <- imputed.BC_updated
  unimp <- dat$pre_impute
  unimp <- unimp[,colnames(unimp)%in%no_neg$run_chan]
  unimp <- unimp[,colnames(unimp)%in%colnames(matrix.sc.batch)]
  mat.sc.imp_NA[is.na(unimp)==T] <- NA # replace imputed with NA
  mat.sc.imp_NA <- cr_norm_log(mat.sc.imp_NA)
  unimputed.BC <- melt(mat.sc.imp_NA); colnames(unimputed.BC)<-c("Genes","run_chan","norm_prot")
  
  dat_fin <- list(unimputed.BC = unimputed.BC, imputed.BC_updated = imputed.BC_updated, no_neg=no_neg, mat.sc.imp_NA=mat.sc.imp_NA)
  return(dat_fin)
}

pD_prepost <- function(data_prefiltering, data_postfiltering){
  pre <- data_prefiltering %>% filter(Ms1.Area>0) %>% filter(Celltype%in%c("NUC_NT","NUC_10","NUC_30","NUC_60")) %>%
    distinct(Genes,run_chan,.keep_all=T) %>% add_count(run_chan) %>% distinct(run_chan,.keep_all=T) %>% dplyr::select("run_chan","n")
  pre$type <- "pre"
  
  post <- data_postfiltering %>% filter(!is.na(norm_prot)) %>%
    distinct(Genes, run_chan, .keep_all=T) %>% add_count(run_chan) %>% distinct(run_chan, .keep_all=T) %>% dplyr::select("run_chan","n")
  post$type <- "post"
  both <- rbind(post,pre)
  medians <- both %>%
    group_by(type) %>%
    dplyr::summarize(median_n = median(n, na.rm = TRUE))
  maxval <- max(both$n)
  median_n_pre <- median(pre$n)
  maxval <- max(pre$n)
  ggplot(pre, aes(x=n)) + geom_histogram() + theme_classic() +
    labs(subtitle=paste0(nrow(pre), " nuclei, Median: ", median_n_pre)) +
    geom_vline(aes(xintercept=median_n_pre), color="red", linetype="dashed", size=1) +
    xlim(0,maxval+50)
  ggsave("Figures/prefilteringnuc.pdf",height=2.75,width=4)
  
  median_n_post <- median(post$n)
  ggplot(post, aes(x=n)) + geom_histogram() + theme_classic() +
    labs(subtitle=paste0(nrow(post), " nuclei, Median: ", median_n_post)) +
    geom_vline(aes(xintercept=median_n_post), color="red", linetype="dashed", size=1)+
    xlim(0,maxval+50)
  ggsave("Figures/postfilteringnuc.pdf",height=2.75,width=4)
  
}

pD_weightedBulk_PCA <- function(dat, DA_fpath, mono, quant="prot.p.adj"){
  bulk_dat <- DA_fpath
  bulk_dat$zscore <- p.to.Z(bulk_dat[[quant]]) #convert pval to zscore
  bulk_dat$zscore[is.na(bulk_dat$zscore)] <- max(bulk_dat$zscore,na.rm=T) 
  bulk_dat$zscore[is.infinite(bulk_dat$zscore)] <- max(bulk_dat$zscore,na.rm=T)
  #get zscore from each indiivudal biorep at bulk-level, take the mean across those bioreps and square it
  bulk_dat <- bulk_dat %>% group_by(Gene) %>% add_count() %>%
    dplyr::mutate("mean_z_bulk" = (sum(sign(log2FC)*zscore,na.rm=T))^2/sqrt(n)) %>% ungroup() %>%
    dplyr::distinct(Gene,.keep_all=T) %>% dplyr::select(Gene, mean_z_bulk) %>% filter(!grepl("KRT",Gene))
  
  min_nonzero <- bulk_dat %>% filter(mean_z_bulk>0)
  min_nonzero <- min(min_nonzero$mean_z_bulk)
  bulk_dat$mean_z_bulk[bulk_dat$mean_z_bulk==0] <- min_nonzero #replace 0 pval with lowest non-zero pval
  
  mat.sc.imp <- dat$imputed.BC_updated
  no_neg <- dat$no_neg
  
  # weight proteins based on bulk data
  X.m <- mat.sc.imp
  bulk_dat <- bulk_dat[bulk_dat$Gene%in%row.names(mat.sc.imp),]
  Gene <- row.names(mat.sc.imp)[row.names(mat.sc.imp)%!in%bulk_dat$Gene]
  Gene <- data.frame(Gene)
  bulk_dat <- rbind.fill(bulk_dat,Gene)
  bulk_dat$mean_z_bulk[is.na(bulk_dat$mean_z_bulk)] <- 1 #replace NA with 1
  bulk_dat <- bulk_dat[order(match(bulk_dat$Gene, row.names(mat.sc.imp))),]
  bulk_dat <- diag(bulk_dat$mean_z_bulk)
  mono1 <- mono$df
  mono_dat <- mono1[mono1$Gene%in%row.names(mat.sc.imp),]
  Gene <- row.names(mat.sc.imp)[row.names(mat.sc.imp)%!in%mono_dat$Gene]
  Gene <- data.frame(Gene)
  mono_dat <- rbind.fill(mono_dat,Gene)
  mono_dat$VNR[is.na(mono_dat$VNR)] <- max(mono_dat$VNR,na.rm=T) #replace NA with poor mono
  mono_dat$VNR <- max(mono_dat$VNR)/mono_dat$VNR #larger values were orginally worse
  mono_dat <- mono_dat[order(match(mono_dat$Gene, row.names(mat.sc.imp))),]
  mono_dat <- diag(mono_dat$VNR)
  
  # # Dot product of each protein correlation vector with itself
  r1<-cor(t(X.m))
  rsum<-rowSums(r1^2)
  # Calculate the weighted data matrix:
  X.m1 <- diag(rsum) %*%  X.m  # weight based on protein cors
  X.m2 <- bulk_dat %*%  X.m1 # weight based on bulk
  X.m2 <- mono_dat %*% X.m2 # weight based on continuity of trends (proteins with more continuous trends probably have better signal)
  row.names(X.m2) <- row.names(X.m)
  pca.imp.cor <- cor(X.m2, use="complete.obs")
  # PCA
  sc.pca<-eigen(pca.imp.cor)
  scx<-as.data.frame(sc.pca$vectors)
  colnames(scx)<-paste0("PC",1:ncol(scx))
  scx$cells<-colnames(pca.imp.cor)
  
  ########### FOR PC LOADINGS:
  prot_pca <- cor(t(X.m2), use="complete.obs")
  # PCA
  prot_loadings<-eigen(prot_pca)
  prot_x<-as.data.frame(prot_loadings$vectors)
  colnames(prot_x)<-paste0("PC",1:ncol(prot_x))
  prot_x$prots<-colnames(prot_pca)
  row.names(prot_x)<-colnames(prot_pca)
  PC_loadings <- prot_x
  
  ###########
  
  # percent of variance explained by each principle component
  pca_var <- sc.pca$values
  percent_var<- pca_var/sum(pca_var)*100
  plot(1:length(percent_var), percent_var, xlab="PC", ylab="% of variance explained")
  # map meta data
  pca.melt <- melt(scx); colnames(pca.melt)<-c("id","pc","value")
  # re map ...
  pca.display <- dcast(pca.melt, id ~ pc, value.var = "value", fill=NA)
  pca.display$celltype<-no_neg$Celltype[match(pca.display$id, no_neg$run_chan)]
  pca.display$label<-no_neg$Label[match(pca.display$id, no_neg$run_chan)]
  
  # PC's to display:
  PCx<-"PC1"
  PCy<-"PC2"
  PCz<-"PC3"
  PCza<-"PC4"
  PCzb<-"PC5"
  
  # Display
  pca.display$lab <- paste0("d",pca.display$label)
  pca.display$celltype <- factor(pca.display$celltype, levels=c("NUC_NT", "NUC_10","NUC_30","NUC_60"))
  pca.display2 <- pca.display %>% left_join(no_neg, by =c("id"="run_chan"))
  
  my_fin_colors <- c("violetred","skyblue1","royalblue1", "midnightblue")
  
  pca_1_2<-ggplot(pca.display,aes(x =PC1, y = PC2, color =celltype)) +
    labs(shape="",x = paste0(PCx,"  (", round(percent_var[1],0),"%)"), y = paste0(PCy,"  (", round(percent_var[2],0),"%)"), fill = "LPS treatment") +
    font("ylab",size=20) +
    font("xlab",size=20) +
    font("xy.text", size=15) +
    geom_point(aes(fill=celltype),alpha=0.7,colour="black",pch=21, size=3) +
    theme_classic() +
    theme(legend.position = "none") +
    scale_fill_manual(values = my_fin_colors)+
    scale_color_manual(values = my_fin_colors)+
    annotate("text", x=0.025, y=-0.045, label=paste0(ncol(mat.sc.imp), " nuclei"), size=3) +
    annotate("text", x=0.025, y=-0.05, label=paste0(dim(mat.sc.imp)[1], " proteins"), size=3)
  
  density_PC1_2 <- ggplot(pca.display, aes(x=PC1, fill=celltype)) + geom_density(alpha=0.5) +
    scale_fill_manual(values = my_fin_colors) + theme_classic() + theme(axis.line=element_blank(),
                                                                        axis.ticks = element_blank(),
                                                                        axis.text = element_blank(),
                                                                        axis.title = element_blank(),
                                                                        legend.position = "top",
                                                                        legend.background = element_rect(color = "transparent", fill = "transparent"))
  
  both2 <- grid.arrange(
    density_PC1_2,
    pca_1_2,
    heights = c(0.3, 0.7))  
  
  ggsave("Figures/PC1_2_nuclei.pdf", both2, width = 8.5, height = 8.5)
  
  PC_scores <- pca.display
  variance <- c(PCx, PCy, PCz, PCza, PCzb)
  dat_fin <- list(pca_1_2 = pca_1_2, PC_loadings = PC_loadings, PC_scores = PC_scores, percent_var=percent_var)
  return(dat_fin)
}

pD_TransportScore <- function(nucs_BC, DA_fpath, frac=0.75, quant="prot.p.adj",FC="log2FC"){
  bulk_dat <- DA_fpath
  bulk_dat$zscore <- p.to.Z(bulk_dat[[quant]]) #convert pval to zscore
  bulk_dat$zscore[is.na(bulk_dat$zscore)] <- max(bulk_dat$zscore,na.rm=T)
  bulk_dat$zscore <- sign(bulk_dat[[FC]])*bulk_dat$zscore #convert pval to zscore
  bulk_dat <- bulk_dat %>% dplyr::group_by(Gene) %>% add_count() %>%
    dplyr::mutate("mean_z_bulk" = sum(zscore)/sqrt(n)) %>% ungroup() %>%
    dplyr::distinct(Gene,.keep_all=T) %>% dplyr::select(Gene, mean_z_bulk)
  
  min_nonzero <- min(abs(bulk_dat$mean_z_bulk))
  bulk_dat$mean_z_bulk[bulk_dat$mean_z_bulk==0] <- min_nonzero #replace 0 zscore with lowest non-zero
  datUse <- nucs_BC$unimputed.BC
  no_neg <- nucs_BC$no_neg
  cells_runchan <- no_neg %>% dplyr::select("run_chan","Celltype", "BioRep")
  dat3 <- datUse %>% left_join(cells_runchan, by =c("run_chan" = "run_chan"))
  dat3 <- dat3[!is.na(dat3$norm_prot),]
  dat_treated <- dat3[!grepl("NT",dat3$Celltype),]
  dat_NT <- dat3[grepl("NT",dat3$Celltype),]
  numcells <- length(unique(dat3$run_chan))
  
  prots_perc <- dat3 %>% dplyr::add_count(Genes) %>% dplyr::filter(n/numcells>frac) %>%
    distinct(Genes)
  
  uni_cell <- unique(dat_treated$run_chan)
  NT_runs <- length(unique(dat_NT$run_chan))
  
  list_zscores <- list()
  for(i in 1:length(uni_cell)){
    temp_cell <- dat_treated[dat_treated$run_chan==uni_cell[i],]
    uni_PGs <- as.character(unique(temp_cell$Genes))
    temp_dat_NT <- dat_NT[dat_NT$Genes%in%uni_PGs,] %>% dplyr::add_count(Genes) %>%
      dplyr::filter(n/NT_runs>0.25) %>% dplyr::select(-c("n"))#protein needs to be in >25% of NT nuclei
    uni_PGs <- as.character(unique(temp_dat_NT$Genes))
    temp_cell <- temp_cell[temp_cell$Genes%in%temp_dat_NT$Genes,]
    temps <- rbind(temp_cell, temp_dat_NT)
    temps <- temps %>% group_by(Genes) %>% 
      dplyr::mutate("norm_prot" = norm_prot-mean(norm_prot,na.rm=T)) %>% ungroup() %>%
      dplyr::group_by(Genes) %>% 
      dplyr::mutate(z_score = scale(norm_prot)) %>% ungroup()
    list_zscores[[i]] <- temps[!grepl("NT",temps$Celltype),]
    
  }
  
  ### here we will weight the contribution of proteins to transport score
  # based on the bulk data p-values (which we convert to z-scores).
  # Proteins which are differentially abundant from NT have higher weighted contribution to 
  # the transport score. We also required the transport score to 
  # involve only proteins present in 75% of the unimputed single nuc data.
  df_zscores <- do.call("rbind", list_zscores)
  df_zscores <- df_zscores %>% left_join(bulk_dat, by =c("Genes"="Gene"))
  df_zscores$Zweighted <- df_zscores$z_score*df_zscores$mean_z_bulk
  df_zscores$Zweighted[is.na(df_zscores$Zweighted)] <- df_zscores$z_score[is.na(df_zscores$Zweighted)]
  df_zscores <- df_zscores %>% dplyr::group_by(Genes) %>%
    dplyr::mutate(mean_zscore = mean(z_score)) %>% ungroup()
  df_tran <- df_zscores %>% filter(Genes%in%prots_perc$Genes) %>% dplyr::group_by(run_chan) %>% 
    dplyr::mutate("Trans_score" = mean(Zweighted)) %>% ungroup() %>% 
    distinct(run_chan,.keep_all=T) %>% dplyr::select("run_chan","Trans_score")
  df_zscores <- df_zscores %>% left_join(df_tran, by =c("run_chan"= "run_chan"))
  df_zscores_uni <- df_zscores %>% distinct(run_chan,.keep_all=T)
  
  NT_df_zscores <-  dat_NT %>% group_by(Genes) %>% 
    dplyr::mutate("norm_prot" = norm_prot-mean(norm_prot,na.rm=T)) %>% ungroup() %>%
    dplyr::group_by(Genes) %>% 
    dplyr::mutate(z_score = scale(norm_prot)) %>% ungroup()
  NT_df_zscores <-NT_df_zscores %>% left_join(bulk_dat, by =c("Genes"="Gene"))
  NT_df_zscores$Zweighted <- NT_df_zscores$z_score*NT_df_zscores$mean_z_bulk
  NT_df_zscores$Zweighted[is.na(NT_df_zscores$Zweighted)] <- NT_df_zscores$z_score[is.na(NT_df_zscores$Zweighted)]
  NT_df_zscores <- NT_df_zscores %>% dplyr::group_by(Genes) %>%
    dplyr::mutate(mean_zscore = mean(z_score)) %>% ungroup() 
  NT_df_tran <- NT_df_zscores  %>% filter(Genes%in%prots_perc$Genes) %>% dplyr::group_by(run_chan) %>% 
    dplyr::mutate("Trans_score" = mean(Zweighted)) %>% ungroup() %>% 
    distinct(run_chan,.keep_all=T) %>% dplyr::select("run_chan","Trans_score")
  NT_df_zscores <- NT_df_zscores %>% left_join(NT_df_tran, by =c("run_chan"= "run_chan"))
  
  NT_T_all <- rbind(NT_df_zscores, df_zscores)
  NT_T_all_uni <- NT_T_all %>% distinct(run_chan,.keep_all=T)
  NT_T_all_uni$Celltype <- factor(NT_T_all_uni$Celltype, levels=c("NUC_NT","NUC_10","NUC_30","NUC_60"))
  my_fin_colors <- c("violetred","skyblue1","royalblue1", "midnightblue")
  
  trans_plot <- ggplot(NT_T_all_uni, aes(x=Celltype, y=Trans_score, fill=Celltype)) + 
    geom_boxplot(alpha=0.3,outlier.shape = NA) + geom_jitter(shape=21,alpha=0.3,width=0.15) + 
    labs(x="LPS condition", y= "Transport score") + theme_classic() + 
    theme(legend.position = "none") + 
    scale_fill_manual(values = my_fin_colors)
  ggsave("Figures/TransportScore.pdf",trans_plot,height=3.4,width=2.6)
  
  
  dat_fin <- list(NT_T_all = NT_T_all, list_zscores = list_zscores, trans_plot=trans_plot)
  return(dat_fin)  
  
}

pD_TransportScore_v2_dif <- function(nucs_BC, DA_fpath, frac=0.75, quant="prot.p.adj",FC="log2FC"){
  bulk_dat <- DA_fpath
  sig_bulk <- bulk_dat[bulk_dat[[quant]]<0.05,]
  bulk_dat$zscore <- p.to.Z(bulk_dat[[quant]]) #convert pval to zscore
  bulk_dat$zscore[is.na(bulk_dat$zscore)] <- max(bulk_dat$zscore,na.rm=T)
  bulk_dat$zscore <- sign(bulk_dat[[FC]])*bulk_dat$zscore #convert pval to zscore
  bulk_dat <- bulk_dat %>% dplyr::group_by(Gene) %>% add_count() %>%
    dplyr::mutate("mean_z_bulk" = sum(zscore)/sqrt(n)) %>% ungroup() %>%
    dplyr::distinct(Gene,.keep_all=T) %>% dplyr::select(Gene, mean_z_bulk)
  
  min_nonzero <- bulk_dat %>% filter(mean_z_bulk>0)
  min_nonzero <- min(min_nonzero$mean_z_bulk)
  bulk_dat$mean_z_bulk[bulk_dat$mean_z_bulk==0] <- min_nonzero #replace 0 pval with lowest non-zero
  datUse <- nucs_BC$unimputed.BC
  no_neg <- nucs_BC$no_neg
  cells_runchan <- no_neg %>% dplyr::select("run_chan","Celltype", "BioRep")
  dat3 <- datUse %>% left_join(cells_runchan, by =c("run_chan" = "run_chan"))
  dat3 <- dat3[!is.na(dat3$norm_prot),]
  dat_treated <- dat3[!grepl("NT",dat3$Celltype),]
  dat_NT <- dat3[grepl("NT",dat3$Celltype),]
  numcells <- length(unique(dat3$run_chan))
  
  prots_perc <- dat3 %>% dplyr::add_count(Genes) %>% dplyr::filter(n/numcells>frac) %>%
    distinct(Genes)
  
  uni_cell <- unique(dat_treated$run_chan)
  NT_runs <- length(unique(dat_NT$run_chan))
  
  list_zscores <- list()
  for(i in 1:length(uni_cell)){
    temp_cell <- dat_treated[dat_treated$run_chan==uni_cell[i],]
    uni_PGs <- as.character(unique(temp_cell$Genes))
    temp_dat_NT <- dat_NT[dat_NT$Genes%in%uni_PGs,] %>% dplyr::add_count(Genes) %>%
      dplyr::filter(n/NT_runs>0.25) %>% dplyr::select(-c("n"))#protein needs to be in >25% of NT nuclei
    uni_PGs <- as.character(unique(temp_dat_NT$Genes))
    temp_cell <- temp_cell[temp_cell$Genes%in%temp_dat_NT$Genes,]
    temps <- rbind(temp_cell, temp_dat_NT)
    temps <- temps %>% group_by(Genes) %>% 
      dplyr::mutate("norm_prot" = norm_prot-mean(norm_prot,na.rm=T)) %>% ungroup() %>%
      dplyr::group_by(Genes) %>% 
      dplyr::mutate(z_score = scale(norm_prot)) %>% ungroup()
    list_zscores[[i]] <- temps[!grepl("NT",temps$Celltype),]
    
  }
  
  ### here we will weight the contribution of proteins to translocalization score
  # based on the bulk data p-values (which we convert to z-scores).
  # Proteins which are differentially abundant from NT have higher weighted contribution to 
  # the transport score. We also required the transport score to 
  # involve only proteins present in 75% of the unimputed single nuc data.
  df_zscores <- do.call("rbind", list_zscores)
  df_zscores <- df_zscores %>% left_join(bulk_dat, by =c("Genes"="Gene"))
  df_zscores$Zweighted <- df_zscores$z_score*df_zscores$mean_z_bulk
  df_zscores$Zweighted[is.na(df_zscores$Zweighted)] <- df_zscores$z_score[is.na(df_zscores$Zweighted)]
  df_zscores <- df_zscores %>% dplyr::group_by(Genes) %>%
    dplyr::mutate(mean_zscore = mean(z_score)) %>% ungroup()
  
  SN_dif_ab <- df_zscores %>% group_by(Genes,Celltype) %>% add_count() %>%
    dplyr::mutate("combined_Zscore" = sum(z_score)/sqrt(n)) %>% 
    dplyr::mutate(mean_Zscore = mean(z_score)) %>% ungroup() %>% dplyr::distinct(Genes,Celltype,.keep_all=T) %>%
    dplyr::mutate(p.val = Z.to.p(combined_Zscore)) %>% dplyr::mutate(p.adj = p.adjust(p.val, method = "BH")) %>%
    dplyr::select(c("Genes","Celltype","mean_Zscore","p.adj")) %>% rename("Gene"="Genes") %>% filter(p.adj<0.05)
  all_sig <- unique(c(unique(SN_dif_ab$Gene),unique(sig_bulk$Gene)))
  df_tran <- df_zscores %>% filter(Genes%in%prots_perc$Genes) %>% filter(Genes%in%all_sig) %>% 
    dplyr::group_by(run_chan) %>% 
    dplyr::mutate("Trans_score" = mean(Zweighted)) %>% ungroup() %>% 
    distinct(run_chan,.keep_all=T) %>% dplyr::select("run_chan","Trans_score")
  df_zscores <- df_zscores %>% left_join(df_tran, by =c("run_chan"= "run_chan"))
  df_zscores_uni <- df_zscores %>% distinct(run_chan,.keep_all=T)
  
  NT_df_zscores <-  dat_NT %>% group_by(Genes) %>% 
    dplyr::mutate("norm_prot" = norm_prot-mean(norm_prot,na.rm=T)) %>% ungroup() %>%
    dplyr::group_by(Genes) %>% 
    dplyr::mutate(z_score = scale(norm_prot)) %>% ungroup()
  NT_df_zscores <-NT_df_zscores %>% left_join(bulk_dat, by =c("Genes"="Gene"))
  NT_df_zscores$Zweighted <- NT_df_zscores$z_score*NT_df_zscores$mean_z_bulk
  NT_df_zscores$Zweighted[is.na(NT_df_zscores$Zweighted)] <- NT_df_zscores$z_score[is.na(NT_df_zscores$Zweighted)]
  NT_df_zscores <- NT_df_zscores %>% dplyr::group_by(Genes) %>%
    dplyr::mutate(mean_zscore = mean(z_score)) %>% ungroup() 
  NT_df_tran <- NT_df_zscores  %>% filter(Genes%in%prots_perc$Genes) %>% dplyr::group_by(run_chan) %>% 
    dplyr::mutate("Trans_score" = mean(Zweighted)) %>% ungroup() %>% 
    distinct(run_chan,.keep_all=T) %>% dplyr::select("run_chan","Trans_score")
  NT_df_zscores <- NT_df_zscores %>% left_join(NT_df_tran, by =c("run_chan"= "run_chan"))
  
  NT_T_all <- rbind(NT_df_zscores, df_zscores)
  NT_T_all_uni <- NT_T_all %>% distinct(run_chan,.keep_all=T)
  NT_T_all_uni$Celltype <- factor(NT_T_all_uni$Celltype, levels=c("NUC_NT","NUC_10","NUC_30","NUC_60"))
  my_fin_colors <- c("violetred","skyblue1","royalblue1", "midnightblue")
  
  trans_plot <- ggplot(NT_T_all_uni, aes(x=Celltype, y=Trans_score, fill=Celltype)) + 
    geom_boxplot(alpha=0.3,outlier.shape = NA) + geom_jitter(shape=21,alpha=0.3,width=0.15) + 
    labs(x="LPS condition", y= "Transport score") + theme_classic() + 
    theme(legend.position = "none") + 
    scale_fill_manual(values = my_fin_colors)
  
  ggsave("Figures/TransportScore.pdf",trans_plot,height=3.4,width=2.6)
  
  dat_fin <- list(NT_T_all = NT_T_all, list_zscores = list_zscores, trans_plot=trans_plot)
  return(dat_fin)  
  
}

pD_SDs <- function(nucs_BC_MLFQ){
  dat <- nucs_BC_MLFQ$unimputed.BC
  m <- nucs_BC_MLFQ$no_neg
  results <- dat %>% left_join(m, by =c("run_chan" = "run_chan")) %>% 
    group_by(Genes,Celltype) %>% 
    dplyr::mutate("sd_withinTime" = sd(norm_prot,na.rm=T)) %>% ungroup() %>%
    group_by(Genes) %>% 
    dplyr::mutate("sd_acrossTime" = sd(norm_prot,na.rm=T)) %>% ungroup() %>%
    distinct(Genes,Celltype,.keep_all=T) %>% 
    dplyr::mutate("within_across_rat" = sd_withinTime/sd_acrossTime) %>% 
    group_by(Genes) %>% dplyr::mutate("MEAN_within_across_rat" = mean(within_across_rat,na.rm=T)) %>% ungroup() %>%
    dplyr::select(Genes,Celltype,sd_withinTime,sd_acrossTime,within_across_rat,MEAN_within_across_rat)
    
}

pD_TransportScore_cor <- function(dat, keep=c("NUC_60"),keepGenes=c("NUP62","PLEC","RAN")){

  df_zscores <- dat$NT_T_all
  df_zscores <- df_zscores %>% dplyr::group_by(Celltype) %>%
    dplyr::mutate("Trans_Score_norm" = Trans_score-median(Trans_score,na.rm=T)) %>% ungroup()
  df_zscores_LPS <- df_zscores[df_zscores$Celltype%in%keep,]

  cors <- df_zscores_LPS %>%
    group_by(Genes) %>% add_count() %>%
    dplyr::mutate(sp_r = cor(z_score, Trans_Score_norm, method="spearman")) %>% 
    dplyr::mutate(pr_r = cor(z_score, Trans_Score_norm, method="pearson")) %>% ungroup() %>% 
    distinct(Genes,.keep_all=T) %>% dplyr::select("Genes","mean_zscore","n","sp_r","pr_r")
  cors_timeSpecific <- df_zscores_LPS %>%
    group_by(Genes,Celltype) %>% add_count() %>%
    dplyr::mutate(sp_r = cor(z_score, Trans_Score_norm, method="spearman")) %>% 
    dplyr::mutate(pr_r = cor(z_score, Trans_Score_norm, method="pearson")) %>% ungroup() %>% 
    distinct(Genes,Celltype,.keep_all=T) %>% dplyr::select("Genes","mean_zscore","n","sp_r","pr_r","Celltype")
  cors_timeSpecific_BR <- df_zscores_LPS %>%
    group_by(Genes,Celltype,BioRep) %>% add_count() %>%
    dplyr::mutate(sp_r = cor(z_score, Trans_Score_norm, method="spearman")) %>% 
    dplyr::mutate(pr_r = cor(z_score, Trans_Score_norm, method="pearson")) %>% ungroup() %>% 
    distinct(Genes,Celltype,BioRep,.keep_all=T) %>% dplyr::select("Genes","mean_zscore","n","sp_r","pr_r","Celltype","BioRep")
  cors_br <- df_zscores_LPS %>%
    group_by(Genes,BioRep) %>% add_count() %>%
    dplyr::mutate(sp_r = cor(z_score, Trans_Score_norm, method="spearman")) %>% 
    dplyr::mutate(pr_r = cor(z_score, Trans_Score_norm, method="pearson")) %>% ungroup() %>% 
    distinct(Genes,BioRep,.keep_all=T) %>% dplyr::select("Genes","BioRep","mean_zscore","n","sp_r","pr_r")
  df_zscores_NT <- df_zscores[df_zscores$Celltype=="NUC_NT",]
  cors_NT <- df_zscores_NT %>%
    group_by(Genes) %>% add_count() %>%
    dplyr::mutate(sp_r = cor(z_score, Trans_Score_norm, method="spearman")) %>% 
    dplyr::mutate(pr_r = cor(z_score, Trans_Score_norm, method="pearson")) %>% ungroup() %>% 
    distinct(Genes,.keep_all=T) %>% dplyr::select("Genes","mean_zscore","n","sp_r","pr_r")
  cors_br_NT <- df_zscores_NT %>%
    group_by(Genes,BioRep) %>% add_count() %>%
    dplyr::mutate(sp_r = cor(z_score, Trans_Score_norm, method="spearman")) %>% 
    dplyr::mutate(pr_r = cor(z_score, Trans_Score_norm, method="pearson")) %>% ungroup() %>% 
    distinct(Genes,BioRep,.keep_all=T) %>% dplyr::select("Genes","BioRep","mean_zscore","n","sp_r","pr_r")

  top <- as.numeric(quantile(cors$pr_r, 0.99))
  bot <- as.numeric(quantile(cors$pr_r, 0.01))
  cors <- cors %>% dplyr::mutate(keep = ifelse(pr_r>top, "TRUE",
                                               ifelse(pr_r<bot,"TRUE","FALSE")))
  cors$keep[cors$pr_r>top,] <- T
  cors$keep[cors$pr_r<bot,] <- T
  
  cors <- cors %>% arrange(desc(pr_r))
  cors$Rank <- as.numeric(row.names(cors))
  
  plot_cors_all <- ggplot(cors, aes(x=Rank,y=pr_r)) + geom_point(size=0.7,alpha=0.4,color="black") +
    geom_point(data = cors %>% filter(Genes%in%keepGenes), size=0.7,alpha=0.4,color="red") +
    #geom_point(aes(y=mean_zscore),color="red") + 
    geom_hline(yintercept=0, linetype="dashed") + theme_classic() + 
    theme(axis.title.y=element_text(size=18)) +
    geom_label_repel(data=cors %>% filter(Genes%in%keepGenes),
                     aes(label=Genes),box.padding = 0.5) + 
    labs(y="Pearson Correlation")
  
  plot_cor_rank <- ggplot(cors, aes(x=Rank,y=pr_r)) + geom_point(size=0.7,alpha=0.4,color="black") +
    geom_point(data = cors %>% filter(Genes%in%keepGenes), size=1,alpha=0.8,color="red")+
    geom_hline(yintercept=0, linetype="dashed") +
    labs(x="Rank", y="Pearson Correlation")+
    geom_label_repel(data=cors %>% filter(Genes%in%keepGenes),
                     aes(label=Genes),box.padding = 0.5,min.segment.length = 0) + theme_bw() + 
    theme(axis.title.y=element_text(size=14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"))+
    geom_hline(yintercept=0, linetype="dashed")
  ggsave("Figures/corTS_ranks.pdf",plot_cor_rank,height=3.8,width=2.4)
  
  for(i in 1:length(keepGenes)){
    temp <- df_zscores_LPS[df_zscores_LPS$Genes%in%paste0(keepGenes[i]),]
    fit <- tls_fit(temp$z_score, temp$Trans_Score_norm)
    TLS_cor_prot <- ggplot(temp, aes(x=z_score, y=Trans_Score_norm)) + 
      geom_pointdensity() +
      scale_color_viridis() + 
      theme_classic() +
      stat_cor(method = "pearson") +
      theme(legend.position = "none") +
      geom_abline(intercept = fit$intercept, slope = fit$slope, color="red", linewidth=1,linetype="dashed") + #TLS fit
      #geom_smooth(method="lm")
      labs(x="NPC Z-score", y="Transport Score", subtitle=paste0(nrow(temp)," nuclei"))
    ggsave(filename=paste0("Figures/",keepGenes[i],"_zscore_TS.pdf"),TLS_cor_prot,height=3.4,width=3.2)
  }
  
  
  dat_fin <- list(plot_cors_all = plot_cors_all, cors_timeSpecific_BR=cors_timeSpecific_BR, cors_br_NT=cors_br_NT, cors_NT=cors_NT, cors=cors, cors_timeSpecific=cors_timeSpecific, cors_br=cors_br)
  return(dat_fin)  
}

pD_TS_cor_plot <- function(dat, NPC_regions, for_NPC_annot, scaffold, genes = c("O43242","Q9Y5Q3","Q9BRJ6","Q93009","Q9NRG0")){
  dat1 <- dat$cors
  dat1 <- dat1 %>% dplyr::mutate("keep" = ifelse(Genes%in%genes, T,F))
  NPC_regs <- data.frame(fread(NPC_regions))
  NPC_scaff <- data.frame(fread(scaffold))
  
  dat1 <- dat1 %>% left_join(NPC_regs, by =c("Genes"="Gene")) %>% left_join(NPC_scaff, by =c("Genes"="gene_name"))
  dat1 <- dat1 %>% arrange(desc(pr_r))
  dat1$Rank <- as.numeric(row.names(dat1))
  plot_cor_rank <- ggplot(dat1, aes(x=Rank,y=pr_r)) + 
    geom_point(data=dat1 %>% filter(keep==F), alpha=0.4,size=0.7,color="black") +
    geom_point(data=dat1 %>% filter(keep==T), size=1.5, alpha=0.8, color="red") +
    geom_hline(yintercept=0, linetype="dashed") +
    labs(x="Rank", y="Spearman Correlation")+
    geom_label_repel(data = dat1 %>% filter(keep==T) %>% filter(Genes%in%for_NPC_annot),
                     aes(label=Genes),box.padding = 0.5) + theme_bw() + 
    theme(axis.title.y=element_text(size=14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"))+
    geom_hline(yintercept=0, linetype="dashed")
  ggsave("Figures/Nucleoporins_corTS.pdf",plot_cor_rank,height=3.4,width=2.5)
  
  dat2 <- dat1
  dat2 <- dat2 %>% filter(!is.na(Localization))
  dat2_summary <- dat2 %>%
    group_by(Localization) %>%
    summarize(
      Mean = mean(pr_r),
      StdDev = sd(pr_r),
      SE = StdDev / sqrt(n()) 
    )
  
  dat2 <- dat2 %>% dplyr::mutate("type_local"=ifelse(grepl("Transmem|Inner|Central",Localization),"Inner","Outer"))
  t_test_result <- t.test(pr_r ~ type_local, data = dat2)
  dat2$Localization <- as.character(dat2$Localization)
  dat2$Localization <- factor(dat2$Localization, levels = c(
    "Nuclear basket", "Cytoplasmic region", "Y-Subcomplex", "Transmembrane", "Central channel","Inner ring"
  ))
  
  NPC_local <- ggplot(dat2) + 
    geom_jitter(aes(x=Localization, y=pr_r, fill=Localization), shape = 21,
                position=position_jitter(width=0.05),size=4,alpha=0.8) +
    geom_errorbar(data = dat2_summary,
                  aes(ymin=Mean-SE, ymax=Mean+SE, x=Localization), width=.2,
                  position=position_dodge(0.05), alpha=0.95) +
    theme_classic() + scale_fill_manual(values=c("grey10","grey25","grey40","grey65","grey80","grey95"))+
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + 
    labs(y="Correlation to transport score",subtitle=paste0("pval=",round(t_test_result$p.value,4)))
  
  ggsave("Figures/Nucleoporins_localization.pdf",NPC_local,height=3.2,width=2.7)
  
  ########################
  # scaffold or peripheral
  
  dat2 <- dat1
  dat2 <- dat2 %>% filter(!is.na(classification)) %>% filter(!grepl("undet",classification))
  dat2_summary <- dat2 %>% 
    group_by(classification) %>%
    summarize(
      Mean = mean(pr_r),
      StdDev = sd(pr_r),
      SE = StdDev / sqrt(n()) 
    )
  
  t_test_result <- t.test(pr_r ~ classification, data = dat2)
  dat2$classification <- as.character(dat2$classification)
  dat2$classification <- factor(dat2$classification, levels = c(
    "peripheral","scaffold"
  ))
  
  NPC_scaffold <- ggplot(dat2) + 
    geom_jitter(aes(x=classification, y=pr_r), shape = 21, fill="grey50",
                position=position_jitter(width=0.15),size=4,alpha=0.8) +
   # geom_beeswarm(aes(x=classification, y=pr_r), shape = 21, fill="grey50",
   #             size=4,alpha=0.8) +
    geom_errorbar(data = dat2_summary,
                  aes(ymin=Mean-SE, ymax=Mean+SE, x=classification), width=.2,
                  position=position_dodge(0.05), alpha=0.95) +
    theme_classic() + 
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + 
    labs(y="Correlation to transport score",subtitle=paste0("pval=",round(t_test_result$p.value,4)))
  ggsave("Figures/Nucleoporins_scaffold.pdf",NPC_scaffold,height=2.9,width=1.8)
  
  fin <- list(plot_cor_rank=plot_cor_rank, NPC_local=NPC_local)
  return(fin)
}

pD_NPC_Zscore_cor <- function(Trans_Score_weight_0.75_MLFQ, nucleoporins){
  dat <- Trans_Score_weight_0.75_MLFQ$NT_T_all
  dat <- dat[!grepl("NT",dat$Celltype),]

  dat <- dat %>% dplyr::group_by(Celltype) %>%
    dplyr::mutate("Trans_Score_norm" = Trans_score-median(Trans_score,na.rm=T)) %>% ungroup()
  NPCs <- dat[dat$Genes%in%nucleoporins,]
  NPCs_sum <- NPCs %>% group_by(run_chan) %>% add_count() %>%
    dplyr::mutate("NPC_ab" = sum(z_score)/sqrt(n)) %>% ungroup() %>% dplyr::distinct(run_chan,.keep_all=T)
  fit <- tls_fit(NPCs_sum$NPC_ab, NPCs_sum$Trans_Score_norm)
  
  NPC_cor_Z_plot_TLS <- ggplot(NPCs_sum, aes(x=NPC_ab, y=Trans_Score_norm)) + 
    geom_pointdensity() +
    scale_color_viridis() + 
    theme_classic() +
    stat_cor(method = "pearson") +
    theme(legend.position = "none") +
    geom_abline(intercept = fit$intercept, slope = fit$slope, color="red", linewidth=1,linetype="dashed") + 
    labs(x="NPC Z-score", y="Transport Score", subtitle=paste0(nrow(NPCs_sum)," nuclei"))
  ggsave("Figures/NPC_zscore_TS.pdf",NPC_cor_Z_plot_TLS,height=3.4,width=3.2)
  
  NPC_cor_Z_plot_allTimepoints <- ggplot(NPCs_sum, aes(x=NPC_ab, y=Trans_score)) + 
    facet_wrap(~Celltype)+
    geom_pointdensity() +
    scale_color_viridis() + 
    theme_classic() +
    stat_cor(method = "pearson") +
    theme(legend.position = "none") +
    labs(x="NPC Z-score", y="Transport Score", subtitle=paste0(nrow(NPCs_sum)," nuclei"))
  ggsave("Figures/Supp/NPC_zscore_TS_alltimepoints.pdf",NPC_cor_Z_plot_allTimepoints,height=3.4,width=8)
  
  # does accounting for the abundance of other proteins better explain the variance in NPC abundance and transport?
  
  NPCs_sum_vals <- dat %>% left_join(NPCs_sum %>% select(run_chan,NPC_ab), by=c("run_chan"="run_chan")) %>%
    filter(Genes%!in%nucleoporins)
  results <- NPCs_sum_vals %>%
    group_by(Genes) %>% 
    do({
      model <- lm(Trans_Score_norm ~ NPC_ab + z_score, data = .)
      predicted_values <- predict(model, .)
      pearson_cor <- cor(.$Trans_Score_norm, predicted_values, method = "pearson")
      spearman_cor <- cor(.$Trans_Score_norm, predicted_values, method = "spearman")
      
      data.frame(
        Genes = unique(.$Genes),
        Intercept = coef(model)[1],
        NPC_ab_coef = coef(model)[2],
        Z_score_coef = coef(model)[3],
        R_squared = summary(model)$r.squared,
        Pearson_Correlation = pearson_cor,
        Spearman_Correlation = spearman_cor
      )
    })
  
  most_cor <- results %>% ungroup() %>% arrange(desc(R_squared)) %>% slice_max(order_by = R_squared, n=1)
  NPCs_sum_vals_mostCor <- NPCs_sum_vals %>% ungroup() %>% filter(Genes%in%most_cor$Genes)
  fit <- tls_fit(NPCs_sum_vals_mostCor$NPC_ab, NPCs_sum_vals_mostCor$Trans_Score_norm)

  sunset <- c("#394d71","#5f97ce","white","#f89b62","#f05056")
  NPCs_sum_vals_mostCor$z_score[NPCs_sum_vals_mostCor$z_score>3] <- 3
  NPCs_sum_vals_mostCor$z_score[NPCs_sum_vals_mostCor$z_score< -3] <- -3
  
  NPC_cor_Z_plot_TLS_mostCor <- ggplot(NPCs_sum_vals_mostCor, aes(x=NPC_ab, y=Trans_Score_norm, fill=z_score)) + 
    scale_fill_gradientn(colors=sunset,limits=c(-3,3))+
    geom_point(shape=21,color="black",size=3) +
    theme_classic() +
    stat_cor(method = "pearson") + labs(fill=paste0(most_cor$Genes[1]))+
    geom_abline(intercept = fit$intercept, slope = fit$slope, color="red", linewidth=1,linetype="dashed") + 
    labs(x="NPC Z-score", y="Transport Score", subtitle=paste0(nrow(NPCs_sum_vals_mostCor)," nuclei"))
  return(list(NPCs_sum=NPCs_sum, NPC_cor_Z_plot_TLS=NPC_cor_Z_plot_TLS))
}

pD_halfLives <- function(cor_transScore_MLFQ, HL_path, annot = c("RANBP2", "NUP153", "NUP62", "NUP214", "NUP93", "NUP98", "NUP88", "NUP50")){
  HL <- data.frame(fread(HL_path))
  row.names(HL) <- HL$gene_name
  HL <- HL[,grepl("half_life|gene_name",colnames(HL))]
  HL_m <- melt(HL)
  HL_m <- HL_m[!grepl("R_sq",HL_m$variable),]
  HL_m$variable <- as.character(HL_m$variable)
  HL_m$cell_rep <- substr(HL_m$variable,1,nchar(HL_m$variable)-4)
  HL_m$cell <- substr(HL_m$variable,1,nchar(HL_m$variable)-6)
  HL_m$cell <- sub("\\..*", "", HL_m$variable)
  HL_m_result <- HL_m %>% group_by(cell,gene_name) %>% summarize(mean_HL = mean(value,na.rm=T)) %>% ungroup()
  HL_m_result_fin <- HL_m_result %>% group_by(gene_name) %>% summarize(mean_HL_all = mean(mean_HL,na.rm=T)) %>% ungroup()
  cors <- cor_transScore_MLFQ$cors
  cors <- cors %>% dplyr::select("Genes","pr_r")
  HL1 <- HL_m_result_fin %>% left_join(cors, by =c("gene_name"="Genes"))
  HL1$r <- as.numeric(HL1$pr_r)
  HL1 <- HL1[!is.na(HL1$r),]
  HL1 <- HL1[!grepl("NaN",HL1$mean_HL_all),]

  HL1_NUP <- HL1[HL1$gene_name%in%nucleoporins,]
  pr <- cor(HL1_NUP$mean_HL_all,HL1_NUP$r, method="pearson")
  sp <- cor(HL1_NUP$mean_HL_all,HL1_NUP$r, method="spearman")

  fit <- tls_fit(HL1_NUP$mean_HL_all, HL1_NUP$r)
  NPC_regs <- data.frame(fread(NPC_regions))
  HL1_NUP <- HL1_NUP %>% left_join(NPC_regs, by =c("gene_name"="Gene"))

  p_black <- ggplot(HL1_NUP, aes(x=(mean_HL_all), y=r)) + 
    theme_classic() + stat_cor(method = "pearson",label.x = 50, label.y = 0.24) +
    labs(x="Half life (hrs)", y ="Correlation to transport score") +
    geom_point(shape=21,fill="grey30",color="black",alpha=0.8, size=3.5)+
    geom_abline(intercept = fit$intercept, slope = fit$slope, color="red", linewidth=1,linetype="dashed") + #TLS fit
    geom_text_repel(aes(label=gene_name),min.segment.length = 0,color="black",size=2.75)
  
  ggsave("Figures/Nucleoporins_halfLife.pdf",p_black,height=3.2,width=3.2)
  
  return(p_black)
}

pD_passive_NPC <- function(Trans_Score_weight_0.75_MLFQ, nucleoporins, SN = c("gJD2507_RL1_1_1795mTRAQ4","gJD3177_B-L10_1_3386mTRAQ4")){
  
  dat <- Trans_Score_weight_0.75_MLFQ$NT_T_all
  SN_dif_ab <- dat %>% filter(Celltype!="NUC_NT") %>% group_by(Genes,Celltype) %>% add_count() %>%
    dplyr::mutate("combined_Zscore" = sum(z_score)/sqrt(n)) %>% 
    dplyr::mutate(mean_Zscore = mean(z_score)) %>% ungroup() %>% dplyr::distinct(Genes,Celltype,.keep_all=T) %>%
    dplyr::mutate(p.val = Z.to.p(combined_Zscore)) %>% dplyr::mutate(p.adj = p.adjust(p.val, method = "BH")) %>%
    dplyr::select(c("Genes","Celltype","mean_Zscore","p.adj")) %>% rename("Gene"="Genes")
  dat <- dat[!grepl("NT",dat$Celltype),] %>% rename("Gene"="Genes")
  NPCs_sum <- dat %>% 
    filter(Gene%in%nucleoporins) %>% group_by(run_chan) %>% add_count() %>%
    dplyr::mutate("NPC_ab" = sum(z_score)/sqrt(n)) %>% ungroup() %>% dplyr::distinct(run_chan,.keep_all=T) %>%
    dplyr::select("run_chan","NPC_ab") %>% arrange(NPC_ab) %>% dplyr::mutate("nuc_number"=row.names(.)) %>%
    dplyr::mutate("nuc_number"=as.numeric(nuc_number))
  
  keep_prots <- SN_dif_ab %>% filter(p.adj<0.05) %>% distinct(Gene,Celltype,.keep_all=T) %>% 
    dplyr::mutate("Gene_Celltype"=paste0(Gene,Celltype)) %>% dplyr::select("Gene_Celltype","p.adj")
  uniprot_masses <- read.delim(masses_fpath) %>% dplyr::mutate(log10_Mass=log10(Mass))
  dat2 <- dat %>% dplyr::mutate("Gene_Celltype"=paste0(Gene,Celltype)) %>% 
    inner_join(keep_prots, by =c("Gene_Celltype"="Gene_Celltype")) %>%
    inner_join(uniprot_masses, by =c("Gene" = "Gene.Names..primary."),relationship = "many-to-many") %>%
    filter(Gene%!in%nucleoporins) %>% #remove nucleoporins that change signficantly, as this is the factor of the x axis, and could be a circular association.
    inner_join(NPCs_sum, by =c("run_chan"="run_chan"))%>%
    dplyr::select("Gene","z_score","log10_Mass","run_chan","Celltype","NPC_ab","nuc_number") 
  dat2$z_score <- abs(dat2$z_score)
  dat2$z_score <- as.numeric(dat2$z_score)
  slopes <- dat2 %>%
    group_by(run_chan) %>%
    do(fit_slope_OLS(.)) %>%
    ungroup() 
  
  slopes_with_info <- slopes %>%
    left_join(dat2 %>%
                group_by(run_chan) %>%
                summarize(NPC_ab = first(NPC_ab),
                          Celltype = first(Celltype),
                          .groups = 'drop'),
              by = "run_chan")
  
  slopes_with_info_cors <- slopes_with_info %>% dplyr::group_by(Celltype) %>%
    dplyr::summarize(cor_slope = cor(slope,NPC_ab,method="pearson"),
                     cor_yint = cor(yint,NPC_ab,method="pearson")) %>% ungroup()
  sp_cor <- cor.test(slopes_with_info$slope,slopes_with_info$NPC_ab,method="spearman")
  pr_cor <- cor.test(slopes_with_info$slope,slopes_with_info$NPC_ab,method="pearson")
  tls <- tls_fit(slopes_with_info$NPC_ab,slopes_with_info$slope)
  NPC_andMassSlope <- ggplot(slopes_with_info, aes(x=NPC_ab,y=slope)) +
    geom_pointdensity() +
    scale_color_viridis() + 
    theme_classic() +
    stat_cor(method = "pearson") +
    theme(legend.position = "none") +
   # geom_smooth(color="red",se=F) +
    geom_hline(yintercept=0, linetype="dashed")+
    geom_point(data=slopes_with_info %>% filter(run_chan%in%SN), aes(x=NPC_ab,y=slope),shape=21,color="black",fill="red2",size=3)+
    labs(x="NPC Z-score", y="Slope (protein ab. vs mass)") +
    geom_abline(intercept = tls$intercept, slope = tls$slope, color="red", linewidth=1,linetype="solid") #TLS fit
  ggsave("Figures/NPC_zscore_andSlope.pdf",NPC_andMassSlope,height=3.4,width=3.2)
  
  plot_select <- ggplot(dat2 %>% filter(run_chan%in%SN), aes(x=log10_Mass,y=z_score,group=run_chan)) +
    geom_point(shape=21,fill="grey40",color="black",alpha=0.8) + facet_wrap(~run_chan,nrow=1)+
    theme_bw() +
    theme(legend.position = "none",            
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black")) +
    geom_smooth(method="lm", color="red", se=F) +
    labs(x="Log10, Protein mass", y="|Protein z-score|") +
    geom_text(data=slopes_with_info %>% filter(run_chan%in%SN), aes(label=paste("Slope:", round(slope, 2))), x=4.3, y=5.3, inherit.aes=FALSE)
  
  ggsave("Figures/NPC_zscoreSelect.pdf",plot_select,height=2.2,width=3.2)
  return(slopes_with_info)
}

pD_SN_difAb <- function(Trans_Score_weight_0.75_MLFQ){
  dat <- Trans_Score_weight_0.75_MLFQ$NT_T_all
  SN_dif_ab <- dat %>% filter(Celltype!="NUC_NT") %>% group_by(Genes,Celltype) %>% add_count() %>%
    dplyr::mutate("combined_Zscore" = sum(z_score)/sqrt(n)) %>% 
    dplyr::mutate(mean_Zscore = mean(z_score)) %>% ungroup() %>% dplyr::distinct(Genes,Celltype,.keep_all=T) %>%
    dplyr::mutate(p.val = Z.to.p(combined_Zscore)) %>% dplyr::mutate(p.adj = p.adjust(p.val, method = "BH")) %>%
    dplyr::select(c("Genes","Celltype","mean_Zscore","p.adj")) %>% rename("Gene"="Genes")
  
  SN_dif_ab_all <- dat %>% filter(Celltype!="NUC_NT") %>% group_by(Genes) %>% add_count() %>%
    dplyr::mutate("combined_Zscore" = sum(z_score)/sqrt(n)) %>% 
    dplyr::mutate(mean_Zscore = mean(z_score)) %>% ungroup() %>% dplyr::distinct(Genes,.keep_all=T) %>%
    dplyr::mutate(p.val = Z.to.p(combined_Zscore)) %>% dplyr::mutate(p.adj = p.adjust(p.val, method = "BH")) %>%
    dplyr::select(c("Genes","mean_Zscore","p.adj")) %>% rename("Gene"="Genes")
  
  return(SN_dif_ab_all)
}

pD_SN_cors_difAb <- function(Trans_Score_weight_0.75_MLFQ, Transport_Score_cor){
  SN_DA <- pD_SN_difAb(Trans_Score_weight_0.75_MLFQ)
  cors <- Transport_Score_cor$cors %>% left_join(SN_DA, by =c("Genes"="Gene")) 
}

fit_slope_OLS <- function(df) {
  model <- lm(z_score ~ log10_Mass, data = df)
  slope <- coef(model)[2] #slope coefficient
  yint <- coef(model)[1] #yint coefficient
  return(data.frame(slope = slope,yint=yint))
}

join_datasets <- function(df_list){
  combined_df <- do.call("rbind", df_list)
  return(combined_df)
}

rm_data <- function(dat, rm=c("FileName")){
  dat <- dat[dat$run_chan%!in%rm,]
  return(dat)
}

pD_impute_bulk <- function(dat, missing.prot.frac=0.7, missing.cell.frac=0.95, k=10, Label=T, BioRep=T, LC_batch=T){
  df <- dat
  if(Label==F){
    df$Label <- 0
  }
  if(BioRep==F){
    df$BioRep <- 0
  }
  if(LC_batch==F){
    df$LC_batch <- 0
  }
  df$batch <- paste0(df$Label,"_",df$BioRep,"_",df$LC_batch)
  df$log2val <- log2(df$value)
  df <- df %>% dplyr::group_by(run_chan) %>% dplyr::mutate(norm = log2val-median(log2val, na.rm=T)) %>% 
      ungroup() %>% dplyr::group_by(Gene, batch) %>% dplyr::mutate(norm = norm-mean(norm, na.rm=T)) %>% ungroup()
  dat_fin1 <- df %>% dplyr::select("Gene", "run_chan", "norm")
  dat_fin1$Gene <- as.factor(dat_fin1$Gene)
  dat_fin1$run_chan <- as.factor(dat_fin1$run_chan)
  dat_fin1$norm <- as.numeric(dat_fin1$norm)
  dat_fin2<-dcast(dat_fin1, Gene~run_chan, value.var = "norm", fill=NA)
  dat3<-as.matrix(dat_fin2[,-1]); row.names(dat3)<-dat_fin2[,1]
  
  ## Impute
  pre_impute <-(filt.mat.rc(dat3,missing.prot.frac,missing.cell.frac))
  sc.imp<-cr_norm_log(pre_impute)
  sc.imp <- 2^sc.imp
  sc.imp <- hknn_weighted_EucNorm(sc.imp, k)
  sc.imp <- log2(sc.imp)
  sc.imp<-cr_norm_log((sc.imp))
  
  dat <- list(pre_impute = pre_impute, sc.imp = sc.imp)
  
  return(dat)
}

pD_bulk_BC <- function(dat, meta_bulk_path, Label=T, BioRep=T, LC_batch=T, Cell_Batch=T,use_mod=T){
  sc.imp <- dat$sc.imp
  m <- fread(meta_bulk_path)
  if(Label==F){
    m$Label <- 0
  }
  if(BioRep==F){
    m$BioRep <- 0
  }
  if(Cell_Batch==F){
    m$Cell_Batch <- 0
  }
  if(LC_batch==F){
    m$LC_batch <- 0
  }
  m$batch <- paste0(m$Label,"_",m$BioRep,"_",m$LC_batch,"_",m$Cell_Batch)
  m$run_chan <- paste0(m$Raw,"mTRAQ",m$Label)
  colnames(sc.imp) <- gsub(".", "-", colnames(sc.imp), fixed = TRUE) #in case column names changed hyphen to dot
  batch.covs <-m$batch[match(colnames(sc.imp), m$run_chan)]
  mod<-data.frame(m$Celltype[match(colnames(sc.imp), m$run_chan)]); colnames(mod)<-"celltype"
  mod<-model.matrix(~as.factor(celltype), data=mod)
  if(use_mod==F){
    imputed.BC <- ComBat(sc.imp, batch=batch.covs)
  } else{
    imputed.BC <- ComBat(sc.imp, batch=batch.covs, mod=mod)
  }
  imputed.BC <- cr_norm_log(imputed.BC)
  
  mat.sc.imp_NA <- imputed.BC
  unimp <- dat$pre_impute
  mat.sc.imp_NA[is.na(unimp)==T] <- NA #replace imputed with NA
  mat.sc.imp_NA <- cr_norm_log(mat.sc.imp_NA)
  
  unimputed.BC <- melt(mat.sc.imp_NA); colnames(unimputed.BC)<-c("Gene","run_chan","norm_prot")
  no_neg <- fread(meta_bulk_path)
  no_neg$run_chan <- paste0(no_neg$Raw,"mTRAQ",no_neg$Label)
  m <- fread(meta_bulk_path)
  m$run_chan <- paste0(m$Raw,"mTRAQ",m$Label)
  
  dat_fin <- list(unimputed.BC = unimputed.BC, imputed.BC = imputed.BC,m=m, no_neg=no_neg)
  return(dat_fin)
}

pD_effect_on_transport <- function(dat_BC_bulk_siRNA1, meta_bulk_path_timsTOFSCP_siRNA, bulk_DA_all_LMM, cor_transScore_MSEmp_BR1236_sum_03012024, titlename="top100",mostdif=50, useall=T){
  dat <- dat_BC_bulk_siRNA1$unimputed.BC
  m <- fread(meta_bulk_path_timsTOFSCP_siRNA)
  m$run_chan <- paste0(m$Raw,"mTRAQ",m$Label)
  dat <- dat %>% inner_join(m, by =c("run_chan"="run_chan"))

  if(useall==T){
    DA <- bulk_DA_all_LMM %>% filter(!grepl("KRT",Gene)) %>%
      filter(Cond == "min_60") %>% filter(p.adj<0.01)
    dat <- dat %>% filter(grepl("NUC",Type))
  } else{
    DA <- bulk_DA_all_LMM %>% filter(!grepl("KRT",Gene)) %>%
      filter(Cond == "min_60") %>% 
      arrange(p.val) %>% 
      slice(1:mostdif)
    dat <- dat %>% filter(grepl("NUC",Type)) %>% filter(Gene%in%DA$Gene)
  }
  
  DA$zscore <- p.to.Z(DA$p.adj) #convert pval to zscore
  DA$zscore[is.na(DA$zscore)] <- max(DA$zscore,na.rm=T)
  weights <- DA %>% dplyr::select("Gene","zscore")
  DA$log2FC <- DA$log2FC*DA$zscore #weight the FC by the z-score (more change for more significant prots)
  dat <- dat %>%
    group_by(Cell_Batch,Gene,Celltype,Cond) %>% 
    dplyr::mutate(mean_norm_prot=mean(norm_prot,na.rm=T)) %>% ungroup() %>% 
    distinct(Cell_Batch,Gene,Celltype,Cond,.keep_all=T)
  dat$gene_Cell_Batch <- paste0(dat$Gene,"_",dat$Cell_Batch,"_",dat$Celltype)
  dat_LPS <- dat[grepl("LPS",dat$Cond),] %>% dplyr::select("gene_Cell_Batch","Cell_Batch","Gene","Celltype","mean_norm_prot")
  dat_NT <- dat[grepl("NT",dat$Cond),] %>% dplyr::select("gene_Cell_Batch","mean_norm_prot") %>% 
    dplyr::rename("norm_prot_NT"="mean_norm_prot")
  dat_both <- dat_LPS %>% full_join(dat_NT, by =c("gene_Cell_Batch"="gene_Cell_Batch"))
  dat_both$log2FC <- dat_both$mean_norm_prot-dat_both$norm_prot_NT
  dat_both$Gene_Batch <- paste0(dat_both$Gene,"_",dat_both$Cell_Batch)
  DA_lim <- DA %>% dplyr::select("Gene","log2FC") %>% dplyr::rename("log2_FC_bulk" = "log2FC")
  dat_both <- dat_both %>% left_join(DA_lim, by =c("Gene"="Gene")) %>% dplyr::mutate("id"=paste0(Celltype,Cell_Batch))
  unicond <- dat_both %>% filter(Celltype!="Neg_kd") %>% distinct(Celltype)
  unicond <- as.vector(unicond$Celltype)
  temp_list <- list()

  for(i in 1:length(unicond)){
    temp <- dat_both[dat_both$Celltype==paste0(unicond[i]),] %>% na.omit()
    uni_batches <- unique(temp$Cell_Batch)
    temp_neg <- dat_both[dat_both$Celltype=="Neg_kd",] %>% filter(Cell_Batch%in%uni_batches)
    temp_dat <- rbind(temp, temp_neg)
    uni_batches <- unique(temp_dat$Cell_Batch)
    temp_list2 <- list()
    for(k in 1:length(uni_batches)){
      temp_dat2 <- temp_dat %>% filter(Cell_Batch==paste0(uni_batches[k]))
      temp_dat2 <- temp_dat2  %>%
        dplyr::rename("mean_log2FC"="log2FC") %>% 
        ungroup() %>% distinct(Gene,Celltype,Cell_Batch,.keep_all=T) %>% na.omit() %>% add_count(Gene) %>% filter(n==2)
      temp_dat2$Celltype <- factor(temp_dat2$Celltype, levels=c(paste0(unicond[i]),"Neg_kd"))
      temp_dat2$Comparison <- paste0(unicond[i],"_Neg_kd")
      temp_dat2 <- temp_dat2 %>% left_join(weights,by=c("Gene"="Gene")) %>% dplyr::mutate("mean_log2FC"=mean_log2FC*zscore) #weight by differential abundance
      sl <- temp_dat2 %>%
        group_by(Comparison, Celltype) %>%
        do(tidy_model = tidy(lm(mean_log2FC ~ log2_FC_bulk, data = .))) %>%
        unnest(tidy_model) %>%
        filter(term == "log2_FC_bulk") %>%
        select(-c(statistic, term)) %>%
        rename(slope = estimate) %>%
        ungroup()
      sl_d <- dcast(sl, Comparison~Celltype, value.var="slope")
      sl$Celltype <- paste0(sl$Celltype,"_SE")
      sl$Celltype <- factor(sl$Celltype, levels=c(paste0(unicond[i],"_SE"),"Neg_kd_SE"))
      
      sl_SE <- dcast(sl, Comparison~Celltype, value.var="std.error")
      both <- sl_d %>% left_join(sl_SE,by=c("Comparison"="Comparison"))
      z <- (both[1,2]-both[1,3])/sqrt(both[1,4]^2+both[1,5]^2)
      
      sl_d <- both %>% select(-"Neg_kd", everything(), "Neg_kd")
      colnames(both)[2] <- "KD"
      colnames(both)[4] <- "KD_SE"
      
      both$delta_slope <- both$KD-both$Neg_kd
      both <- both %>% dplyr::select("delta_slope","Comparison","KD","Neg_kd","KD_SE","Neg_kd_SE") %>% distinct(Comparison,.keep_all=T)
      both$z <- z
      temp_dat2$Celltype = substr(temp_dat2$Celltype,1,nchar(as.character(temp_dat2$Celltype))-3)
      temp_dat2 <- temp_dat2 %>% dplyr::select("Comparison","Celltype","Gene","Cell_Batch","log2_FC_bulk","mean_log2FC") %>% left_join(both, by=c("Comparison"="Comparison"))
      temp_list2[[k]] <- temp_dat2
    }
    all_reps <- do.call("rbind", temp_list2)
    
    temp_list[[i]] <- all_reps
  }
  all_KDs1 <- do.call("rbind", temp_list)
  all_KDs <- all_KDs1 %>% filter(!grepl("Neg",Celltype))%>%distinct(Celltype,Cell_Batch,.keep_all=T) %>% group_by(Comparison) %>% add_count() %>% dplyr::mutate(z_sum=sum(z)) %>%
    dplyr::mutate("mean_delta_slope"=mean(delta_slope)) %>%
    ungroup() %>% distinct(Comparison,.keep_all=T) %>% 
    dplyr::mutate(z_combined=z_sum/sqrt(n)) %>%
    dplyr::mutate("pval"=Z.to.p(z_combined)) %>% dplyr::mutate("pval.adj"=p.adjust(pval, method = "BH")) %>% 
    dplyr::select("Comparison","Celltype","mean_delta_slope","z_combined","pval","pval.adj") %>% arrange(mean_delta_slope)
  
  damps <- c("PLEC","SND1","IMMT","HSD17B4","RPSA","VIM","TPR")
  enhancers_KDs <- all_KDs1[apply(sapply(damps, function(str) !grepl(str, all_KDs1$Comparison)), 1, all), ]
  damps_KD <- all_KDs1[all_KDs1$Comparison%!in%enhancers_KDs$Comparison,]
  damps_KD$Celltype <- factor(damps_KD$Celltype, levels=c("PLEC","SND1","IMMT","HSD17B4","RPSA","VIM","TPR","Neg"))
  enhancers_KDs$Celltype <- factor(enhancers_KDs$Celltype, levels=c("NOP56","MYBBP1A","BCLAF1","DHX15","RRS1","BOP1","NUP205","RBM6","MDN1","Neg"))
  damps_KD$siRNA <- "Potential dampener"
  enhancers_KDs$siRNA <- "Potential enhancer"
  enhance_damp <- rbind(enhancers_KDs,damps_KD)
  enhance_damp <- enhance_damp %>%
    group_by(Comparison) %>%  # Group by 'Comparison'
    mutate(Cell_Batch_new = as.numeric(as.factor(as.character(Cell_Batch)))) %>% ungroup()
  enhance_damp$Comparison_new = sub("_.*", "", enhance_damp$Comparison)
  enhance_damp <- enhance_damp %>% dplyr::mutate("siRNA" = ifelse(grepl("Neg",Celltype), "Negative control", siRNA))
  enhance_damp$Comparison_new <- factor(enhance_damp$Comparison_new, levels = unique(all_KDs$Celltype))
  allInd_plot <- ggplot(enhance_damp, aes(x=log2_FC_bulk, y=mean_log2FC, color=siRNA, group=Celltype, fill=siRNA)) +
    geom_point(alpha=0.2)+geom_smooth(method="lm", formula=y~x, se=TRUE, linetype="solid")+
    geom_text(aes(x=-0.5, y=2, label=round(delta_slope,2)), inherit.aes=FALSE) +
    theme_classic() +
    facet_grid(Comparison_new~ Cell_Batch_new, scales="free_x",drop=TRUE) +
   # coord_cartesian(xlim=c(-1.5, 3), ylim=c(-2, 3))+
    scale_color_manual(values=c("#666666","#C567EE","#95D3C9")) +
    scale_fill_manual(values=c("#666666","#C567EE","#95D3C9"))
  ggsave(paste0("Figures/siRNA_allIndividual_deltaSlope_",titlename,".pdf"),allInd_plot,height=12.5,width=9)
  
  #####################
  enhancers_KDs <- all_KDs[apply(sapply(damps, function(str) !grepl(str, all_KDs$Comparison)), 1, all), ]
  damps_KD <- all_KDs[all_KDs$Comparison%!in%enhancers_KDs$Comparison,]
  damps_KD$Celltype <- factor(damps_KD$Celltype, levels=c("PLEC","SND1","IMMT","HSD17B4","RPSA","VIM","TPR","Neg"))
  enhancers_KDs$Celltype <- factor(enhancers_KDs$Celltype, levels=c("NOP56","MYBBP1A","BCLAF1","DHX15","RRS1","BOP1","NUP205","RBM6","MDN1","Neg"))
  damps_KD$Type <- "Potential Dampener"
  enhancers_KDs$Type <- "Potential Enhancer"
  
  enhance_damp <- rbind(enhancers_KDs,damps_KD)
  enhance_damp$Genes <- enhance_damp$Celltype
  
  pval_deltaSlope <- ggplot() + geom_point(data=enhance_damp %>% filter(pval.adj<0.05), aes(x=mean_delta_slope, y=-log10(pval.adj),fill=Type), shape=21, size=3.5,color="black") +
    geom_point(data=enhance_damp %>% filter(pval.adj>=0.05), aes(x=mean_delta_slope,y=-log10(pval.adj),fill=Type), shape=21, size=3.5, color="black",alpha=0.5) +
    scale_fill_manual(values=c("#7C65A8","#95D3C9")) + theme_classic() + geom_hline(yintercept=-log10(0.05),linetype="dashed")+
    geom_label_repel(data=enhance_damp%>% filter(pval.adj<0.05), 
                     aes(label=Genes, x=mean_delta_slope, y=-log10(pval.adj)),box.padding = 0.5,min.segment.length = 0,max.overlaps=Inf,size=2.5) +
    theme(legend.position = "none")
  ggsave(paste0("Figures/siRNA_pval_deltaSlope_",titlename,".pdf"),pval_deltaSlope,height=4,width=3.5)
  
  #######
  cors <- cor_transScore_MSEmp_BR1236_sum_03012024$cors %>% dplyr::mutate("Celltype" = paste0(Genes))
  all_KDs_lim <- all_KDs %>% left_join(cors, by =c("Celltype"="Celltype")) %>% distinct(Celltype,.keep_all=T) %>% 
    dplyr::select(c("Celltype","pr_r","mean_delta_slope")) %>% filter(!grepl("Neg",Celltype))
  fit <- tls_fit(all_KDs_lim$mean_delta_slope, all_KDs_lim$pr_r)
  all_KDs_lim$Genes = all_KDs_lim$Celltype
  
  # Create the plot
  delta_slope_cor <- ggplot() + 
    geom_point(data=all_KDs_lim %>% filter(Celltype%in%enhancers_KDs$Celltype), aes(x=mean_delta_slope, y=pr_r),color="black",fill="#95D3C9",shape=21,size=3.5,alpha=0.8) +
    geom_point(data=all_KDs_lim %>% filter(Celltype%in%damps_KD$Celltype), aes(x=mean_delta_slope, y=pr_r),color="black",fill="#4C467E",shape=21,size=3.5,alpha=0.8) +
    theme_classic() +
    theme(legend.position = "none") +
    geom_hline(yintercept = 0,linetype="dashed")+
    geom_vline(xintercept = 0,linetype="dashed")+
    geom_abline(intercept = fit$intercept, slope = fit$slope, color="red", linewidth=1,linetype="solid") + #TLS fit
    #geom_smooth(method="lm")
    labs(x="Mean delta slope from knock-down", y="Correlation to transport score")+
    geom_label_repel(data=all_KDs_lim, aes(label=Genes, x=mean_delta_slope, y=pr_r),box.padding = 0.5,min.segment.length = 0,max.overlaps=Inf,size=2.5) +
    stat_cor(data=all_KDs_lim, aes(x=mean_delta_slope, y=pr_r),method = "pearson") 
  ggsave(paste0("Figures/siRNA_corTS_deltaSlope_",titlename,".pdf"),delta_slope_cor,height=3.9,width=4.3)
  
  ########################## non parametric cors BR
  ############################
  cors_BR <- cor_transScore_MSEmp_BR1236_sum_03012024$cors_br %>%
    group_by(BioRep) %>%
    arrange(pr_r) %>%
    dplyr::mutate(
      BR_rank = row_number(), #rank within each BR
      BR_total = n(),
      Rankit = qnorm((BR_rank - 0.5) / BR_total)
    ) %>%
    ungroup() %>%
    group_by(Genes) %>% add_count(name="num_BRs") %>%
    dplyr::mutate(
      NP_Z = sum(Rankit)/sqrt(num_BRs)  #Stouffers method combining rankit values
    ) %>% ungroup() %>% distinct(Genes,.keep_all=T) %>%
    dplyr::select("Genes","NP_Z")
  cors_BR$pval <- Z.to.p(cors_BR$NP_Z)
  cors_BR <- cors_BR %>% dplyr::mutate(pval_adj = p.adjust(pval, method = "BH"))
  
  Zs <- cors_BR %>% inner_join(enhance_damp, by =c("Genes"="Genes"))
  fit <- tls_fit(Zs$z_combined, Zs$NP_Z)
  cor.test(Zs$NP_Z, Zs$z_combined, method="pearson")
  
  delta_slope_cor_Zscores <- ggplot() + 
    geom_point(data=Zs %>% filter(Celltype%in%enhancers_KDs$Celltype), aes(x=z_combined, y=NP_Z),color="black",fill="#95D3C9",shape=21,size=3.5,alpha=0.8) +
    geom_point(data=Zs %>% filter(Celltype%in%damps_KD$Celltype), aes(x=z_combined, y=NP_Z),color="black",fill="#4C467E",shape=21,size=3.5,alpha=0.8) +
    
    theme_classic() +
    theme(legend.position = "none") +
    geom_hline(yintercept = 0,linetype="dashed")+
    geom_vline(xintercept = 0,linetype="dashed")+
    geom_abline(intercept = fit$intercept, slope = fit$slope, color="red", linewidth=1,linetype="solid") + #TLS fit
    #geom_smooth(method="lm")s
    labs(x="Z-score, effect of KD on transport", y="Z-score, correlation of protein abundance to transport")+
    geom_label_repel(data=Zs, aes(label=Genes, x=z_combined, y=NP_Z),box.padding = 0.5,min.segment.length = 0,max.overlaps=Inf,size=2.5) +
    stat_cor(data=Zs, aes(x=z_combined, y=NP_Z),method = "pearson") 
  ggsave(paste0("Figures/siRNA_Zscores_corTS_deltaSlope_",titlename,".pdf"),delta_slope_cor_Zscores,height=3.9,width=4.3)
  
  fin <- list(delta_slope_cor_Zscores=delta_slope_cor_Zscores, delta_slope_cor=delta_slope_cor, pval_deltaSlope=pval_deltaSlope, allInd_plot=allInd_plot,Zs=Zs, cors_BR=cors_BR)
  return(fin)
}

pD_bulk_DA_all <- function(dat, collapse_conds=F){
  df <- dat$unimputed.BC
  m <- dat$m
  df <- df %>% left_join(m, by =c("run_chan" = "run_chan"))
  unirunchan <- unique(df$run_chan)
  if(collapse_conds==T){
    df <- df %>% group_by(Celltype,Label,Cell_Batch,Gene) %>%
      dplyr::mutate("norm_prot" = mean(norm_prot,na.rm=T)) %>%
      ungroup() %>% dplyr::mutate("run_chan"=paste0(Celltype,Label,Cell_Batch)) %>%
      distinct(run_chan,Gene,.keep_all=T)
  }
  list_pvals <- list()
  uni_conds <- unique(df$Celltype)
  uni_conds <- uni_conds[uni_conds!="Neg_kd"]
  for(i in 1:length(uni_conds)){
    df_temp <- df[df$Celltype=="Neg_kd"|df$Celltype==uni_conds[i],]
    df_temp$Celltype <- factor(df_temp$Celltype, levels=c(paste0(uni_conds[i]),"Neg_kd"))
    
    df_temp <- df_temp %>% group_by(run_chan) %>% dplyr::mutate("norm" = norm_prot-median(norm_prot,na.rm=T)) %>% ungroup() %>%
      group_by(Gene) %>% dplyr::mutate("norm" = norm-mean(norm,na.rm=T)) %>% ungroup() %>% na.omit() %>%
      add_count(Gene,Celltype) %>% filter(n>2) %>% 
      dplyr::select(-c("n"))
    keep_prots <- df_temp %>% distinct(Celltype,Gene) %>% add_count(Gene) %>% filter(n==2)
    df_temp <- df_temp[df_temp$Gene%in%keep_prots$Gene,]
    
    vals <- df_temp %>% dplyr::group_by(Gene, Celltype) %>% dplyr::mutate(mean_val=mean(norm,na.rm=T)) %>% 
      ungroup() %>% distinct(Gene, Celltype, .keep_all=T)
    vals <- dcast(vals, Gene~Celltype, value.var = "mean_val")
    vals$log2_FC <- vals[,2]-vals[,3]
    vals <- vals[,c(1,4)]
    
    df_pvals <-  df_temp %>% group_by(Gene) %>%
      t_test(norm ~ Celltype) %>%
      adjust_pvalue(method = "BH") %>%
      add_significance()
    df_pvals <- df_pvals %>% left_join(vals, by =c("Gene"="Gene"))
    list_pvals[[i]] <- df_pvals
    print(paste0(i,"/",length(uni_conds)," done"))
  }
  pvals <- do.call("rbind", list_pvals)
  return(pvals)
}

pD_plot_DA_v2 <- function(bulk_DA, xlab="Log2, 30 min LPS / NT", titleLab = "Nuclear fraction"){
  uni_conds <- unique(bulk_DA$group1)
  plot_list <- list()
  for(i in 1:length(uni_conds)){
    bulk_DA_lim <- bulk_DA %>% filter(group1==paste0(uni_conds[i])) %>% filter(group2=="Neg_kd")
    specify <- substr(uni_conds[i],1,nchar(uni_conds[i])-3)
    result <- bulk_DA_lim %>% 
      dplyr::mutate("direction" = ifelse(log2_FC>0 & p.adj<0.05, "Increased",
                                         ifelse(log2_FC<0 & p.adj<0.05, "Decreased", "Not significant")))
    
    volc_plot <- ggplot(result, aes(x=log2_FC, y = -log10(p.adj), color=direction,fill=direction)) + 
      geom_point(size=1, alpha=0.6) + 
      geom_point(data= result %>% filter(Gene%in%specify), shape=21,size=2, alpha=0.6,color="black") + 
      scale_color_manual(values=c("Decreased"="steelblue","Increased"="tomato1","Not significant"="gray30")) +
      scale_fill_manual(values=c("Decreased"="steelblue","Increased"="tomato1","Not significant"="gray30")) + 
      geom_label_repel(
        data = result %>% filter(Gene%in%specify),
        aes(label = Gene), color="black",
        size = 3,
        min.segment.length = unit(0, 'lines'),
        point.padding = unit(0.3, "lines"),
        segment.color = 'black') +
      labs(x=xlab, y="-Log10, q-value", color="", title = paste0(specify," knock-down")) +
      theme(legend.position="none")
    plot_list[[i]] <- volc_plot
  }
  
  fin_list <- list(plot_list=plot_list)
  return(fin_list)
}

pD_PSEA_bulk_siRNA <- function(bulk,keep_KDs,min_prot=3,min_fraction = 0.05, ontology_sources =c("GO")){

  bulk1 <- bulk %>% dplyr::select("log2_FC","Gene","group1") %>% 
    filter(group1%in%keep_KDs) %>% dplyr::arrange(desc(log2_FC)) 
  uni_cond <- unique(bulk1$group1)
  empty_list <- list()
  for(i in 1:length(uni_cond)){
    bulk1_temp <- bulk1[bulk1$group1==paste0(uni_cond[i]),]
    res <- gprofiler2::gost(bulk1_temp$Gene, organism="hsapiens", ordered_query=T, 
                            user_threshold=0.01, sources=ontology_sources, correction_method='fdr',evcodes=T)
    res2 <- res$result
    res2 <- res2[res2$intersection_size >= min_prot,] #require certain number of proteins per GO term
    res2 <- res2[(1-res2$precision) >= min_fraction,] #require certain fraction of GO term coverage
    res_meta <- res$meta
    res2$med_ab <- 0 #place holder
    res2$mean_ab <- 0 #place holder
    colnum <- ncol(res2)
    for(k in 1:nrow(res2)){
      res2_temp <- res2[k,]
      res2_temp <- unlist(strsplit(res2_temp$intersection,split=","))
      dat_temp_abs <- bulk1_temp[bulk1_temp$Gene%in%res2_temp,]
      dat_temp_abs_med <- median(dat_temp_abs$log2_FC) #get median value
      dat_temp_abs_mean <- mean(dat_temp_abs$log2_FC) #get mean value
      
      res2[k,colnum-1] <- dat_temp_abs_med
      res2[k,colnum] <- dat_temp_abs_mean
      res2$Cond <- paste0(uni_cond[i])
    }
    empty_list[[i]] <- res2
    print(paste0(i,"/",length(uni_cond)," finished")) 
  }
  df <- do.call("rbind", empty_list)
  df <- convert_list_columns(df)
  
  return(df)
}

pD_plot_PSEA_siRNA_bulk <- function(PSEA_KDs_WC, bulk, Trans_Score_siRNA1_slopes_MSEmp_ind, GOterms){

  df_uniqueGeneTerms <- PSEA_KDs_WC %>%
    separate_rows(intersection) %>% distinct(term_id, intersection)
  uni_conds <- unique(PSEA_KDs_WC$Cond)
  empty_list <- list()
  for(i in 1:length(uni_conds)){
    temp <- PSEA_KDs_WC %>% distinct(term_id, .keep_all=T)
    bulk_DA_temp <- bulk[bulk$group1==paste0(uni_conds[i]),]
    temp$med_ab <- 0 
    temp$mean_ab <- 0 
    colnum <- ncol(temp)
    for(k in 1:nrow(temp)){
      df_uniqueGeneTerms_temp <- df_uniqueGeneTerms[df_uniqueGeneTerms$term_id==temp$term_id[k],]
      dat_temp_abs <- bulk_DA_temp[bulk_DA_temp$Gene%in%df_uniqueGeneTerms_temp$intersection,]
      dat_temp_abs_med <- median(dat_temp_abs$log2_FC) 
      dat_temp_abs_mean <- mean(dat_temp_abs$log2_FC) 
      
      temp[k,colnum-2] <- dat_temp_abs_med
      temp[k,colnum-1] <- dat_temp_abs_mean
      temp$Cond <- paste0(uni_conds[i])
    }
    empty_list[[i]] <- temp
    print(paste0(i,"/",length(uni_conds)," finished")) 
    
  }
  df <- do.call("rbind", empty_list)
  df_min0 <- df %>% distinct(term_id, .keep_all=T)
  df_min0$Cond <- "Neg_kd"
  df_min0$med_ab <- 0
  df_min0$mean_ab <- 0
  all <- rbind(df, df_min0)
  all <- all %>% dplyr::group_by(term_id) %>% dplyr::mutate("norm_ab" = mean_ab-mean(mean_ab,na.rm=T)) %>%
    dplyr::mutate(sd=sd(norm_ab)) %>% ungroup() %>% filter(intersection_size>3)
  Zs_neg <- data.frame("kd"="Neg_kd","z_combined"=0)
  Zs <- Trans_Score_siRNA1_slopes_MSEmp_ind$Zs %>% 
    dplyr::mutate("kd" = paste0(Genes,"_kd")) %>% dplyr::select("kd","z_combined") 
  Zs <- rbind(Zs,Zs_neg) %>% arrange(z_combined)
  all <- all %>% left_join(Zs, by =c("Cond"="kd"))  %>% group_by(term_name) %>% dplyr::mutate(z_ab=scale(norm_ab)) %>% ungroup()
  cors_PSEA <- all %>% group_by(term_name)  %>% 
    dplyr::summarize(term_id=first(term_id),
                     cor = cor(norm_ab,z_combined),
                     intersection_size=first(intersection_size)
    ) %>% ungroup() %>% arrange(cor)
  
  filt <- all %>% filter(term_id%in%GOterms) %>% select(term_name,norm_ab,z_combined,Cond) %>% 
     left_join(cors_PSEA, by =c("term_name"="term_name")) %>% arrange(cor)
  filt$Cond <- factor(filt$Cond, levels=unique(Zs$kd))
  filt$term_name <- factor(filt$term_name, levels=unique(filt$term_name))
  
  med_data <- filt %>%
    filter(Cond %in% c("RBM6_kd", "MYBBP1A_kd", "NOP56_kd", "MDN1_kd", "BCLAF1_kd", "NUP205_kd")) %>%
    group_by(term_name) %>%
    summarize(mean_value = mean(norm_ab))
  # order by mean_value in descending order
  filt$term_name <- factor(filt$term_name, levels = med_data$term_name[order(med_data$mean_value, decreasing = TRUE)])
  
  filt <- filt %>% dplyr::mutate("med_rel_adj"=ifelse(norm_ab< -0.25,-0.25,
                                                      ifelse(norm_ab>0.25,0.25,norm_ab)))
  
  PSEA_plot_adj_med <- ggplot(filt, aes(x=Cond,y=term_name,fill=med_rel_adj)) + geom_tile(color="black") + theme_classic() + 
    scale_fill_gradientn(colors = c("#394d71", "#5f97ce", "white", "#f89b62", "#f05056"),
                         values = scales::rescale(c(min(filt$med_rel_adj), 
                                                    mean(range(filt$med_rel_adj)), 
                                                    max(filt$med_rel_adj))),
                         na.value = "grey50",
                         guide = "colourbar") + labs(fill="Log2ab") +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.y = element_blank())
  
  ggsave("Figures/Supp/siRNA_PSEA_WCs.pdf",PSEA_plot_adj_med,height=5.5,width=7.8)
  
  return(filt)
}

pD_plot_PSEA_siRNA_bulk_cors_order <- function(PSEA_KDs_WC, bulk, Trans_Score_siRNA1_slopes_MSEmp_ind, GOterms, Transport_Score_cor, KD_genes = c(dampeners,enhancers)){
  
  df_uniqueGeneTerms <- PSEA_KDs_WC %>%
    separate_rows(intersection) %>% distinct(term_id, intersection)
  uni_conds <- unique(PSEA_KDs_WC$Cond)
  empty_list <- list()
  for(i in 1:length(uni_conds)){
    temp <- PSEA_KDs_WC %>% distinct(term_id, .keep_all=T)
    bulk_DA_temp <- bulk[bulk$group1==paste0(uni_conds[i]),]
    temp$med_ab <- 0 
    temp$mean_ab <- 0 
    colnum <- ncol(temp)
    for(k in 1:nrow(temp)){
      df_uniqueGeneTerms_temp <- df_uniqueGeneTerms[df_uniqueGeneTerms$term_id==temp$term_id[k],]
      dat_temp_abs <- bulk_DA_temp[bulk_DA_temp$Gene%in%df_uniqueGeneTerms_temp$intersection,]
      dat_temp_abs_med <- median(dat_temp_abs$log2_FC) 
      dat_temp_abs_mean <- mean(dat_temp_abs$log2_FC) 
      
      temp[k,colnum-2] <- dat_temp_abs_med
      temp[k,colnum-1] <- dat_temp_abs_mean
      temp$Cond <- paste0(uni_conds[i])
    }
    empty_list[[i]] <- temp
    print(paste0(i,"/",length(uni_conds)," finished")) 
    
  }
  df <- do.call("rbind", empty_list)
  df_min0 <- df %>% distinct(term_id, .keep_all=T)
  df_min0$Cond <- "Neg_kd"
  df_min0$med_ab <- 0
  df_min0$mean_ab <- 0
  all <- rbind(df, df_min0)
  all <- all %>% dplyr::group_by(term_id) %>% dplyr::mutate("norm_ab" = mean_ab-mean(mean_ab,na.rm=T)) %>%
    dplyr::mutate(sd=sd(norm_ab)) %>% ungroup() %>% filter(intersection_size>3)
  Zs_neg <- data.frame("kd"="Neg_kd","z_combined"=0)
  Zs <- Trans_Score_siRNA1_slopes_MSEmp_ind$Zs %>% 
    dplyr::mutate("kd" = paste0(Genes,"_kd")) %>% dplyr::select("kd","z_combined") 
  Zs <- rbind(Zs,Zs_neg) %>% arrange(z_combined)
  all <- all %>% left_join(Zs, by =c("Cond"="kd"))  %>% group_by(term_name) %>% dplyr::mutate(z_ab=scale(norm_ab)) %>% ungroup()
  cors_PSEA <- all %>% group_by(term_name)  %>% 
    dplyr::summarize(term_id=first(term_id),
                     cor = cor(norm_ab,z_combined),
                     intersection_size=first(intersection_size),
                     p_val_minimum = min(p_value)
    ) %>% ungroup() %>% arrange(cor)
  
  filt <- all %>% filter(term_id%in%GOterms) %>% select(term_name,norm_ab,z_combined,Cond) %>% 
    left_join(cors_PSEA, by =c("term_name"="term_name")) %>% arrange(cor)
  
  all_lim <- all %>% select(term_name,norm_ab,z_combined,Cond) %>% 
    left_join(cors_PSEA, by =c("term_name"="term_name")) %>% arrange(cor)
  
  filt$Cond <- factor(filt$Cond, levels=unique(Zs$kd))
  filt$term_name <- factor(filt$term_name, levels=unique(filt$term_name))
  med_data <- filt %>%
    filter(Cond %in% c("RBM6_kd", "MYBBP1A_kd", "NOP56_kd", "MDN1_kd", "BCLAF1_kd", "NUP205_kd","RRS1_kd","DHX15_kd","BOP1")) %>%
    group_by(term_name) %>%
    summarize(mean_value = mean(norm_ab))
  filt$term_name <- factor(filt$term_name, levels = med_data$term_name[order(med_data$mean_value, decreasing = TRUE)])
  
  # order by mean_value in descending order
  cors <- Transport_Score_cor$cors
  cors$Genes <- paste0(cors$Genes,"_kd")
  cors_neg <- cors[1,]
  cors_neg$pr_r <- 0 #just setting the negative control to 0 for ordering reasons
  cors_neg$Genes[1] <- "Neg_kd"
  cors <- rbind(cors, cors_neg)  %>% filter(Genes%in%c(paste0(KD_genes,"_kd"),"Neg_kd"))
  filt$Cond <- factor(filt$Cond, levels = cors$Genes[order(cors$pr_r, decreasing = TRUE)])
  
  filt <- filt %>% dplyr::mutate("med_rel_adj"=ifelse(norm_ab< -0.25,-0.25,
                                                      ifelse(norm_ab>0.25,0.25,norm_ab)))
  
  PSEA_plot_adj_med <- ggplot(filt, aes(x=Cond,y=term_name,fill=med_rel_adj)) + geom_tile(color="black") + theme_classic() + 
    scale_fill_gradientn(colors = c("#394d71", "#5f97ce", "white", "#f89b62", "#f05056"),
                         values = scales::rescale(c(min(filt$med_rel_adj), 
                                                    mean(range(filt$med_rel_adj)), 
                                                    max(filt$med_rel_adj))),
                         na.value = "grey50",
                         guide = "colourbar") + labs(fill="Log2ab") +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.y = element_blank())
  
  ggsave("Figures/Supp/siRNA_PSEA_WCs_cors_ordered.pdf",PSEA_plot_adj_med,height=5.5,width=7.8)
  
  return(all_lim)
}

pD_plot_DA_WC_NUC <- function(NUC, WC, SelectCond="min_60"){
 
  NUC$Condtype <- paste0(NUC$Cond,"Nuc")
  NUC$type <- "Nuclei"
  WC$Condtype <- paste0(WC$Cond,"WC")
  WC$type <- "Whole-cell"
  result <- rbind(NUC,WC)
  result_counts  <- result %>% group_by(Condtype) %>% dplyr::mutate(total = n()) %>% 
    dplyr::mutate(sig = sum(p.adj < 0.05)) %>%
    dplyr::mutate(total_frac = total/total) %>%
    distinct(Condtype, .keep_all = TRUE) %>% ungroup() %>% select(Cond,type, total,sig) %>% melt()
  result_counts$Cond <- factor(result_counts$Cond, levels = c("min_10","min_30","min_60"))
  result_counts$y_pos <- interaction(result_counts$Cond, result_counts$variable, sep = ": ")
  
  result_counts <- result_counts %>% group_by(Cond,type) %>% dplyr::mutate("total"=max(value)) 
  difprots <- ggplot(result_counts %>% filter(!grepl("total",y_pos)), aes(y = value/total, x = y_pos, group = variable)) +
    #geom_segment(aes(xend = y_pos, yend = 0), size = 0.85) +
    geom_col(aes(fill = type),position="dodge2",color="black")+
    #geom_point(shape=21, size = 3.5, aes(fill = type),color="black") +
    ylim(0,0.18)+
    scale_fill_manual(values=c("grey","grey20"))+
    facet_wrap(~Cond, nrow = 1, scales = "free_x") +
    geom_text_repel(aes(label = round(value/total,3)), color="black", size = 3.5,
                    segment.color = 'black')+
    labs(y = "Fraction of differential proteins", x = "Condition", color = "Variable") +
    theme_bw()+  theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),legend.position = "none")
  ggsave("Figures/Bulk_difprots.pdf",difprots,width=5.75,height=4.5)
  
  #####
  result_filt <- result %>% filter(p.adj<0.05)
  range_log2FC <- max(result_filt$log2FC) - min(result_filt$log2FC)
  zero <- abs(min(result_filt$log2FC))/range_log2FC
  result <- result %>% dplyr::mutate("p.adj"=ifelse(-log10(p.adj)>15,10^-15,p.adj)) #set ceiling to 15
  result$Condtype <- factor(result$Condtype, levels = c("min_10WC","min_30WC","min_60WC","min_10Nuc","min_30Nuc","min_60Nuc"))
  
  a <- ggplot(result) + facet_wrap(~Condtype,nrow=2)+
    geom_hline(yintercept=-log10(0.05),linetype="dashed")+
    geom_point(data=result %>% filter(p.adj<0.05), aes(x=log2FC, y=-log10(p.adj), fill=log2FC), shape=21, alpha=0.9,color="gray30",stroke=0.2) +
    scale_fill_gradientn(colors=c("#394d71","#5f97ce","white","#f89b62","#f05056"),
                         values=c(0,zero-0.05,zero,zero+0.05,1), 
                         na.value="grey30", guide="colourbar") +
    geom_point(data=result %>% filter(p.adj>=0.05), aes(x=log2FC, y=-log10(p.adj)), shape=21, alpha=0.2,color="grey50",fill="grey50") +
    geom_text_repel(data = result %>% filter(p.adj < 0.05),
                    aes(x=log2FC, y=-log10(p.adj), label = Gene), color="black", size = 1.5,
                    nudge_x = 0.01,
                    nudge_y = 0.1,
                    segment.color = 'black') +
    labs(x="log2FC", y="-Log10(p.adj)") +
    theme_bw()+  theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),legend.position = "none")
  ggsave("Figures/Supp/Bulk_Volcano.pdf",a,width=6.5,height=4.3)
  
  min60 <- result %>% filter(Cond%in%SelectCond)
  result_filt <- min60 %>% filter(p.adj<0.05)
  # Calculate the range of data
  range_log2FC <- max(result_filt$log2FC) - min(result_filt$log2FC)
  zero <- abs(min(result_filt$log2FC))/range_log2FC
  b <- ggplot(min60) + facet_wrap(~Condtype,nrow=1)+
    geom_hline(yintercept=-log10(0.05),linetype="dashed")+
    geom_point(data=min60 %>% filter(p.adj<0.05), aes(x=log2FC, y=-log10(p.adj), fill=log2FC), shape=21, alpha=0.9,color="gray30",stroke=0.2) +
    scale_fill_gradientn(colors=c("#394d71","#5f97ce","white","#f89b62","#f05056"),
                         values=c(0,zero-0.05,zero,zero+0.05,1), 
                         na.value="grey30", guide="colourbar") +
    geom_point(data=min60 %>% filter(p.adj>=0.05), aes(x=log2FC, y=-log10(p.adj)), shape=21, alpha=0.2,color="grey50",fill="grey50") +
    geom_text_repel(data = min60 %>% filter(p.adj < 0.05),
                    aes(x=log2FC, y=-log10(p.adj), label = Gene), color="black", size = 3.5,
                    min.segment.length = unit(0, 'lines'),
                    segment.color = 'black') +
    labs(x="log2FC", y="-Log10(p.adj)") +
    theme_bw()+  theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),legend.position = "none")
  ggsave("Figures/Bulk_Volcano.pdf",b,width=6.5,height=4.3)
  
  return(list(a,b,difprots))
}

tls_fit <- function(x, y) {
  fit <- prcomp(~ x + y)
  intercept <- fit$center[2] - fit$rotation[2,1] / fit$rotation[1,1] * fit$center[1]
  slope <- fit$rotation[2,1] / fit$rotation[1,1]
  list(intercept = intercept, slope = slope)
}

pD_spline_correct_forLC <- function(SN_maxLFQprep_select){
  SN_maxLFQprep_select_1 <- SN_maxLFQprep_select %>% dplyr::mutate(log2val = log2(val+0.00001)) %>% #avoid log2 of zero
    dplyr::group_by(run_chan) %>% dplyr::mutate("med_log2val" = median(log2val,na.rm=T)) %>%
    dplyr::mutate("log2val"=log2val-med_log2val) %>% ungroup() 
  SN_maxLFQprep_select_1$run <- as.numeric(as.character(SN_maxLFQprep_select_1$run))
  
  results <- SN_maxLFQprep_select_1 %>%
    group_by(seqcharge) %>% 
    mutate(n = n()) %>%
    filter(n>100) %>% #require 100 precs to do the fit
    do({
      fit <- smooth.spline(.$run, .$log2val)
      predicted <- predict(fit, .$run)
      residuals <- .$log2val - predicted$y
      RSS <- sum((.$log2val - predicted$y)^2)
      TSS <- sum((.$log2val - mean(.$log2val))^2)
      R_squared <- 1 - (RSS / TSS)
      
      data.frame(run_chan=.$run_chan, run = .$run, adjusted_log2val = residuals, log2val=.$log2val, degrees=fit$df, R_squared = rep(R_squared, .$n[1]))
    }) %>%
    ungroup()
  
  results_keep <- results %>% filter(R_squared>0.2) %>%
    dplyr::select("run_chan","seqcharge","adjusted_log2val")
  
  SN_maxLFQprep_select_1 <- SN_maxLFQprep_select_1 %>% left_join(results_keep, by =c("run_chan"="run_chan","seqcharge"="seqcharge")) %>%
    dplyr::mutate("adj_val" = 2^(adjusted_log2val+med_log2val)) %>%
    dplyr::mutate("adj_val" = ifelse(is.na(adj_val)==F,adj_val,val)) #if its not NA, keep the value, else replace with the original value
    
  return(SN_maxLFQprep_select_1)
  
  
    
  
  #   
  # 
  # 
  # 
  # # Compute mean R_squared for each seqcharge (since R_squared is the same for all rows in a group, taking mean is just a formality)
  # seqcharge_rsq <- results %>%
  #   group_by(seqcharge) %>% add_count() %>% filter(n>300) %>% 
  #   summarise(R_squared = mean(R_squared)) %>%
  #   ungroup() %>%
  #   arrange(desc(R_squared))  # Sort by descending order of R_squared
  # 
  # # Identify the top 5 best and worst fits
  # best_fits <- head(seqcharge_rsq, 6)
  # worst_fits <- tail(seqcharge_rsq, 2)
  # 
  # # Filter results for best and worst fits
  # best_worst_fits <- results %>%
  #   filter(seqcharge %in% c(best_fits$seqcharge))
  # 
  # # Prepare data for plotting
  # plot_data <- best_worst_fits %>%
  #   group_by(seqcharge) %>%
  #   do({
  #     fit <- smooth.spline(.$run, .$log2val)
  #     predicted <- predict(fit, .$run)
  #     data.frame(run = .$run, log2val = .$log2val, predicted = predicted$y, degrees=fit$df, group = unique(.$seqcharge))
  #   }) %>%
  #   ungroup()
  # 
  # # Plotting the top 5 best and worst fits
  # ggplot(plot_data, aes(x = run, y = log2val, group = group)) +
  #   geom_point() +
  #  # geom_line(aes(y = log2val, colour = "Actual"), size = 1) +
  #  # geom_line(aes(y = predicted, colour = "Predicted"), size = 1, linetype = "dashed") +
  #   facet_wrap(~ group, scales = "free_y") +
  #   labs(title = "best spline fits",
  #        x = "Run", y = "Log2 Value",
  #        colour = "Type") +
  #   #theme_minimal() +
  #   scale_colour_manual(values = c("Actual" = "blue", "Predicted" = "red")) + 
  #   geom_smooth(method = lm, formula = y ~ splines::bs(x, 5), se = FALSE)
  # 
  # ggplot(plot_data, aes(x = run, y = log2val-predicted, group = group)) +
  #   geom_point() +
  #   # geom_line(aes(y = log2val, colour = "Actual"), size = 1) +
  #   # geom_line(aes(y = predicted, colour = "Predicted"), size = 1, linetype = "dashed") +
  #   facet_wrap(~ group, scales = "free_y") +
  #   labs(title = "best spline fits",
  #        x = "Run", y = "Log2 Value",
  #        colour = "Type") +
  #   #theme_minimal() +
  #   scale_colour_manual(values = c("Actual" = "blue", "Predicted" = "red")) + 
  #   geom_smooth(method = lm, formula = y ~ splines::bs(x, 4), se = FALSE)
  # 
  # hist(seqcharge_rsq$R_squared)
  
}

save_plots <- function(dat, myfilename,mywidth=7,myheight=7,myunits="in"){
  ggsave(filename=paste0(myfilename), dat,width=mywidth,height=myheight, units=myunits)
}

save_plots_list <- function(dat, myfilename,mywidth=7,myheight=7,myunits="in"){
  for(i in 1:length(dat)){
    ggsave(filename=paste0(myfilename,"_",i,".pdf"), dat[[i]],width=mywidth,height=myheight, units=myunits)
    
  }
}

pD_plot_cors_DE <- function(cors,Enhancers=c("MDN1","BOP1","RBM6"),Dampeners=c("TPR","PLEC","IMMT","VIM")){
  cors2 <- cors %>% dplyr::mutate("type_movement" = ifelse(Genes%in%Enhancers, "Enhancers", 
                                                           ifelse(Genes%in%Dampeners, "Dampeners", "other")))
  both <- c(Enhancers,Dampeners)
  plot_cors_all_static_colored <- ggplot(cors2, aes(x=Rank,y=pr_r)) + 
    geom_point(data = cors2 %>% filter(type_movement=="other"), size=0.7,alpha=0.5, color="black") +
    geom_point(data = cors2 %>% filter(type_movement=="Dampeners"), size=3,alpha=0.8, fill="#4C467E",shape=21) +
    geom_point(data = cors2 %>% filter(type_movement=="Enhancers"), size=3,alpha=0.8, fill="#95D3C9",shape=21) +
    geom_hline(yintercept=0, linetype="dashed") + theme_classic() + 
    theme(axis.title.y=element_text(size=18),
          legend.position = "none") +
    geom_label_repel(data=cors2 %>% filter(Genes%in%both),
                     aes(label=Genes),box.padding = 0.5,min.segment.length = 0,max.overlaps=Inf) + 
    labs(y="Spearman Correlation")
  ggsave("Figures/Cors_damp_enhance.pdf",plot_cors_all_static_colored,height=4.25*1.5,width=2.3*1.5)
  
  return(plot_cors_all_static_colored)
}

overlap_integral <- function(mu_A, sigma_A, mu_B, sigma_B) {
  #density functions
  dA <- function(x) dnorm(x, mean = mu_A, sd = sigma_A)
  dB <- function(x) dnorm(x, mean = mu_B, sd = sigma_B)
  
  # Define a function for the minimum of the two densities
  min_density <- function(x) pmin(dA(x), dB(x))
  
  # Adjust the integration bounds to cover a broader range if needed
  lower_bound <- min(mu_A - 4*sigma_A, mu_B - 4*sigma_B)
  upper_bound <- max(mu_A + 4*sigma_A, mu_B + 4*sigma_B)
  
  # Numerical integration of the minimum of two densities
  integral <- integrate(min_density, lower_bound, upper_bound)
  return(integral$value)
}


median.quartile <- function(x){
  out <- quantile(x, probs = c(0.5))
  names(out) <- c("y")
  return(out) 
}

pD_overlap <- function(nucs_BC_MLFQ, prot=c("NUP205","VIM")){
  m <- nucs_BC_MLFQ$no_neg
  dat <-nucs_BC_MLFQ$unimputed.BC %>% left_join(m, by =c("run_chan"="run_chan"))
  dat <- dat[!is.na(dat$norm_prot),] %>% group_by(Genes,Celltype,BioRep) %>% add_count() %>% filter(n>10)
  dat_d <- dcast(dat, Genes+run_chan+BioRep~Celltype, value.var = "norm_prot")
  
  uniprots <- unique(dat_d$Genes)
  uniBRs <- unique(dat_d$BioRep)
  
  fin_list <- list()
  for(i in 1:length(unique(dat_d$Genes))){
    specific <- dat_d %>% dplyr::filter(Genes%in%paste0(uniprots[i]))
    NUC_NT <- specific$NUC_NT %>% na.omit()
    NUC_10 <- specific$NUC_10 %>% na.omit()
    NUC_30 <- specific$NUC_30 %>% na.omit()
    NUC_60 <- specific$NUC_60 %>% na.omit()
    
    if(length(NUC_10)>10 & length(NUC_NT)>10){
      ov_10 <- overlap(NUC_NT,NUC_10)
    } else(ov_10 <- NA)
    
    if(length(NUC_30)>10 & length(NUC_NT)>10){
      ov_30 <- overlap(NUC_NT,NUC_30)
    } else(ov_30 <- NA)
    
    if(length(NUC_60)>10 & length(NUC_NT)>10){
      ov_60 <- overlap(NUC_NT,NUC_60)
    } else(ov_60 <- NA)
    
    tempdf <- data.frame(Gene=paste0(uniprots[i]),ov_10,ov_30,ov_60)
    fin_list[[i]] <- tempdf
  }
  df <- do.call("rbind",fin_list)
  df_m <- melt(df)
  
  df_m <- df_m %>% dplyr::mutate("Celltype" = ifelse(grepl("10",variable),"NUC_10",
                                                 ifelse(grepl("10",variable),"NUC_30","NUC_60")))
  OC <- ggplot(df_m, aes(x=variable,y=value)) + 
    geom_violin(fill='grey80',color="black")+#, draw_quantiles = c(0.5)) + 
    geom_boxplot(width=0.2,outlier.shape = NA)+
    #ylim(0,1) + 
    theme_classic() +
    labs(x= "LPS duration (min)", y= "Overlap coefficient")
  ggsave("Figures/OverlapCoeff_distribution.pdf",OC,height=5.1,width=3.9)
  
  my_fin_colors <- c("skyblue1","royalblue1", "midnightblue","violetred")
  
  prot=c("NUP205","VIM")
  specific <- dat %>% dplyr::filter(Genes%in%prot)
  specific_NT_10 <- specific[grepl("NT",specific$Celltype),] %>% dplyr::mutate("comparison" = "NUC_10")
  specific_NT_30 <- specific[grepl("NT",specific$Celltype),] %>% dplyr::mutate("comparison" = "NUC_30")
  specific_NT_60 <- specific[grepl("NT",specific$Celltype),] %>% dplyr::mutate("comparison" = "NUC_60")
  specific <- specific[!grepl("NT",specific$Celltype),]
  specific$comparison <- specific$Celltype
  specific <- rbind(specific,specific_NT_10,specific_NT_30,specific_NT_60)
  
  OC_ind <- ggplot(specific) + geom_density(aes(x=norm_prot,group=Celltype,fill=Celltype),alpha=0.5) +
    theme_classic() + 
    facet_wrap(~comparison+Genes,ncol=2)+
    scale_fill_manual(values=my_fin_colors)
  ggsave("Figures/OverlapCoeff_IndProts.pdf",OC_ind,height=5.1,width=5.1)
  
}


pD_SN_PSEA_plot <- function(PSEA_TS_cor, GOs=c("use these GOs")){
  
  cors <- PSEA_TS_cor %>% arrange(desc(med_ab))
  cors$Rank <- as.numeric(row.names(cors))
  cors <- cors %>% dplyr::mutate(p_value = ifelse(-log10(p_value)>100, 10^-100,p_value))
  
  plot_cor_rank <- ggplot(cors, aes(x=Rank,y=med_ab,size=-log10(p_value))) + geom_point(alpha=0.2,color="black") +
    geom_point(data = cors %>% filter(term_id%in%GOs),alpha=0.8,color="red")+
    geom_hline(yintercept=0, linetype="dashed") +
    labs(x="Rank", y="Median pearson correlation")+
    geom_label_repel(data=cors %>% filter(term_id%in%GOs),
                     aes(label=term_name),size=3,box.padding = 0.5,min.segment.length = 0,max.overlaps = Inf) + theme_bw() + 
    theme(axis.title.y=element_text(size=14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"))+
    geom_hline(yintercept=0, linetype="dashed")
  ggsave("Figures/corTS_PSEA.pdf",plot_cor_rank,height=5.5,width=4.7)
  
  return(plot_cor_rank)
}


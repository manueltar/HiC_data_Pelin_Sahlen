suppressMessages(library("plyr", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("data.table", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("crayon", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("withr", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("ggplot2", lib.loc = "/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("farver", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("labeling", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("optparse", lib.loc = "/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("dplyr", lib.loc = "/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("withr", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("backports", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("broom", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("rstudioapi", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("cli", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("tzdb", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("tidyverse", lib.loc="/nfs/team151/software/manuel_R_libs_4_1//"))
suppressMessages(library("svglite", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("ggeasy", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("sandwich", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("digest", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("reshape2", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("BiocGenerics", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("S4Vectors", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("IRanges", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("GenomeInfoDb", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("GenomicRanges", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("Biobase", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("AnnotationDbi", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("GenomicFeatures", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("OrganismDbi", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("GO.db", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("org.Hs.eg.db", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("TxDb.Hsapiens.UCSC.hg19.knownGene", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("Homo.sapiens", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("gwascat", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("rtracklayer", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("liftOver",lib.loc = "/nfs/team151/software/manuel_R_libs_4_1/"))




opt = NULL

Data_wrangling = function(option_list)
{
  
  
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("OUT\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  ### Read input variants ----
  
  Input_rds<-readRDS(opt$Input_rds)
  
  
  cat("Input_rds\n")
  str(Input_rds)
  cat("\n")
 
  #### READ and transform Cell_Type ----
  
  Cell_Type_array = unlist(strsplit(opt$SELECTED_CT, split=','))
  
  cat("Cell_Type_array\n")
  cat(str(Cell_Type_array))
  cat("\n")
 
 
  #### Read input datasets ----
  
  Gather<-data.frame()
  
  for(i in 1:length(Cell_Type_array))
  {
    Cell_Type_array_sel<-Cell_Type_array[i]
    
    cat("------->\t")
    cat(sprintf(as.character(Cell_Type_array_sel)))
    cat("\n")
    
    explore_path<-paste(out,Cell_Type_array_sel,'/',sep='')
    
    file_list <- list.files(path=explore_path, include.dirs = FALSE)
    
    cat("file_list_0\n")
    cat(str(file_list))
    cat("\n")
    
    if(length(file_list) >0)
    {
      indexes_sel <- grep("Long_matrix_findings_[^\\.]+\\.tsv$",file_list)
      
      file_list_sel <- file_list[indexes_sel]
      
      cat("file_list_sel_0\n")
      cat(str(file_list_sel))
      cat("\n")
      
      indexes_del <- grep("Long_matrix_findings_genes_[^\\.]+\\.tsv$",file_list_sel)
      
      if(length(indexes_del) >0)
      {
        file_list_sel <- file_list_sel[-indexes_del]
        
      }#length(indexes_del) >0
      
      cat("file_list_sel_0\n")
      cat(str(file_list_sel))
      cat("\n")
      
      
      
      setwd(explore_path)
      
      Long_matrix<-as.data.frame(fread(file=file_list_sel,sep="\t",header=T, stringsAsFactors = F))
      
      Long_matrix$Cell_Type<-Cell_Type_array_sel
      
      cat("Long_matrix_1\n")
      str(Long_matrix)
      cat("\n")
      
      
      #### check duplicated lines ----
      
      Long_matrix.dt<-data.table(Long_matrix, key=c("VAR","IntGroup","InteractorName","InteractorAnnotation","FeatureID","FeatureAnnotation","Cell_Type","variable","FDR"))
      
      cat("Long_matrix.dt_2\n")
      str(Long_matrix.dt)
      cat("\n")
      
      sum_DUP<-sum(duplicated(Long_matrix.dt))
      
      cat("sum_DUP_0\n")
      str(sum_DUP)
      cat("\n")
      
      
      check_DUP<-duplicated(Long_matrix.dt)[which(duplicated(Long_matrix.dt) == TRUE)]
      
      cat("check_DUP_0\n")
      str(check_DUP)
      cat("\n")
      

      Long_matrix_unique<-as.data.frame(unique(Long_matrix.dt, by=key(Long_matrix.dt)), stringsAsFactors=F)
      
      cat("Long_matrix_unique\n")
      str(Long_matrix_unique)
      cat("\n")
      
      Gather<-rbind(Long_matrix_unique,Gather)
      
    }#length(file_list) >0
    
  }#i in 1:length(Cell_Type_array)
  
 
  
  if(dim(Gather)[1] >0)
  {
    Gather$FDR<-as.character(Gather$FDR)
    
    Gather$Cell_Type<-factor(Gather$Cell_Type,
                             levels=rev(c("K562","Molm1","THP1","HEKa","CMK","GM12878")),
                             ordered=T)
    
    Gather$variable<-factor(Gather$variable,
                             levels=rev(c("Normal","CarboplatinTreated","GemcitabineTreated","LPSTreated","R1","R2","R3")),
                             ordered=T)
    
    Gather$FDR<-factor(Gather$FDR,
                            levels=c("0.1","0.01","0.001"),
                            ordered=T)
    
    
    cat("Gather\n")
    str(Gather)
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(Gather$Cell_Type))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(Gather$Cell_Type)))))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(Gather$variable))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(Gather$variable)))))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(Gather$FDR))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(Gather$FDR)))))
    cat("\n")
    

 
    
    
    # quit(status = 1)
    
    Gather_wide<-unique(as.data.frame(pivot_wider(Gather,
                                                  id_cols=c("VAR","IntGroup","InteractorName","InteractorAnnotation","FeatureID","FeatureAnnotation"),
                                                  names_from = c("Cell_Type","FDR","variable"),
                                                  names_sep='|',
                                                  values_from = value), stringsAsFactors=F))
    
    
    # #### SAVE ----
    
    if(dim(Gather_wide)[1] >0)
    {
      cat("Gather_wide_0\n")
      str(Gather_wide)
      cat("\n")
      
      Gather_wide<-merge(Input_rds,
                         Gather_wide,
                         by="VAR",
                         all.y=T)
      
      cat("Gather_wide_1\n")
      str(Gather_wide)
      cat("\n")
      
      setwd(out)
      
      write.table(Gather_wide, file=paste("Wide_matrix_findings",".tsv",sep=''),
                  sep="\t", quote = F, row.names = F)
      
    }# dim(Gather_wide)[1] >0
    
    
  }# dim(Gather)[1] >0

}



Subsetter = function(option_list)
{
  
  
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("OUT\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  
  ##### Read Wide_matrix_findings -----
  
  setwd(out)
  
  Wide_matrix_findings<-as.data.frame(fread(file=paste("Wide_matrix_findings",".tsv",sep=''),sep="\t", header = T), stringsAsFactors=F)
  
  cat("Wide_matrix_findings_0\n")
  str(Wide_matrix_findings)
  cat("\n")
  cat(str(unique(Wide_matrix_findings$VAR)))
  cat("\n")
  
  #### READ SUpp TABLE 4----
  
  
  Supp4_Table_CURATED_PLUS_PHENOTYPES<-readRDS(file=opt$Supp4_Table_CURATED_PLUS_PHENOTYPES)
  
  
  cat("Supp4_Table_CURATED_PLUS_PHENOTYPES_0:\n")
  cat(str(Supp4_Table_CURATED_PLUS_PHENOTYPES))
  cat("\n")
  cat(str(unique(Supp4_Table_CURATED_PLUS_PHENOTYPES$VAR)))
  cat("\n")
  
  indx.int<-c(which(colnames(Supp4_Table_CURATED_PLUS_PHENOTYPES) == "VAR"),which(colnames(Supp4_Table_CURATED_PLUS_PHENOTYPES) == "rs"),which(colnames(Supp4_Table_CURATED_PLUS_PHENOTYPES) == "Mechanistic_Class"),
              which(colnames(Supp4_Table_CURATED_PLUS_PHENOTYPES) == "Manual_curation"),which(colnames(Supp4_Table_CURATED_PLUS_PHENOTYPES) == "Candidate_effector"))
  
  
  Supp4_Table_subset<-unique(Supp4_Table_CURATED_PLUS_PHENOTYPES[,indx.int])
  
  
  cat("Supp4_Table_subset_0:\n")
  cat(str(Supp4_Table_subset))
  cat("\n")
  cat(str(unique(Supp4_Table_subset$VAR)))
  cat("\n")
  
  Supp4_Table_subset<-merge(Supp4_Table_subset,
                            Wide_matrix_findings,
                            by="VAR")
  
  if(dim(Supp4_Table_subset)[1] >0)
  {
    cat("Supp4_Table_subset_1\n")
    str(Supp4_Table_subset)
    cat("\n")
    cat(str(unique(Supp4_Table_subset$VAR)))
    cat("\n")
   
    
    setwd(out)
    
    write.table(Supp4_Table_subset, file=paste("Wide_matrix_findings_Supp4",".tsv",sep=''),
                sep="\t", quote = F, row.names = F)
    
  }# dim(Supp4_Table_subset)[1] >0
  
  
  
  
  #### Read Federica's file----
  
  Fedes_variants<-as.data.frame(fread(file=opt$Fedes_variants,sep=",", header=F) , stringsAsFactors=F)
  
  colnames(Fedes_variants)<-c("rsid","chr","pos38","ref","alt")
  
  Fedes_variants$chr<-paste('chr',Fedes_variants$chr,sep='')
  
  Fedes_variants$VAR_38<-paste(Fedes_variants$chr,Fedes_variants$pos38,Fedes_variants$ref,Fedes_variants$alt, sep='_')
  
  cat("Fedes_variants_0\n")
  cat(str(Fedes_variants))
  cat("\n")
  
  Fedes_variants<-merge(Fedes_variants,
                            Wide_matrix_findings,
                            by=c("VAR_38","chr","pos38","ref","alt"))
  
  if(dim(Fedes_variants)[1] >0)
  {
    cat("Fedes_variants_1\n")
    str(Fedes_variants)
    cat("\n")
    cat(str(unique(Fedes_variants$VAR_38)))
    cat("\n")
    
    
    setwd(out)
    
    write.table(Fedes_variants, file=paste("Wide_matrix_findings_Fedes",".tsv",sep=''),
                sep="\t", quote = F, row.names = F)
    
  }# dim(Fedes_variants)[1] >0
}

Data_wrangling_genes = function(option_list)
{
  
  
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("OUT\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  ### Read input variants ----
  
  Input_rds<-readRDS(opt$Input_rds)
  
  
  cat("Input_rds\n")
  str(Input_rds)
  cat("\n")
  
  #### READ and transform Cell_Type ----
  
  Cell_Type_array = unlist(strsplit(opt$SELECTED_CT, split=','))
  
  cat("Cell_Type_array\n")
  cat(str(Cell_Type_array))
  cat("\n")
  
  
  #### Read input datasets ----
  
  Gather<-data.frame()
  
  for(i in 1:length(Cell_Type_array))
  {
    Cell_Type_array_sel<-Cell_Type_array[i]
    
    cat("------->\t")
    cat(sprintf(as.character(Cell_Type_array_sel)))
    cat("\n")
    
    explore_path<-paste(out,Cell_Type_array_sel,'/',sep='')
    
    file_list <- list.files(path=explore_path, include.dirs = FALSE)
    
    cat("file_list_0\n")
    cat(str(file_list))
    cat("\n")
    
    if(length(file_list) >0)
    {
      indexes_sel <- grep("Long_matrix_findings_genes_[^\\.]+\\.tsv$",file_list)
      
      if(length(indexes_sel) >0)
      {
        
        file_list_sel <- file_list[indexes_sel]
        
        cat("file_list_sel_0\n")
        cat(str(file_list_sel))
        cat("\n")
        
        
        
        
        setwd(explore_path)
        
        Long_matrix<-as.data.frame(fread(file=file_list_sel,sep="\t",header=T, stringsAsFactors = F))
        
        Long_matrix$Cell_Type<-Cell_Type_array_sel
        
        cat("Long_matrix_1\n")
        str(Long_matrix)
        cat("\n")
        
        
        #### check duplicated lines ----
        
        Long_matrix.dt<-data.table(Long_matrix, key=c("IntGroup","InteractorName","InteractorAnnotation","FeatureID","FeatureAnnotation","Cell_Type","variable","FDR"))
        
        cat("Long_matrix.dt_2\n")
        str(Long_matrix.dt)
        cat("\n")
        
        sum_DUP<-sum(duplicated(Long_matrix.dt))
        
        cat("sum_DUP_0\n")
        str(sum_DUP)
        cat("\n")
        
        
        check_DUP<-duplicated(Long_matrix.dt)[which(duplicated(Long_matrix.dt) == TRUE)]
        
        cat("check_DUP_0\n")
        str(check_DUP)
        cat("\n")
        
        
        Long_matrix_unique<-as.data.frame(unique(Long_matrix.dt, by=key(Long_matrix.dt)), stringsAsFactors=F)
        
        cat("Long_matrix_unique\n")
        str(Long_matrix_unique)
        cat("\n")
        
        Gather<-rbind(Long_matrix_unique,Gather)
      }#length(indexes_sel) >0
    }#length(file_list) >0
    
  }#i in 1:length(Cell_Type_array)
  
  cat("Gather\n")
  str(Gather)
  cat("\n")
  
  if(dim(Gather)[1] >0)
  {
    Gather_wide<-unique(as.data.frame(pivot_wider(Gather,
                                                  id_cols=c("IntGroup","InteractorName","InteractorAnnotation","FeatureID","FeatureAnnotation"),
                                                  names_from = c("Cell_Type","FDR","variable"),
                                                  names_sep='|',
                                                  values_from = value), stringsAsFactors=F))
    
    
    # #### SAVE ----
    
    if(dim(Gather_wide)[1] >0)
    {
      cat("Gather_wide_0\n")
      str(Gather_wide)
      cat("\n")
      
     
      setwd(out)
      
      write.table(Gather_wide, file=paste("Wide_matrix_findings_genes",".tsv",sep=''),
                  sep="\t", quote = F, row.names = F)
      
    }# dim(Gather_wide)[1] >0
    
    
  }# dim(Gather)[1] >0
  
}

printList = function(l, prefix = "    ") {
  list.df = data.frame(val_name = names(l), value = as.character(l))
  list_strs = apply(list.df, MARGIN = 1, FUN = function(x) { paste(x, collapse = " = ")})
  cat(paste(paste(paste0(prefix, list_strs), collapse = "\n"), "\n"))
}


#### main script ----

main = function() {
  cmd_line = commandArgs()
  cat("Command line:\n")
  cat(paste(gsub("--file=", "", cmd_line[4], fixed=T),
            paste(cmd_line[6:length(cmd_line)], collapse = " "),
            "\n\n"))
  traceback()
  options(show.error.locations = TRUE)
  
  option_list <- list(
    make_option(c("--Input_rds"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--SELECTED_CT"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Fedes_variants"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Supp4_Table_CURATED_PLUS_PHENOTYPES"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
       make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
        make_option(c("--out"), type="character", default=NULL, 
                metavar="filename", 
                help="Path to tab-separated input file listing regions to analyze. Required.")
    
    
    
  )
  parser = OptionParser(usage = "140__Rscript_v106.R
                        --subset type
                        --TranscriptEXP FILE.txt
                        --cadd FILE.txt
                        --ncboost FILE.txt
                        --type type
                        --out filename",
                        option_list = option_list)
  opt <<- parse_args(parser)
  
 
 Data_wrangling(opt)
 Subsetter(opt)
 Data_wrangling_genes(opt)
}




###########################################################################

system.time( main() )

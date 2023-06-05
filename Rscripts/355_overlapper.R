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

OVerlapper_SNPS = function(option_list)
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
  
  test<-as.data.frame(cbind("chr10","C","G",'DUMMY',26991395,'DUMMY','DUYMMY'), stringsAsFactors=F)
  colnames(test)<-colnames(Input_rds)
  
  cat("test\n")
  str(test)
  cat("\n")
  
  # Input_rds<-rbind(Input_rds,test)
  # 
  # cat("Input_rds_POST\n")
  # str(Input_rds)
  # cat("\n")
  
    
  gr_VARS <-unique(GRanges(
    seqnames = as.character(Input_rds$chr),
    name2=rep("VAR",length(Input_rds$VAR)),
    ranges=IRanges(
      start=as.numeric(Input_rds$pos),
      end=as.numeric(Input_rds$pos),
      names = Input_rds$VAR)))
  
  cat("gr_VARS\n")
  str(gr_VARS)
  cat("\n")
  
  #### READ and transform tag ----
  
  tag = unlist(strsplit(opt$tag, split="_"))
  
  cat("tag\n")
  cat(sprintf(as.character(tag)))
  cat("\n")
  
 
  CT<-tag[1]
  
  cat("CT\n")
  cat(sprintf(as.character(CT)))
  cat("\n")
  
  FDR<-tag[2]
  
  cat("FDR\n")
  cat(sprintf(as.character(FDR)))
  cat("\n")
 
  #### Read input dataset ----
  
  explore_path<-out
  
  file_list <- list.files(path=explore_path, include.dirs = FALSE)
  
  cat("file_list_0\n")
  cat(str(file_list))
  cat("\n")
  
  indexes_sel <- grep("[^\\.]+\\.txt$",file_list)
  
  file_list_sel <- file_list[indexes_sel]
  
  cat("file_list_sel_0\n")
  cat(str(file_list_sel))
  cat("\n")
  
  indexes_sel_2 <- grep(CT,file_list_sel)
  
  file_list_sel <- file_list_sel[indexes_sel_2]
  
  cat("file_list_sel_1\n")
  cat(str(file_list_sel))
  cat("\n")
  
  indexes_sel_3 <- grep(FDR,file_list_sel)
  
  file_list_sel <- file_list_sel[indexes_sel_3]
  
  cat("file_list_sel_2\n")
  cat(str(file_list_sel))
  cat("\n")
  
  
 
  setwd(out)
  HiC<-as.data.frame(fread(file=file_list_sel,sep="\t",header=T, stringsAsFactors = F))
  
  cat("HiC\n")
  str(HiC)
  cat("\n")
  
  Condition_DEBUG <- 1
  
 #  #### Retrieve only the normal results ----
 #  
 #  Non_treated_CT<-c("K562","CMK","Molm1","HEKa","THP1")
 #  
 #  Replicates_CT<-c("GM12878")
 #  
 #  
 #  
 #  FLAG_treatment<-NA
 #  
 #  FLAG_treatment<-Non_treated_CT[which(CT%in%Non_treated_CT)]
 #  
 #  if(Condition_DEBUG == 1)
 #  {
 #    cat("FLAG_treatment_0\n")
 #    str(FLAG_treatment)
 #    cat("\n")
 # 
 #  }
 #  
 #  FLAG_replicates<-NA
 #  
 #  
 #  if(length(FLAG_treatment) == 0)
 #  {
 #    FLAG_replicates<-Replicates_CT[which(CT%in%Replicates_CT)]
 #    
 #    if(Condition_DEBUG == 1)
 #    {
 #      cat("FLAG_replicates\n")
 #      str(FLAG_replicates)
 #      cat("\n")
 #      
 #    }
 #    
 #    HiC_subset<-HiC[which(HiC$R1 == 1 &
 #                            HiC$R2 == 1),] 
 #    
 #    cat("HiC_subset\n")
 #    str(HiC_subset)
 #    cat("\n")
 #    
 #    
 #  }else{
 #    
 #    HiC_subset<-HiC[which(HiC$Normal == 1),] 
 #    
 #    cat("HiC_subset\n")
 #    str(HiC_subset)
 #    cat("\n")
 #    
 #    }#length(FLAG_treatment) == 0
 #  
 #  
 # #### Retrieve only PD and PP type of interactions ---- 
 #  
 #  HiC_subset_PDPP<-HiC_subset[which(HiC_subset$IntGroup%in%c("PD","PP")),]
 #  
 #  cat("HiC_subset_PDPP\n")
 #  str(HiC_subset_PDPP)
 #  cat("\n")
 
  HiC_subset_PDPP<-HiC
  
  gr_Interactor <-unique(GRanges(
    seqnames = as.character(HiC_subset_PDPP$Interactor_Chr),
    name2=rep("Hic_Feature_",length(HiC_subset_PDPP$InteractorID)),
    ranges=IRanges(
      start=as.numeric(HiC_subset_PDPP$Interactor_Start),
      end=as.numeric(HiC_subset_PDPP$Interactor_End),
      names = HiC_subset_PDPP$InteractorID)))
  
  
  if(Condition_DEBUG == 1)
  {
    cat("gr_Interactor\n")
    str(gr_Interactor)
    cat("\n")
    
  }
  
 
 
  # 
  ################## find regions that overlap -------------------
  
  m <- findOverlaps(gr_Interactor,gr_VARS)
  
  if(Condition_DEBUG == 1)
  {
    cat("m\n")
    str(m)
    cat("\n")
   
  }
  
  # Add gc.name to subject GRanges (i.e. left join)
  mcols(gr_VARS)$gc.name <- names(gr_VARS)
  mcols(gr_Interactor)$gc.name <- names(gr_Interactor)
  mcols(gr_VARS)[subjectHits(m), "gc.name"] = mcols(gr_Interactor)[queryHits(m), "gc.name"]
  
  
 
  df2 <- data.frame(InteractorID=gr_VARS$gc.name,
                    OV_chr=seqnames(gr_VARS),
                    OV_start=as.integer(start(gr_VARS)-1),
                    OV_end=as.integer(end(gr_VARS)),
                    VAR=names(gr_VARS), stringsAsFactors = F)
  
  cat("df2_PRE\n")
  cat(str(df2))
  cat("\n")
  
  
  indx.overlaps<-which(df2$InteractorID != df2$VAR)
  
  cat("indx.overlaps\n")
  cat(str(indx.overlaps))
  cat("\n")
  
  Selected_overlaps<-df2[indx.overlaps,]
  
  cat("Selected_overlaps_0\n")
  cat(str(Selected_overlaps))
  cat("\n")
  cat(str(unique(Selected_overlaps$VAR)))
  cat("\n")
  
  
  
  # #### SAVE ----
  
  if(dim(Selected_overlaps)[1] >0)
  {
    path5<-paste(out,CT,'/', sep='')
    
    if (file.exists(path5)){
      
      
      
      
    } else {
      dir.create(file.path(path5))
      
    }
    
    path6<-paste(path5,FDR,'/', sep='')
    
    if (file.exists(path6)){
      
      
      
      
    } else {
      dir.create(file.path(path6))
      
    }
    
    Selected_overlaps<-merge(Selected_overlaps,
                             HiC_subset_PDPP,
                             by="InteractorID")
    
    cat("Selected_overlaps_1\n")
    cat(str(Selected_overlaps))
    cat("\n")
    cat(str(unique(Selected_overlaps$VAR)))
    cat("\n")
    
    # ################
    # quit(status = 1)
    
    cat("Hello_world\t")
    cat(sprintf(as.character(tag)))
    cat("\n")
    
    setwd(path6)
    
    write.table(Selected_overlaps, file=paste("Overlaps_",paste(tag, collapse="_"),".tsv",sep=''),
                sep="\t", quote = F, row.names = F)
    
  }# dim(Selected_overlaps)[1] >0

  #file=paste("overlaps_",tag,".tsv", sep=''),
}

Overlapper_genes = function(option_list)
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
  
  input_gene_names<-unlist(strsplit(opt$input_gene_names, split=","))
  
  
  cat("input_gene_names\n")
  str(input_gene_names)
  cat("\n")
  
  
  
  #### READ and transform tag ----
  
  tag = unlist(strsplit(opt$tag, split="_"))
  
  cat("tag\n")
  cat(sprintf(as.character(tag)))
  cat("\n")
  
  
  CT<-tag[1]
  
  cat("CT\n")
  cat(sprintf(as.character(CT)))
  cat("\n")
  
  FDR<-tag[2]
  
  cat("FDR\n")
  cat(sprintf(as.character(FDR)))
  cat("\n")
  
  #### Read input dataset ----
  
  explore_path<-out
  
  file_list <- list.files(path=explore_path, include.dirs = FALSE)
  
  cat("file_list_0\n")
  cat(str(file_list))
  cat("\n")
  
  indexes_sel <- grep("[^\\.]+\\.txt$",file_list)
  
  file_list_sel <- file_list[indexes_sel]
  
  cat("file_list_sel_0\n")
  cat(str(file_list_sel))
  cat("\n")
  
  indexes_sel_2 <- grep(CT,file_list_sel)
  
  file_list_sel <- file_list_sel[indexes_sel_2]
  
  cat("file_list_sel_1\n")
  cat(str(file_list_sel))
  cat("\n")
  
  indexes_sel_3 <- grep(FDR,file_list_sel)
  
  file_list_sel <- file_list_sel[indexes_sel_3]
  
  cat("file_list_sel_2\n")
  cat(str(file_list_sel))
  cat("\n")
  
  
  
  setwd(out)
  HiC<-as.data.frame(fread(file=file_list_sel,sep="\t",header=T, stringsAsFactors = F))
  
  cat("HiC\n")
  str(HiC)
  cat("\n")
  
  Selected_overlaps<-HiC[which(HiC$InteractorName%in%input_gene_names),]
  
  cat("Selected_overlaps_0\n")
  cat(str(Selected_overlaps))
  cat("\n")
  
  
  if(dim(Selected_overlaps)[1] >0)
  {
    path5<-paste(out,CT,'/', sep='')
    
    if (file.exists(path5)){
      
      
      
      
    } else {
      dir.create(file.path(path5))
      
    }
    
    path6<-paste(path5,FDR,'/', sep='')
    
    if (file.exists(path6)){
      
      
      
      
    } else {
      dir.create(file.path(path6))
      
    }
    
   
    cat("Selected_overlaps_1\n")
    cat(str(Selected_overlaps))
    cat("\n")
    
    setwd(path6)
    
    write.table(Selected_overlaps, file=paste("Overlaps_genes_",paste(tag, collapse="_"),".tsv",sep=''),
                sep="\t", quote = F, row.names = F)
    
  
    cat("\n")
  } # dim(Selected_overlaps)[1] >0
    
    
 

  
  Condition_DEBUG <- 1
  
  
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
    make_option(c("--input_gene_names"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--tag"), type="character", default=NULL,
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
  
 
  OVerlapper_SNPS(opt)
  Overlapper_genes(opt)
  
}




###########################################################################

system.time( main() )

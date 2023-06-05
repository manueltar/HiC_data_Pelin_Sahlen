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

Put_together_SNPS = function(option_list)
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
  
 
  #### READ and transform Cell_Type ----
  
  Cell_Type = opt$Cell_Type
  
  cat("Cell_Type\n")
  cat(sprintf(as.character(Cell_Type)))
  cat("\n")
 
 
  #### Read input dataset ----
  
  explore_path<-paste(out,Cell_Type,'/',sep='')
  
  dir_list <- list.dirs(path=explore_path, recursive = FALSE)
  
  cat("dir_list_0\n")
  cat(str(dir_list))
  cat("\n")
  
  Gather<-data.frame()
  
  if(length(dir_list) >0)
  {
    for(i in 1:length(dir_list))
    {
      dir_list_sel<-dir_list[i]
      
      cat("------->\t")
      cat(sprintf(as.character(dir_list_sel)))
      cat("\n")
      
      FDR_sel<-gsub(out,"",dir_list_sel)
      
      cat("---FDR1---->\t")
      cat(sprintf(as.character(FDR_sel)))
      cat("\n")
      
      FDR_sel<-gsub(paste(Cell_Type,'/',sep=''),"",FDR_sel)
      
      cat("---FDR2---->\t")
      cat(sprintf(as.character(FDR_sel)))
      cat("\n")
      
      FDR_sel<-gsub(paste('/','FDR',sep=''),"",FDR_sel)
      
      cat("---FDR3---->\t")
      cat(sprintf(as.character(FDR_sel)))
      cat("\n")
      
      file_list <- list.files(path=dir_list_sel, include.dirs = FALSE)
      
      cat("file_list_0\n")
      cat(str(file_list))
      cat("\n")
      
      indexes_sel <- grep("[^\\.]+\\.tsv$",file_list)
      
      file_list_sel <- file_list[indexes_sel]
      
      cat("file_list_sel_0\n")
      cat(str(file_list_sel))
      cat("\n")
      
      indexes_del <- grep("genes",file_list_sel)
      
      if(length(indexes_del) >0)
      {
        file_list_sel <- file_list_sel[-indexes_del]
      }#length(indexes_del) >0
      else{
        
      }
      
      
      cat("file_list_sel_0.5\n")
      cat(str(file_list_sel))
      cat("\n")
      
      
      setwd(dir_list_sel)
      
      intersections<-as.data.frame(fread(file=file_list_sel,sep="\t",header=T, stringsAsFactors = F))
      
      cat("intersections\n")
      str(intersections)
      cat("\n")
      
      indx.int<-c(which(colnames(intersections) == "VAR"),which(colnames(intersections) == "IntGroup"),
                  which(colnames(intersections) == "InteractorName"),which(colnames(intersections) == "InteractorAnnotation"),
                  which(colnames(intersections) == "RefSeqName"),which(colnames(intersections) == "Annotation"))
      
      cat("indx.int\n")
      str(indx.int)
      cat("\n")
      
      anchor_point_1<-which(colnames(intersections) == "IntGroup")+1
      
      
      cat("anchor_point_1\n")
      str(anchor_point_1)
      cat("\n")
      
      anchor_point_2<-dim(intersections)[2]-1
      
      
      cat("anchor_point_2\n")
      str(anchor_point_2)
      cat("\n")
      
      indx.int2<-c(anchor_point_1:anchor_point_2)
      
      cat("indx.int2\n")
      str(indx.int2)
      cat("\n")
      
      intersections_subset<-unique(intersections[,c(indx.int,indx.int2)])
      
      cat("intersections_subset_0\n")
      str(intersections_subset)
      cat("\n")
      
      colnames(intersections_subset)[which(colnames(intersections_subset) == "RefSeqName")]<-"FeatureID"
      colnames(intersections_subset)[which(colnames(intersections_subset) == "Annotation")]<-"FeatureAnnotation"
      
      cat("intersections_subset_1\n")
      str(intersections_subset)
      cat("\n")
      
      anchor_point_1<-which(colnames(intersections_subset) == "FeatureAnnotation")+1
      
      
      cat("anchor_point_1\n")
      str(anchor_point_1)
      cat("\n")
      
      anchor_point_2<-dim(intersections_subset)[2]
      
      
      cat("anchor_point_2\n")
      str(anchor_point_2)
      cat("\n")
      
      indx.int3<-c(anchor_point_1:anchor_point_2)
      
      cat("indx.int3\n")
      str(indx.int3)
      cat("\n")
      
      intersections_subset.m<-melt(intersections_subset, id.vars=colnames(intersections_subset)[1:which(colnames(intersections_subset) == "FeatureAnnotation")])
      
      intersections_subset.m$variable<-as.character(intersections_subset.m$variable)
      
      intersections_subset.m$FDR<-FDR_sel
            
      cat("intersections_subset.m_0\n")
      str(intersections_subset.m)
      cat("\n")
      
      Gather<-rbind(Gather,intersections_subset.m)
      
      
    }#i in 1:length(dir_list)
    
  }# length(dir_list) >0
  
  cat("Gather_0\n")
  str(Gather)
  cat("\n")
  
 
  
  Condition_DEBUG <- 1
  

  
  
  
  # #### SAVE ----
  
  if(dim(Gather)[1] >0)
  {
    path5<-paste(out,Cell_Type,'/', sep='')
    
    if (file.exists(path5)){
      
      
      
      
    } else {
      dir.create(file.path(path5))
      
    }
    
    
    
    setwd(path5)
    
    write.table(Gather, file=paste("Long_matrix_findings_",Cell_Type,".tsv",sep=''),
                sep="\t", quote = F, row.names = F)
    
  }# dim(Gather)[1] >0

  #file=paste("overlaps_",Cell_Type,".tsv", sep=''),
}


Put_together_genes = function(option_list)
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
  
  
  #### READ and transform Cell_Type ----
  
  Cell_Type = opt$Cell_Type
  
  cat("Cell_Type\n")
  cat(sprintf(as.character(Cell_Type)))
  cat("\n")
  
  
  #### Read input dataset ----
  
  explore_path<-paste(out,Cell_Type,'/',sep='')
  
  dir_list <- list.dirs(path=explore_path, recursive = FALSE)
  
  cat("dir_list_0\n")
  cat(str(dir_list))
  cat("\n")
  
  Gather<-data.frame()
  
  if(length(dir_list) >0)
  {
    for(i in 1:length(dir_list))
    {
      dir_list_sel<-dir_list[i]
      
      cat("------->\t")
      cat(sprintf(as.character(dir_list_sel)))
      cat("\n")
      
      FDR_sel<-gsub(out,"",dir_list_sel)
      
      cat("---FDR1---->\t")
      cat(sprintf(as.character(FDR_sel)))
      cat("\n")
      
      FDR_sel<-gsub(paste(Cell_Type,'/',sep=''),"",FDR_sel)
      
      cat("---FDR2---->\t")
      cat(sprintf(as.character(FDR_sel)))
      cat("\n")
      
      FDR_sel<-gsub(paste('/','FDR',sep=''),"",FDR_sel)
      
      cat("---FDR3---->\t")
      cat(sprintf(as.character(FDR_sel)))
      cat("\n")
      
      file_list <- list.files(path=dir_list_sel, include.dirs = FALSE)
      
      cat("file_list_0\n")
      cat(str(file_list))
      cat("\n")
      
      indexes_sel <- grep("genes",file_list)
      
      if(length(indexes_sel) >0)
      {
        file_list_sel <- file_list[indexes_sel]
        
        cat("file_list_sel_0\n")
        cat(str(file_list_sel))
        cat("\n")
        
        
        
        setwd(dir_list_sel)
        
        intersections<-as.data.frame(fread(file=file_list_sel,sep="\t",header=T, stringsAsFactors = F))
        
        cat("intersections\n")
        str(intersections)
        cat("\n")
        
        indx.int<-c(which(colnames(intersections) == "VAR"),which(colnames(intersections) == "IntGroup"),
                    which(colnames(intersections) == "InteractorName"),which(colnames(intersections) == "InteractorAnnotation"),
                    which(colnames(intersections) == "RefSeqName"),which(colnames(intersections) == "Annotation"))
        
        cat("indx.int\n")
        str(indx.int)
        cat("\n")
        
        anchor_point_1<-which(colnames(intersections) == "IntGroup")+1
        
        
        cat("anchor_point_1\n")
        str(anchor_point_1)
        cat("\n")
        
        anchor_point_2<-dim(intersections)[2]-1
        
        
        cat("anchor_point_2\n")
        str(anchor_point_2)
        cat("\n")
        
        indx.int2<-c(anchor_point_1:anchor_point_2)
        
        cat("indx.int2\n")
        str(indx.int2)
        cat("\n")
        
        intersections_subset<-unique(intersections[,c(indx.int,indx.int2)])
        
        cat("intersections_subset_0\n")
        str(intersections_subset)
        cat("\n")
        
        colnames(intersections_subset)[which(colnames(intersections_subset) == "RefSeqName")]<-"FeatureID"
        colnames(intersections_subset)[which(colnames(intersections_subset) == "Annotation")]<-"FeatureAnnotation"
        
        cat("intersections_subset_1\n")
        str(intersections_subset)
        cat("\n")
        
        anchor_point_1<-which(colnames(intersections_subset) == "FeatureAnnotation")+1
        
        
        cat("anchor_point_1\n")
        str(anchor_point_1)
        cat("\n")
        
        anchor_point_2<-dim(intersections_subset)[2]
        
        
        cat("anchor_point_2\n")
        str(anchor_point_2)
        cat("\n")
        
        indx.int3<-c(anchor_point_1:anchor_point_2)
        
        cat("indx.int3\n")
        str(indx.int3)
        cat("\n")
        
        intersections_subset.m<-melt(intersections_subset, id.vars=colnames(intersections_subset)[1:which(colnames(intersections_subset) == "FeatureAnnotation")])
        
        intersections_subset.m$variable<-as.character(intersections_subset.m$variable)
        
        intersections_subset.m$FDR<-FDR_sel
        
        cat("intersections_subset.m_0\n")
        str(intersections_subset.m)
        cat("\n")
        
        Gather<-rbind(Gather,intersections_subset.m)
        
      }#length(indexes_sel) >0
      
    
      
      
    }#i in 1:length(dir_list)
    
  }# length(dir_list) >0
  
  cat("Gather_0\n")
  str(Gather)
  cat("\n")
  
  
  
  Condition_DEBUG <- 1
  
  
  
  
  
  # #### SAVE ----
  
  if(dim(Gather)[1] >0)
  {
    path5<-paste(out,Cell_Type,'/', sep='')
    
    if (file.exists(path5)){
      
      
      
      
    } else {
      dir.create(file.path(path5))
      
    }
    
    
    
    setwd(path5)
    
    write.table(Gather, file=paste("Long_matrix_findings_genes_",Cell_Type,".tsv",sep=''),
                sep="\t", quote = F, row.names = F)
    
  }# dim(Gather)[1] >0
  
  #file=paste("overlaps_",Cell_Type,".tsv", sep=''),
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
    make_option(c("--Cell_Type"), type="character", default=NULL,
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
  
 
 Put_together_SNPS(opt)
 Put_together_genes(opt)
  
}




###########################################################################

system.time( main() )

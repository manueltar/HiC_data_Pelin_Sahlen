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
  
 
  
  #### Read Input_list_of_variants file ----
  
  Input_list_of_variants<-as.data.frame(fread(file=opt$Input_list_of_variants,sep=",", header=F) , stringsAsFactors=F)
  
  colnames(Input_list_of_variants)<-c("rsid","chr","pos38","ref","alt")
  
  Input_list_of_variants$chr<-paste('chr',Input_list_of_variants$chr,sep='')
  
  Input_list_of_variants$VAR_38<-paste(Input_list_of_variants$chr,Input_list_of_variants$pos38,Input_list_of_variants$ref,Input_list_of_variants$alt, sep='_')
  
  cat("Input_list_of_variants_0\n")
  cat(str(Input_list_of_variants))
  cat("\n")
 
 
  
  VAR_38_df<-unique(data.frame(chr=Input_list_of_variants$chr,
                            pos38=Input_list_of_variants$pos38,
                            ref=Input_list_of_variants$ref,
                            alt=Input_list_of_variants$alt,
                            VAR_38=Input_list_of_variants$VAR_38,
                            stringsAsFactors = F))
  
  Condition_DEBUG <- 1
  
  if(Condition_DEBUG == 1)
  {
    cat("VAR_df_\n")
    str(VAR_38_df)
    cat("\n")
    cat(str(unique(VAR_38_df$VAR_38)))
    cat("\n")
  }
  
  gr_VARS <-unique(GRanges(
    seqnames = as.character(gsub("^chr","",VAR_38_df$chr)),
    name2=rep("VAR",length(VAR_38_df$VAR_38)),
    ranges=IRanges(
      start=as.numeric(VAR_38_df$pos38),
      end=as.numeric(VAR_38_df$pos38),
      names = VAR_38_df$VAR_38)))
  
  
  if(Condition_DEBUG == 1)
  {
    cat("gr_VARS_\n")
    str(gr_VARS)
    cat("\n")
   
  }
 
  
  ##### LiftOver 38_to_37 ----
  
  path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
  ch = import.chain(path)
  
  seqlevelsStyle(gr_VARS) = "UCSC"  # necessary
  gr_VARS19 = liftOver(gr_VARS, ch)
  gr_VARS19 = unlist(gr_VARS19)
  genome(gr_VARS19) = "hg19"
  
  if(Condition_DEBUG == 1)
  {
    cat("------------->gr_VARS19_\n")
    str(gr_VARS19)
    cat("\n")
  }
  
 
  
  if(length(gr_VARS19) >0)
  {
    
    chr_19<-as.character(seqnames(gr_VARS19))
    names_19<-as.character(names(gr_VARS19))
    
    ref_VAR19<-gsub("^chr[^_]+_[0-9]+_","",names_19)
    ref_VAR19<-gsub("_.+$","",ref_VAR19)
    
    if(Condition_DEBUG == 1)
    {
      cat("ref_VAR19\n")
      cat(sprintf(as.character(ref_VAR19)))
      cat("\n")
    }
    
    alt_VAR19<-gsub("^chr[^_]+_[0-9]+_[^_]+_","",names_19)
    # alt_VAR19<-gsub("_.+$","",alt_VAR19)
    
    if(Condition_DEBUG == 1)
    {
      cat("alt_VAR19\n")
      cat(sprintf(as.character(alt_VAR19)))
      cat("\n")
    }
    
    
    VAR_df<-data.frame(chr=as.character(seqnames(gr_VARS19)),
                       pos=start(gr_VARS19),
                       ref=ref_VAR19,
                       alt=alt_VAR19,
                       VAR_38=names(gr_VARS19),
                       stringsAsFactors = F)
    
   
    VAR_df$VAR<-paste(VAR_df$chr,VAR_df$pos,VAR_df$ref,VAR_df$alt,sep='_')
    
    if(Condition_DEBUG == 1)
    {
      cat("VAR_df_0\n")
      str(VAR_df)
      cat("\n")
    }
    
    
    if(Condition_DEBUG == 1)
    {
      cat("-------REMEMBER ------>VAR38_df_\n")
      str(VAR_38_df)
      cat("\n")
    }
    
    VAR_df<-unique(merge(VAR_df,
                             VAR_38_df,
                             by=c("chr","ref","alt","VAR_38"),
                             all=T))
    
    VAR_df$VAR[is.na(VAR_df$VAR)]<-"ABSENT"
    
    if(Condition_DEBUG == 1)
    {
        cat("VAR_DEF_df_POST_MERGE\n")
        str(VAR_df)
        cat("\n")
    }
  
    
    
    check.ABSENT<-VAR_df[which(VAR_df$VAR == "ABSENT"),]
    #
    if(Condition_DEBUG == 1)
    {
      cat("check.ABSENT\n")
      str(check.ABSENT) # 151 absent variants
      cat("\n")
    }
    
    VAR_df<-VAR_df[which(VAR_df$VAR != "ABSENT"),]
    
    if(Condition_DEBUG == 1)
    {
      cat("VAR_DEF_df_NO ABSENT\n")
      str(VAR_df)
      cat("\n")
    }
    
    
  }else{
    
    
    stop("NO_LIFT_OVER\n")
    
  }# length(gr_VARS19) >0 
  
  
  #### Read ALL_dB file ----
  
  ALL_dB<-as.data.frame(fread(file=opt$ALL_dB,sep="\t") , stringsAsFactors=F)
  
  cat("ALL_dB_0\n")
  cat(str(ALL_dB))
  cat("\n")
  
  
  indx.int<-c(which(colnames(ALL_dB) == "chr"),which(colnames(ALL_dB) == "ref"),which(colnames(ALL_dB) == "alt"),which(colnames(ALL_dB) == "pos37"),which(colnames(ALL_dB) == "VAR"))
  
  ALL_dB_subset<-unique(ALL_dB[,indx.int])
  
  colnames(ALL_dB_subset)[which(colnames(ALL_dB_subset) == "pos37")]<-"pos"
  
  cat("ALL_dB_subset\n")
  cat(str(ALL_dB_subset))
  cat("\n")
  
  ###### LiftOver 37_to_38 ----
  
  
  gr_VARS <- GRanges(
    seqnames = as.character(gsub("chr","",ALL_dB_subset$chr)),
    ranges=IRanges(
      start=as.numeric(ALL_dB_subset$pos),
      end=as.numeric(ALL_dB_subset$pos),
      name=ALL_dB_subset$VAR))
  
  # cat("gr_VARS\n")
  # str(gr_VARS)
  # cat("\n")
  
  VAR_df2<-data.frame(chr=as.character(paste('chr',seqnames(gr_VARS), sep='')),
                     pos=start(gr_VARS),
                     ref=ALL_dB_subset$ref,
                     alt=ALL_dB_subset$alt,
                     VAR=ALL_dB_subset$VAR,
                     stringsAsFactors = F)
  
  if(Condition_DEBUG == 1)
  {
    cat("VAR_df2_\n")
    str(VAR_df2)
    cat("\n")
  }
  
  #path = system.file(package="liftOver", "extdata", "Hg19Tohg38.over.chain")
  ch = import.chain("/nfs/team151/software/manuel_R_ext_data_4_1/hg19ToHg38.over.chain")
  
  seqlevelsStyle(gr_VARS) = "UCSC"  # necessary
  gr_VARS38 = liftOver(gr_VARS, ch)
  gr_VARS38 = unlist(gr_VARS38)
  genome(gr_VARS38) = "hg38"
  
  if(length(gr_VARS38) >0)
  {
    
    chr_38<-as.character(seqnames(gr_VARS38))
    names_38<-as.character(names(gr_VARS38))
    
    ref_VAR38<-gsub("^chr[^_]+_[0-9]+_","",names_38)
    ref_VAR38<-gsub("_.+$","",ref_VAR38)
    
    
    # cat("ref_VAR38\n")
    # cat(sprintf(as.character(ref_VAR38)))
    # cat("\n")
    
    alt_VAR38<-gsub("^chr[^_]+_[0-9]+_[^_]+_","",names_38)
    # alt_VAR38<-gsub("_.+$","",alt_VAR38)
    
    
    # cat("alt_VAR38\n")
    # cat(sprintf(as.character(alt_VAR38)))
    # cat("\n")
    
    
    
    
    VAR_38_df<-data.frame(chr=as.character(seqnames(gr_VARS38)),
                          pos38=start(gr_VARS38),
                          ref=ref_VAR38,
                          alt=alt_VAR38,
                          VAR=names(gr_VARS38),
                          stringsAsFactors = F)
    
    VAR_38_df$VAR_38<-paste(VAR_38_df$chr,VAR_38_df$pos38,VAR_38_df$ref,VAR_38_df$alt,sep='_')
    
    if(Condition_DEBUG == 1)
    {
      cat("VAR_38_df_1\n")
      str(VAR_38_df)
      cat("\n")
    }
    
    
    VAR_df2<-unique(merge(VAR_df2,
                                                           VAR_38_df,
                                                           by=c("chr","ref","alt","VAR"),
                                                           all=T))
    
    VAR_df2$VAR_38[is.na(VAR_df2$VAR_38)]<-"ABSENT"
    
    if(Condition_DEBUG == 1)
    {
      cat("VAR_df2_2\n")
      str(VAR_df2)
      cat("\n")
    }
    
    check.ABSENT<-VAR_df2[which(VAR_df2$VAR_38 == "ABSENT"),]
    #
    if(Condition_DEBUG == 1)
    {
      cat("check.ABSENT\n")
      str(check.ABSENT) #90
      cat("\n")
    }
    
    VAR_df2<-VAR_df2[which(VAR_df2$VAR_38 != "ABSENT"),]
    
    if(Condition_DEBUG == 1)
    {
      cat("VAR_DEF_df_NO ABSENT\n")
      str(VAR_df2)
      cat("\n")
    }
    
    
    
   
    
  }else{
    
    
    stop("NO_LIFT_OVER\n")
    
  }# length(gr_VARS38) >0
  
  
 #### rbind ----
  
  DEF<-rbind(VAR_df,VAR_df2)
  
  cat("DEF\n")
  str(DEF)
  cat("\n")
  
  ##### save -----
  
  setwd(out)
  
  saveRDS(file="conversion.rds",DEF)
   
  write.table(DEF, file=paste("conversion",".tsv",sep=''),
              sep="\t", quote = F, row.names = F)

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
    make_option(c("--Input_list_of_variants"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ALL_dB"), type="character", default=NULL,
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
  
}




###########################################################################

system.time( main() )

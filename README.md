
# Bash script command

 bash ~/Scripts/Wraper_scripts/316_Pelin_Sahlen_HiC_capture_intersect.sh /lustre/scratch126/humgen/teams/soranzo/users/mt19/HiC_Sahlen/ CCR2,EEFSEC,GATA2 4000 1 normal

# Rclone commands 

rclone copy Wide_matrix_findings_Fedes.tsv mt19_g_drive:/Project_WetLab_Projects/HiC_data_Sahlen/ --drive-shared-with-me;rclone copy Wide_matrix_findings_genes.tsv mt19_g_drive:/Project_WetLab_Projects/HiC_data_Sahlen/ --drive-shared-with-me;rclone copy Wide_matrix_findings_Supp4.tsv mt19_g_drive:/Project_WetLab_Projects/HiC_data_Sahlen/ --drive-shared-with-me;rclone copy Wide_matrix_findings.tsv mt19_g_drive:/Project_WetLab_Projects/HiC_data_Sahlen/ --drive-shared-with-me


# This HiC data is obtained through a collaboration with Pelin Shalen SciLifeLab Stockholm

# Datasets

  K562 datasets:

  K562.hg19.AllInteractions.SP4.FDR0.001.txt
  K562.hg19.AllInteractions.SP4.FDR0.01.txt
  K562.hg19.AllInteractions.SP4.FDR0.1.txt

  - Different threshold of stringency in the analysis FDR 0.1, 0.01 and 0.001

  - Header. Always values refer to FDR 0.1
  
	-RefSeqName: Region contacted name. Normally a gene name but also rsid and chr:pos. rsid contacts (28846) are classified as IntGroup DD. chr contacts (4229) are classified as IntGroup DD. 13582 genes protein coding and non-coding. Gene with most interactions NUP214 (1614). Genes with interactions are classified as  205617 PD and 56790 PP. PP interactions have a gene and transript name in these fields and in the Interactor fields

	-TranscriptName: The majority (12706) 1 transcript per gene and usually the RefSeq Id but sometimes ensembl (ENST with point of version). Some genes (876) like FOXP1 have more than one transcript. For rsid and chr:pos it is a repetition of the RefSeqName.

11 FOXP1-IT1       ENST00000498714.1
1 FOXP1   NM_001244808
1 FOXP1   NM_001244812
569 FOXP1   NM_001244813
2 FOXP1   NM_001244814
1 FOXP1   NM_001244815
1 FOXP1   NM_001244816
6 FOXP1   NM_001349342

	-Feature_Chr: autosomes and chromosome X
	-Feature_Start: start of the transcripts. GRCh37
	-Annotation: 262407 1 for genes protein coding and non-coding. 33075 2 for rsid and chrom:pos. 1 Promoter and 2 targeted SNP or candidate enhnacer region
	-Strand: for genes protein coding and non-coding 130651 + and 131756 - . for rsid and chrom:pos 33075 +
	-InteractorName: Again rsid chrom and gene data for the interactor. Gene data  56790 PP. rsid data some distal elements identified with the rsid that intersects the region 4861 PD. for the chr  33075 DD  200756 PD

	-InteractorID: rsid and transcript id for the 4861 PD and the 56790 PP. chrom:pos for for the rest
	-Interactor_Chr: chromosome including chrX, chrU and chrM
	-Interactor_Start: INT. GRCh37
	-Interactor_End: INT. GRCh37
	-InteractorAnnotation: 1-> Promoter (56790)  2-> Targeted SNP (5638)  3->Distal region (233054) 


	-distance: Interactor_Start-Feature_Start. INT. Positive if the feature is downstream and negative if it is upstream
	-KG1_SuppPairs: INT 1-200 more frequent in the first numbers. KG1 -> K562 with Gemcitabine replica 1
	-KG2_SuppPairs: INT 1-200 more frequent in the first numbers. KG2 -> K562 with Gemcitabine replica 2
	-KC1_SuppPairs: INT 1-200 more frequent in the first numbers. KC1 -> K562 with Carboplatin replica 1
	-KC2_SuppPairs: INT 1-200 more frequent in the first numbers. KC2 -> K562 with Carboplatin replica 2
	-KN1_SuppPairs: INT 1-200 more frequent in the first numbers. KN1 -> K562 non-treated replica 1
	-KN2_SuppPairs: INT 1-200 more frequent in the first numbers. KN2 -> K562 non-treated replica 2
20	-KG1_p_value: FLOAT pvalue  KG1 -> K562 with Gemcitabine replica 1 # The pvalues are not corrected (so nominal) but the filtering was done using the adjusted pvals. I understand that at least 1 has to be significant
	-KG2_p_value: FLOAT pvalue  KG2 -> K562 with Gemcitabine replica 2
	-KC1_p_value: FLOAT pvalue  KC1 -> K562 with Carboplatin replica 1
	-KC2_p_value: FLOAT pvalue  KC2 -> K562 with Carboplatin replica 2
	-KN1_p_value: FLOAT pvalue  KN1 -> K562 with non-treated replica 1
25	-KN2_p_value: FLOAT pvalue  KN2 -> K562 with non-treated replica 2
	-IntGroup: PD (Promoter, Distal), PP (Promoter,Promoter) and  DD (Distal,Distal)
		   205617 PD
		    56790 PP
		    33075 DD

	-Normal: LOGICAL. 1-> significant in both samples of the non-treated 0-> non-significant in at least one of the two
	-CarboplatinTreated:  LOGICAL. 1-> significant in both samples of the Carboplatin 0-> non-significant in at least one of the two
	-GemcitabineTreated: LOGICAL. 1-> significant in both samples of the Gemcitabine 0-> non-significant in at least one of the two
30	-NofInts: Number of samples with the signficant interaction:

	183273 1
         63984 2
         48225 3


- CMK datasets:

  CMK.hg19.AllInteractions.SP4.FDR0.001.txt
  CMK.hg19.AllInteractions.SP4.FDR0.01.txt
  CMK.hg19.AllInteractions.SP4.FDR0.1.txt

  - Different threshold of stringency in the analysis FDR 0.1, 0.01 and 0.001

  - Header. Always values refer to FDR 0.1

  RefSeqName      TranscriptName  Feature_Chr     Feature_Start
  Annotation      Strand
  InteractorName  InteractorID    Interactor_Chr  Interactor_Start        Interactor_End  InteractorAnnotation    distance
  CG1_SuppPairs   CG2_SuppPairs   CC1_SuppPairs   CC2_SuppPairs   CN1_SuppPairs   CN2_SuppPairs # CG1 -> CMK treated with Gemcitabine replica1
  CG1_p_value  CG2_p_value     CC1_p_value     CC2_p_value     CN1_p_value     CN2_p_value
  IntGroup
  Normal  CarboplatinTreated      GemcitabineTreated      NofInts

- Molm1 datasets

Molm1.hg19.AllInteractions.SP4.FDR0.001.txt  Molm1.hg19.AllInteractions.SP4.FDR0.01.txt  Molm1.hg19.AllInteractions.SP4.FDR0.1.txt

  - Different threshold of stringency in the analysis FDR 0.1, 0.01 and 0.001

  - Header. Always values refer to FDR 0.1

RefSeqName      TranscriptName
Feature_Chr     Feature_Start   Annotation      Strand
InteractorName  InteractorID    Interactor_Chr  Interactor_Start        Interactor_End  InteractorAnnotation    distance
MG1_SuppPairs   MG2_SuppPairs   MC1_SuppPairs   MC2_SuppPairs   MN1_SuppPairs   MN2_SuppPairs #MG1 -> Molm1 treated with Gemcitabine replica1
MG1_p_value  MG2_p_value     MC1_p_value     MC2_p_value     MN1_p_value     MN2_p_value
IntGroup
Normal  CarboplatinTreated      GemcitabineTreated      NofInts

- THP1 datasets

THP1run.hg19.AllInteractions.SP4.FDR0.001.txt  THP1run.hg19.AllInteractions.SP4.FDR0.01.txt  THP1run.hg19.AllInteractions.SP4.FDR0.1.txt

 - Header

RefSeqName      TranscriptName
Feature_Chr     Feature_Start   Annotation      Strand
InteractorName  InteractorID    Interactor_Chr  Interactor_Start        Interactor_End  InteractorAnnotation    distance
THP1.nLPS.rep1_SuppPairs        THP1.nLPS.rep3_SuppPairs        THP1.wLPS.rep1_SuppPairs        THP1.wLPS.rep3_SuppPairs
THP1.nLPS.rep1_p_value  THP1.nLPS.rep3_p_value  THP1.wLPS.rep1_p_value  THP1.wLPS.rep3_p_value # THP1.wLPS.rep3-> THP1 with LPS 2 replica 3.
IntGroup
Normal  LPSTreated      NofInts

- GM12878 datasets

GM12878.hg19.AllInteractions.SP4.FDR0.001.txt  GM12878.hg19.AllInteractions.SP4.FDR0.01.txt  GM12878.hg19.AllInteractions.SP4.FDR0.1.txt

- Header

RefSeqName      TranscriptName
Feature_Chr     Feature_Start   Annotation      Strand
InteractorName  InteractorID    Interactor_Chr  Interactor_Start        Interactor_End  InteractorAnnotation    distance
GM12878_rep1_SuppPairs  GM12878_rep2_SuppPairs  GM12878_rep3_SuppPairs  GM12878_rep1_p_value    GM12878_rep2_p_value GM12878_rep3_p_value
IntGroup
R1      R2      R3      NofReps # Replicate 3 is not sequenced as deep, many interactions are not expected to replicate in R3

- HEKa datasets

HEKa.hg19.AllInteractions.SP4.FDR0.1.txt, etc
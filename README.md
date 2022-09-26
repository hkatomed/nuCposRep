# nuCposRep
Codes for reproducing the results

##############################################################
# [1] FILE PREPARATION
# 
# Before running the R scripts, 
# appropriate files marked with (!) in the next section must be prepared.
# Also, some packages must be installed.
# Large data (~27 GB) will be stored in HOMEDIR.
#  
## These three files can be obtained from indicated URLs.  (bug fixed for the yeast genome)
#   Schizosaccharomyces_pombe.ASM294v2.dna.toplevel.fa.gz
#     URL: ftp://ftp.ensemblgenomes.org/pub/fungi/release-42/fasta/schizosaccharomyces_pombe/dna/
#   S288C_reference_sequence_R64-1-1_20110203.fsa
#     URL: http://sgd-archive.yeastgenome.org/sequence/S288C_reference/genome_releases/S288C_reference_genome_R64-1-1_20110203.tgz
#   Mus_musculus.NCBIM37.67.dna.toplevel.fa.gz
#     URL: ftp://ftp.ensembl.org/pub/release-67/fasta/mus_musculus/dna/
#   -> Store them in HOMEDIR/genome
#     
## These two files must be produced from the original files.
#   chr19_Chemical_NCPscore.txt
#   chr19_unique.map_95pc.txt
#   
#   The original files can be obtained from indicated URL.
#     GSM2183909_Chemical_NCPscore.txt.gz (3.6 Gb)
#     GSM2183909_unique.map_95pc.txt.gz   (35.4 Mb)
#       URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE82127
#   
#   Extraction of Chr19 data can be done as follows.
#     gunzip -dk GSM2183909_Chemical_NCPscore.txt.gz
#     gunzip -dk GSM2183909_unique.map_95pc.txt.gz
#     grep chr19 GSM2183909_Chemical_NCPscore.txt > chr19_Chemical_NCPscore.txt
#     awk 'BEGIN{FS=" ";OFS="\t"}{if($1=="19"){$1="chr19"; print ($0)}}' GSM2183909_unique.map_95pc.txt > chr19_unique.map_95pc.txt
#     -> Store them in HOMEDIR/mm_chemicalmap
# 
## The nucleosome data for mouse "chem_mm9_NucDNA.RData,"
#   "nature11142_s2_147.RData" and "sd01_147.RData"
#   can be generated in the process of nuCpos parameter construction,
#   of which scripts are available at https://doi.org/10.5281/zenodo.3362065. 
#   These three data are required for the dinucleotide frequency analysis.
#   -> Store them in HOMEDIR/RData
#
## The nucleosome data for yeasts, "nature11142_s2_sacCer3.RData,"
#   "nature11142_s3_sacCer3.RData", "sd01.RData" and "sd02.RData"
#   can be found at https://doi.org/10.5281/zenodo.3362065. 
#   -> Store them in HOMEDIR/RData
# 
# Parameters (sysdata.rda) of NuPoP and nuCpos can be obtained from indicated URLs.
#   NuPoP:  http://52.71.54.154/packages/release/bioc/html/https://doi.org/doi:10.18129/B9.bioc.NuPoP
#   nuCpos: http://52.71.54.154/packages/release/bioc/html/https://doi.org/doi:10.18129/B9.bioc.nuCpos
# 
#   They (sysdata.rda) will be found in the R directories of NuPoP and nuCpos.
#     Change the names as follows.
#       cd NuPoP/R
#       mv sysdata.rda sysdata_NuPoP.rda
#       cd nuCpos/R
#       mv sysdata.rda nuCpos_NuPoP.rda
#       -> Store them in HOMEDIR/RData
# 
## The Fortran subroutine for HBA calculation must be compiled as follows.
#    cd HOMEDIR
#    cd src
#    R CMD SHLIB HBA_3.f90
#    R CMD SHLIB HBA_3sp.f90
# 
## Chereji's +1 and -1 nucleosomes can be found at indicated URL.
#    13059_2018_1398_MOESM2_ESM.xlsx
#    URL: URL: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1398-0
#    -> Store it in HOMEDIR/Chereji
#
## To check if the above requirements are ready, see Requirements.png.
# 
## Hardware and R packages we used:
#    iMac Pro (2.3 GHz Intel Xeon W, 256 GB 2666 MHz DDR4, Mojave 10.14.5)
#    Biostrings (2.52.0)
#    openxlsx   (4.1.0.1)
#    nuCpos     (1.2.0)
#    NuPoP      (1.34.0)


##############################################################
# [2] WHAT WILL THE SCRIPTS DO?
# 
## sp_genome.R
# loads
#   genome/Schizosaccharomyces_pombe.ASM294v2.dna.toplevel.fa.gz(!)
# outpots
#   sp_genome/predNuPoP_sc/sp_genome_chr*.fasta
#   sp_genome/predNuPoP_sp/sp_genome_chr*.fasta
#   sp_genome/predNuCpos_sc/sp_genome_chr*.fasta
#   sp_genome/predNuCpos_sp/sp_genome_chr*.fasta
# 
## sc_genome.R
# loads
#   genome/S288C_reference_sequence_R64-1-1_20110203.fsa(!) (bug fixed)
# outputs
#   sc_genome_2011/predNuPoP_sc/sc_genome_2011_chr*.fasta
#   sc_genome_2011/predNuPoP_sp/sc_genome_2011_chr*.fasta
#   sc_genome_2011/predNuCpos_sc/sc_genome_2011_chr*.fasta
#   sc_genome_2011/predNuCpos_sp/sc_genome_2011_chr*.fasta
# 
## mm_genome.R
# loads
#   genome/Mus_musculus.NCBIM37.67.dna.toplevel.fa.gz(!)
# outputs
#   mm_genome/predNuPoP_mm/mm_genome_chr19.fasta
#   mm_genome/predNuCpos_mm/mm_genome_chr19.fasta
#   RData/mm9_unmasked.RData
# 
## yeasts_dHMM2wig.R
# loads
#   scripts/dHMM2wig.R
# runs dHMM2wig() to perform predNuPoP and predNuCpos at sc_genome_2011 and sp_genome
# outputs (* for Sc: 1-16, * for Sp: 1-3)
#   sc_genome_2011/predNuPoP_sc/NuPoPScSc_147bp_chr*.txt
#   sc_genome_2011/predNuPoP_sc/NuPoPScSc_147bp_chr*_Affinity.wig
#   sc_genome_2011/predNuPoP_sc/NuPoPScSc_147bp_chr*_Occupancy.wig
#   sc_genome_2011/predNuPoP_sc/NuPoPScSc_147bp_chr*_P_dyad.wig
#   sc_genome_2011/predNuPoP_sc/NuPoPScSc_147bp_chr*_P_start.wig
#   sc_genome_2011/predNuPoP_sc/NuPoPScSc_147bp_chr*_Viterbi.wig
#   sc_genome_2011/predNuPoP_sp/NuPoPScSp_147bp_chr*.txt
#   sc_genome_2011/predNuPoP_sp/NuPoPScSp_147bp_chr*_Affinity.wig
#   sc_genome_2011/predNuPoP_sp/NuPoPScSp_147bp_chr*_Occupancy.wig
#   sc_genome_2011/predNuPoP_sp/NuPoPScSp_147bp_chr*_P_dyad.wig
#   sc_genome_2011/predNuPoP_sp/NuPoPScSp_147bp_chr*_P_start.wig
#   sc_genome_2011/predNuPoP_sp/NuPoPScSp_147bp_chr*_Viterbi.wig
#   sc_genome_2011/predNuCpos_sc/nuCposScSc_147bp_chr*.txt
#   sc_genome_2011/predNuCpos_sc/nuCposScSc_147bp_chr*_Affinity.wig
#   sc_genome_2011/predNuCpos_sc/nuCposScSc_147bp_chr*_Occupancy.wig
#   sc_genome_2011/predNuCpos_sc/nuCposScSc_147bp_chr*_P_dyad.wig
#   sc_genome_2011/predNuCpos_sc/nuCposScSc_147bp_chr*_P_start.wig
#   sc_genome_2011/predNuCpos_sc/nuCposScSc_147bp_chr*_Viterbi.wig
#   sc_genome_2011/predNuCpos_sp/nuCposScSp_147bp_chr*.txt
#   sc_genome_2011/predNuCpos_sp/nuCposScSp_147bp_chr*_Affinity.wig
#   sc_genome_2011/predNuCpos_sp/nuCposScSp_147bp_chr*_Occupancy.wig
#   sc_genome_2011/predNuCpos_sp/nuCposScSp_147bp_chr*_P_dyad.wig
#   sc_genome_2011/predNuCpos_sp/nuCposScSp_147bp_chr*_P_start.wig
#   sc_genome_2011/predNuCpos_sp/nuCposScSp_147bp_chr*_Viterbi.wig
#   sp_genome/predNuPoP_sc/NuPoPSpSc_147bp_chr*.txt
#   sp_genome/predNuPoP_sc/NuPoPSpSc_147bp_chr*_Affinity.wig
#   sp_genome/predNuPoP_sc/NuPoPSpSc_147bp_chr*_Occupancy.wig
#   sp_genome/predNuPoP_sc/NuPoPSpSc_147bp_chr*_P_dyad.wig
#   sp_genome/predNuPoP_sc/NuPoPSpSc_147bp_chr*_P_start.wig
#   sp_genome/predNuPoP_sc/NuPoPSpSc_147bp_chr*_Viterbi.wig
#   sp_genome/predNuPoP_sp/NuPoPSpSp_147bp_chr*.txt
#   sp_genome/predNuPoP_sp/NuPoPSpSp_147bp_chr*_Affinity.wig
#   sp_genome/predNuPoP_sp/NuPoPSpSp_147bp_chr*_Occupancy.wig
#   sp_genome/predNuPoP_sp/NuPoPSpSp_147bp_chr*_P_dyad.wig
#   sp_genome/predNuPoP_sp/NuPoPSpSp_147bp_chr*_P_start.wig
#   sp_genome/predNuPoP_sp/NuPoPSpSp_147bp_chr*_Viterbi.wig
#   sp_genome/predNuCpos_sc/nuCposSpSc_147bp_chr*.txt
#   sp_genome/predNuCpos_sc/nuCposSpSc_147bp_chr*_Affinity.wig
#   sp_genome/predNuCpos_sc/nuCposSpSc_147bp_chr*_Occupancy.wig
#   sp_genome/predNuCpos_sc/nuCposSpSc_147bp_chr*_P_dyad.wig
#   sp_genome/predNuCpos_sc/nuCposSpSc_147bp_chr*_P_start.wig
#   sp_genome/predNuCpos_sc/nuCposSpSc_147bp_chr*_Viterbi.wig
#   sp_genome/predNuCpos_sp/nuCposSpSp_147bp_chr*.txt
#   sp_genome/predNuCpos_sp/nuCposSpSp_147bp_chr*_Affinity.wig
#   sp_genome/predNuCpos_sp/nuCposSpSp_147bp_chr*_Occupancy.wig
#   sp_genome/predNuCpos_sp/nuCposSpSp_147bp_chr*_P_dyad.wig
#   sp_genome/predNuCpos_sp/nuCposSpSp_147bp_chr*_P_start.wig
#   sp_genome/predNuCpos_sp/nuCposSpSp_147bp_chr*_Viterbi.wig
# 
## yeasts_wig2ROC.R
#   (*: NuPoPScSc, NuPoPScSp, nuCposScSc, nuCposScSp, 
#       NuPoPSpSc, NuPoPSpSp, nuCposSpSc, nuCposSpSp)
# loads
#   scripts/getPredTable.R
#   scripts/getNearestDyad.R
#   scripts/nearestDyadPred.R
#   scripts/getNearestDyadVit.R
#   scripts/NearestDyadVit.R
#   scripts/matchChem.R
#   scripts/MultiMatchChem.R
#   scripts/pred.summary.R
#   scripts/wig2ROC.R
#   scripts/sc_chemical.R
#     this script loads 
#       RData/nature11142_s2_sacCer3.RData (!)
#       RData/nature11142_s3_sacCer3.RData (!)
#   scripts/sp_chemical.R
#     this script loads 
#       RData/sd01.RData (!)
#       RData/sd02.RData (!)
# runs wig2ROC() at sc_genome_2011 and sp_genome
#   wig2ROC()
#     runs getPredTable() and generates *.allpred
#     runs MultiMatchChem() and revise *.allpred
#     runs pred.summary() and generate *.summary
#     runs getNearestDyad() and generate nature11142_s2.*.allpred, etc.
#     runs getNearestDyadVit() and generate nature11142_s2.*, etc.
# outputs 
#   RData/*_allpred.RData
#   RData/*_summary.RData
#   RData/*_perf_UniqueW0.RData
#   RData/*_perf_RedundantW0.RData
#   RData/nature11142_s2_*_allpred.RData
#   RData/sd01_*_allpred.RData
#   RData/nature11142_s3_*_allpred.RData
#   RData/sd02_*_allpred.RData
#   RData/nature11142_s2_*.RData
#   RData/sd01_*.RData
#   RData/nature11142_s3_*.RData
#   RData/sd02_*.RData
# 
## mouse_dHMM2wig.R
# loads
#   scripts/dHMM2wig.R
# runs dHMM2wig() to perform predNuPoP and predNuCpos at mm_genome
# outputs
#   mm_genome/predNuPoP_mm/NuPoPMmMm_147bp_chr19.txt
#   mm_genome/predNuPoP_mm/NuPoPMmMm_147bp_chr19_Affinity.wig
#   mm_genome/predNuPoP_mm/NuPoPMmMm_147bp_chr19_Occupancy.wig
#   mm_genome/predNuPoP_mm/NuPoPMmMm_147bp_chr19_P_dyad.wig
#   mm_genome/predNuPoP_mm/NuPoPMmMm_147bp_chr19_P_start.wig
#   mm_genome/predNuPoP_mm/NuPoPMmMm_147bp_chr19_Viterbi.wig
#   mm_genome/predNuCpos_mm/nuCposMmMm_147bp_chr19.txt
#   mm_genome/predNuCpos_mm/nuCposMmMm_147bp_chr19_Affinity.wig
#   mm_genome/predNuCpos_mm/nuCposMmMm_147bp_chr19_Occupancy.wig
#   mm_genome/predNuCpos_mm/nuCposMmMm_147bp_chr19_P_dyad.wig
#   mm_genome/predNuCpos_mm/nuCposMmMm_147bp_chr19_P_start.wig
#   mm_genome/predNuCpos_mm/nuCposMmMm_147bp_chr19_Viterbi.wig
# 
## mouse_ChemMap.R
# loads 
#   mm_chemicalmap/chr19_Chemical_NCPscore.txt (!)
#   mm_chemicalmap/chr19_unique.map_95pc.txt (!)
# outputs
#   RData/mm9_chr19_Red30.RData
#   RData/mm9_chr19_Red.RData
#   RData/mm9_chr19_Uni.RData
# 
## mouse_wig2ROC.R
#   (*: NuPoPMmMm, nuCposMmMm)
# loads
#   scripts/getPredTableMm.R
#   scripts/getNearestDyadMm.R
#   scripts/nearestDyadPredMm.R
#   scripts/getNearestDyadVitMm.R
#   scripts/NearestDyadVitMm.R
#   scripts/matchChemMm.R
#   scripts/MultiMatchChemMm.R
#   scripts/pred.summaryMm.R
#   scripts/wig2ROCMm.R
#   scripts/mm_chemical.R
#     this script loads 
#       RData/mm9_chr19_Uni.RData
#       RData/mm9_chr19_Red30.RData
# runs wig2ROCMm() at mm_genome
#   wig2ROCMm()
#     runs getPredTableMm() and generates *.allpred
#     runs MultiMatchChemMm() and revise *.allpred
#     runs pred.summaryMm() and generate *.summary
#     NOT runs getNearestDyadMm() and generate mm9_chr19_Uni.*.allpred, etc.
#     NOT runs getNearestDyadVitMm() and generate mm9_chr19_Red30.*, etc.
# outputs 
#   RData/*_allpred.RData
#   RData/*_summary.RData
#   RData/*_perf_UniqueW0.RData
#   RData/*_perf_RedundantW0.RData
#   NOT RData/mm9_chr19_Uni_*_allpred.RData
#   NOT RData/mm9_chr19_Red30_*_allpred.RData
#   NOT RData/mm9_chr19_Uni_*.RData
#   NOT RData/mm9_chr19_Red30_*.RData
# 
## yeast_genomes.R
# loads
#   genome/S288C_reference_sequence_R64-1-1_20110203.fsa (bug fixed)
#   genome/Schizosaccharomyces_pombe.ASM294v2.dna.toplevel.fa.gz
# outputs
#   RData/sc_genome_2011.RData
#   RData/sc_genome.RData
#   RData/sp_genome.RData
# 
## dHMMexamination.R
# processes data for Figures 1A, 1B, 1C, 7C, 7D, S2, S8 and Table S1
# outputs
#   RData/Rate_Redundant_yeasts.RData
#   RData/Rate_Unique_yeasts.RData
#   RData/AUC_summary_yeasts.RData
#   RData/num_summary_yeasts.RData
#   RData/nature11142_s3.NuPoPScSc.PredVsChem.plot400.RData
#   RData/nature11142_s3.nuCposScSc.PredVsChem.plot400.RData
#   RData/Rate_Redundant_mouse.RData
#   RData/Rate_Unique_mouse.RData
#   RData/AUC_summary_mouse.RData
#   RData/num_summary_mouse.RData
# 
## HBA_calc_inVitro.R
# processes data for Figures 2A, 2B and 2C
# outputs
#   RData/ori601HBA_m_mean.RData
#   RData/ori601HBA_c_mean.RData
#   RData/MMTVHBA_m_mean.RData
#   RData/MMTVHBA_c_mean.RData
#   RData/MMTV.RData
#   RData/MMTV_DYAD.RData
#   RData/rDNAHBA_m_mean.RData
#   RData/rDNAHBA_c_mean.RData
#   RData/rDNA.RData
#   RData/rDNA_DYAD.RData
# 
## HBA_calc_inVivo.R
# processes data for Figures 3A and 3B
# outputs
#   RData/Chereji_plus1_HBAm.RData
#   RData/Chereji_plus1_HBAc.RData
#   RData/Chereji_minus1_HBAm.RData
#   RData/Chereji_minus1_HBAc.RData
#   RData/Chereji_plus1_AT.RData
#   RData/Chereji_minus1_AT.RData
# 
## dHMM_locus.R
# processes data for Figures 3C and S3
# outputs
#   RData/results_nuCpos_Sc_TRP1ARS1.RData
#   RData/results_NuPoP_Sc_TRP1ARS1.RData
#   RData/TRP1ARS1_AT.RData
#   RData/results_NuPoP_ura4.RData
#   RData/results_nuCpos_ura4.RData
#   RData/ura4_dyads.RData
#   RData/seqUra4_AT.RData
#   
## dHMM_mutants.R
# processes data for Figures 4 and S4
# outputs
#   RData/mutNuCpos_BAR1.RData
#   RData/mutNuCpos_TALS.RData
#   
## lHBA_inVitro.R
# processes data for Figures 5 and S6
# outputs
#   RData/Widom_lHBA.RData
#   RData/lHBA_Widom601.RData
#   RData/lHBA_NCP147.RData
#   RData/lHBA_Widom601LR.RData
#   RData/lHBA_Widom601L.RData
#   RData/lHBA_Widom601R.RData
#   
## lHBA_inVivo.R
# processes data for Figures 6A and 6B
# outputs
#   RData/ura4_seqs.RData
#   RData/lHBA_ura4_WT.RData
#   RData/lHBA_ura4_Dyad.RData
#   RData/lHBA_ura4_Linker.RData
#   RData/lHBA_ura4_Int.RData
# 
## dinucFreq.R
# processes data for Figures 7A, 7B and S7
# outputs
#   RData/dinuc.RData
#   RData/dinucNorm.RData
#   RData/dinucNorm2.RData
# 
## Figure.R
# draws all the figures in the paper.
# outputs
#   RData/Rplots.pdf


##############################################################
# [3] RUNNING THE SCRIPTS FOR DATA PROCESSING
# 
## Run the scripts below to perform prediction and produce data for examination.
# cd HOMEDIR
# Rscript scripts/sp_genome.R
# Rscript scripts/sc_genome.R
# Rscript scripts/mm_genome.R		# takes long time
# Rscript scripts/yeasts_dHMM2wig.R	# takes long time
# Rscript scripts/yeasts_wig2ROC.R	# takes long time
# Rscript scripts/mouse_dHMM2wig.R	# takes long time
# Rscript scripts/mouse_ChemMap.R
# Rscript scripts/mouse_wig2ROC.R	# takes long time
# Rscript scripts/yeast_genomes.R
# Rscript scripts/dHMMexamination.R
# Rscript scripts/HBA_calc_inVitro.R
# Rscript scripts/HBA_calc_inVivo.R	# takes long time
# Rscript scripts/dHMM_locus.R
# Rscript scripts/dHMM_mutants.R
# Rscript scripts/lHBA_inVitro.R
# Rscript scripts/lHBA_inVivo.R
# Rscript scripts/dinucFreq.R
# Rscript scripts/Figures.R		# outputs all the figures

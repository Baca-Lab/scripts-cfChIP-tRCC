#############
README to remake the analysis from the article tRCC_cf-ChIP
#############

R packages required are described in the scripts
R version used was 4.2.0
Version of other softawre used are present in the manuscript methods section

### DiffBind_RCC.R (expected time 2 hours)
### Script to perform the differential peak enrichment using DiffBind R package (v3.10.1), between ccRCC and tRCC cell lines, to identify tRCC-up peaks 
### Example set with the H3K4me3 data from RCC cell lines
- Download FASTQ files and narrowpeaks from cell lines at GSE280708
- Perform bam alignment and peak calling with ChiLin pipeline (chilin simple -p narrow -r histone)
- Set paths to bam and narrowpeaks in the manifestDiffbindRCC.csv
- Run the script



### cfChIP.R
### Script to process cf-ChIP samples, and estimate the coverage at sites identified in cell lines normalised to the 10k highest common peaks
- Download FASTQ files and narrowpeaks from plasma samples at GSE280708, GSE266530 (expected time 24 hours)
- Generate the coverage RDS files for each plasma sample (using convert_BAM_to_RDS/bam_to_frag.R script - expected time 24 hours)
- Fill the manifest with (data_plasmas_RCC.csv) with path to coverage RDS files
- Turn the reports from DiffBind into bed files in hg19, to use as sites. top_common and blacklist are provided
- Run the script (expected time 3 hours)




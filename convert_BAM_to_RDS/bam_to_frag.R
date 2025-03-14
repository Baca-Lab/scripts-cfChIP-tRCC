# convert bam files to bed files and granges objects with fragment locations

library('stringr')
library('rtracklayer')

bam_dir='Bam'
frag_dir='frags'
rds_dir='rds'
dir.create(frag_dir)
dir.create(rds_dir)

bam_list = list.files(bam_dir, full.names = TRUE)

count = 0
for(bam_file in bam_list) {
  count = count + 1
  bed_file = str_replace(bam_file, bam_dir, frag_dir) %>%
    str_replace('.bam', '_fragments.bed')
  rds_file = str_replace(bam_file, bam_dir, rds_dir) %>%
    str_replace('.bam', '.RDS')
  message('processing bam file ', count, ' of ', length(bam_list))
  if(file.exists(paste0(bed_file, '.gz'))) {
    message(bed_file, ' already exists, skipping')
  } else {
    message('extracting fragments from ', bam_file)
    cmd = paste0("module load samtools bedtools\n 
               samtools view  -b -f1 ", bam_file, " | 
               samtools sort -n | bedtools bamtobed -bedpe 2> /dev/null | 
               awk 'BEGIN{{OFS=\"\t\";FS=\"\t\"}}($1==$4){{print $1, $2, $6}}' > ", bed_file,
                 "\n gzip ", bed_file)
 
    #samtools sort -n | (to use before bedtools if files not sorted)
       
    #for single end bam files use the cmd below instead:
#    cmd = paste0("module load samtools bedtools\n 
#               samtools view  -b ", bam_file, " | 
#               bedtools bamtobed -i - > ", bed_file,
#                 "\n gzip ", bed_file)
    
    
    system(cmd)
  }
  
  if(file.exists(rds_file)) { #file.exists(rds_file)
    message(rds_file, ' already exists, skipping')
  } else {
    message('making granges object for ', bam_file)
    tmp = import(paste0(bed_file, '.gz'))
    tmp$frag_width = width(tmp)
    tmp = resize(tmp, width = 1, fix='center')
    saveRDS(tmp, file = rds_file)
  }
}


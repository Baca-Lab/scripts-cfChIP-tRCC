
library(stringr)
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(dplyr)

# Get counts from cfChIP granges objects at a set of sites.
# Given a character vector with the names of rds files containing granges objects with fragment locations
# and a granges object (or bed file) with targets, get normalized fragment counts at the targets.
frag_counts = function(frag_files,
                       sample_names,
                       targets, # bed file or granges object with target region coordinates
                       normalize = 'total_frags',  # 'total_frags' or a granges object or bed file name with sites for normalization
                       verbose = TRUE) {  
  
  
#  stopifnot(class(normalize) == 'GRanges' | normalize == 'total_frags' | grepl('bed$', normalize))
  
  if(class(targets) == 'character') {
    targets = import(targets)
  }
  
  if(class(normalize) == 'character' && grepl('bed$', normalize)) {
    normalize = import(normalize)
  }
  
  target_names = paste0(seqnames(targets), ':', start(targets), '-', end(targets))
  counts_mat = matrix(NA, length(targets), length(frag_files))
  rownames(counts_mat) = target_names
  colnames(counts_mat) = sample_names
  
  count = 0
  for (f in frag_files) {
    count = count + 1
    
    if(verbose) {
      message('processing sample ', count, ' of ', length(frag_files))
      message('reading in ', f)
    }
    
    frag = readRDS(f)
    # ensure frags have been reduced to 1 bp
    if(!(all(width(frag) == 1))) {
      frag = resize(frag, width = 1, fix = 'center')
    }
    
    if (class(normalize) == 'character' && normalize == 'total_frags') {
      norm_denominator = length(frag)
      
    } else {
      norm_denominator = sum(countOverlaps(normalize, frag)) %>% suppressWarnings()
    }    
    counts_mat[, count] = countOverlaps(targets, frag) / norm_denominator * 1e6
    
  }
  counts_mat
}

# write regions of the form "chrN:start-end" to a bed file
write_regions_to_bed = function(regions, out_file) {
  regions %>% str_split('[:-]', simplify = TRUE) %>%
    write.table(
      file = out_file,
      quote = F,
      row.names = F,
      col.names = F,
      sep = '\t'
    )
}



# function to tile across sites
tile_sites = function(sites_file_list,
                      bin_width = 40,
                      flank_bp = 3000,
                      exclude_regions = NA) {
  # regions to be excluded - eg, blacklisted sites, or peaks in healthy cfChIP data
  # filter sites and get sites names
  return_sites_list = list()
  for (sites_file in sites_file_list) {
    sites = import(sites_file)
    if (class(exclude_regions) == 'GRanges') {
      sites = sites[!(sites %over% exclude_regions)]
    }
    
    #make sure ROIs are uniform size and expand to include flanks
    sites = resize(sites, width = 1, fix = 'center')
    sites = sites + flank_bp
    sites = resize(sites, flank_bp * 2, fix = 'start')
    
    # remove any sites that overlap with each other
    sites = GenomicRanges::reduce(sites)
    sites = subset(sites, width == flank_bp * 2)
    
    # create tiles
    tiles = tile(sites, width = bin_width) %>% unlist()
    
    # map eachtile back to the site it came from
    ol = findOverlaps(tiles, sites)
    agg = tiles[from(ol)]
    tmp = sites[to(ol)]
    agg$site_start = start(tmp)
    
    agg$bin = start(agg) - agg$site_start - flank_bp
    return_sites_list = c(return_sites_list, agg)
  }
  return_sites_list
}


# function to get signal at sites
signal_at_sites = function(frags_file,
                           sites_list,
                           sites_names,
                           remove_top = 0.05,
                           remove_peaks_wider_than = NA) {
  
  frags = readRDS(frags_file)
  
  # add back fragment end and start based on width if needed
  if (all(width(frags)) == 1 && !is.null(frags$frag_width)) {
    frags = resize(frags, width = frags$frag_width, fix = 'center')
  }
  
  return_df = NULL
  
  count = 1
  for (sites in sites_list) {

    if(!is.na(remove_peaks_wider_than)) {
      sites = sites[width(sites) <= remove_peaks_wider_than]
    }
    
    sites$counts = suppressWarnings(countOverlaps(sites, frags)) 
    
    # remove top quantile of sites by fragment count
    tmp = sites %>% as.data.frame() %>% group_by(site_start) %>% summarize(tot = sum(counts))
    quantile_cutoff = quantile(tmp$tot, 1 - remove_top)
    sites_to_remove = tmp$site_start[tmp$tot > quantile_cutoff]
    
    sites = subset(sites, !(site_start %in% sites_to_remove))
    counts = sites %>%
      as.data.frame() %>%
      group_by(bin) %>%
      summarize(
        sites = sites_names[count],
        frag_counts = sum(counts),
      )
    counts$total_counts = length(frags)
    return_df = rbind(return_df, counts)
    count = count + 1
    
  }
  return_df
}


# function to plot profiles output from signal_at_sites
plot_signal_at_sites = function(toplot,
                                normalize_to_these_sites = 'none', # normalize to signal at these sites - must be included in toplot$sites. set to none for normalization to total counts
                                subtract_shoulder = TRUE,
                                group_by_this = NA,
                                main_group = NA, #which group should be colored in the plots (if not specified, it picks the group based on highest signal)
                                auc_boxplot = FALSE, # generate a boxplot of area under the curve, rather than showing the signal curves themselves?
# not yet supported:             min_normalization_counts = NA, # require at least this many fragments in the normalization sites for the sample to be included
                                out_file = NA # tsv file to write scores to, if desired. Only used if auc_boxplot is TRUE
                                ) {
  
  stopifnot(group_by_this %in% colnames(toplot))
  
  
  if(subtract_shoulder) {
    # normalize counts to shoulders
    shoulders = 2800

    
    shoulder_norm = toplot %>% 
      subset(bin < -shoulders | bin > shoulders) %>% 
      group_by(study_name, sites) %>% summarize(shoulder_norm_factor = median(frag_counts))
    
    toplot = merge(toplot, shoulder_norm)
    toplot$frag_counts = toplot$frag_counts - toplot$shoulder_norm_factor
    toplot$frag_counts[toplot$frag_counts < 0] = 0
  }
  
  if(all(normalize_to_these_sites == 'none')) {
    toplot$counts_norm = toplot$frag_counts / toplot$total_counts
    
  } else{

    # allows normalizing to multiple sets of sites: (TODO: remove)
    toplot$counts_norm = toplot$frag_counts
    for(these_sites in normalize_to_these_sites) {
      # drop the insample_norm_factor if it already exists
      if('insample_norm_factor' %in% colnames(toplot))
        toplot = toplot %>% dplyr::select(-insample_norm_factor)
      insample_norm = toplot %>% subset(sites == these_sites) %>%
        group_by(study_name) %>%
#        summarize(insample_norm_factor = mean(frag_counts))   
        summarize(insample_norm_factor = mean(counts_norm))  
      toplot = merge(toplot, insample_norm)
 #     toplot$counts_norm = toplot$frag_counts / toplot$insample_norm_factor
      toplot$counts_norm = toplot$counts_norm / toplot$insample_norm_factor  
    }
  }
  
  
    
  #  
  
#  # get an area under the curve score for each sample
#  auc_score = toplot %>% group_by(study_name, sites) %>% 
#    summarize(auc_score = sum(counts_norm),
#              auc_score_insample_norm = sum(counts_insample_norm))
  
#  toplot = merge(toplot, auc_score)
  
  
#  if(FALSE) { # remove samples with low in sample norm factors (likely due to low)
#    insample_norm_factor_cutoff = 2^(mean(log(toplot$insample_norm_factor,2)) - 
#                                       sd(log(toplot$insample_norm_factor,2)) * 2)
    
#    #    insample_norm_factor_cutoff = 100
#    toplot = subset(toplot, insample_norm_factor > insample_norm_factor_cutoff)
    
#    #  insample_norm_cut = unique(toplot$insample_norm_factor) %>% quantile(probs = 0.1)
#    #  toplot = subset(toplot, insample_norm_factor > insample_norm_cut)
    
#  } 
  
  
  #group_by = 'additional_annotation'
  

  if(is.na(main_group)) {
  # assign colors for plots. Color the subtype with the max value at on-target sites:

  target_sites_name = unique(toplot$sites)[unique(toplot$sites) != normalize_to_these_sites][1]
  
  tmp = subset(toplot, sites == target_sites_name) %>% 
    group_by(study_name) %>% 
    summarize(counts_agg = sum(counts_norm)) %>%
    slice_max(order_by = counts_agg, n = 1)
  type1 = toplot[toplot$study_name == tmp$study_name, group_by_this] %>% unique()
  type2 = unique(toplot[,group_by_this])[unique(toplot[,group_by_this]) != type1]
  type3 = unique(toplot[,group_by_this])[unique(toplot[,group_by_this]) != c(type2, type1)]
  } else {
    type1 = main_group
    type2 = unique(toplot[,group_by_this])[unique(toplot[,group_by_this]) != type1]
    type3 = unique(toplot[,group_by_this])[unique(toplot[,group_by_this]) != c(type2, type1)]
  }
  toplot[,group_by_this] = factor(toplot[,group_by_this], levels = c(type1, type2, type3))  
  cols = setNames(c('#198CCE', 'darkorange', 'darkgreen'), c(type1, type2, type3))
  
  if(!auc_boxplot) {
  

    p = ggplot(toplot, aes_string(y = 'counts_norm', group = group_by_this, x = 'bin', col = group_by_this)) +
      theme_classic() +

      stat_smooth(method="loess", span=0.1, se=TRUE) +
#      stat_summary(fun = mean, geom = 'line', size = 1) +
 #KS     geom_line(aes_string(y = 'counts_norm', group = paste0('interaction(', group_by_this,', study_name)'), x = 'bin', col = group_by_this), 
 #KS              stat="smooth", span = 0.1, alpha = 0.4) +
      geom_vline(xintercept = 0, linetype = 'dashed') +
      scale_color_manual(values = cols) +
      facet_wrap(~sites, scales = 'free')
    
    
    
# if(FALSE){
#     ###tmp testing!
#     mean_signal = toplot %>% group_by(sites, bin, cancer_subtype) %>% summarize(agg_counts = sum(counts_norm))
#     p = ggplot(mean_signal, aes_string(y = 'agg_counts', group = group_by_this, x = 'bin', col = group_by_this)) +
#       theme_classic() +
#       #      stat_smooth(method="loess", span=0.05, se=TRUE, aes_string(fill=group_by_this), alpha=0.2) +
#       stat_smooth(method="loess", span=0.05, se=FALSE) +
# #      geom_line(aes_string(y = 'counts_norm', group = paste0('interaction(', group_by_this,', study_name)'), x = 'bin', col = group_by_this), 
# #                stat="smooth", span = 0.05, alpha = 0.2) +
#       geom_vline(xintercept = 0, linetype = 'dashed') +
#       facet_wrap(~sites, scales = 'free')
# }    
#     
      
    ###
    
#    p = ggplot(toplot, aes_string(y = 'counts_norm', group = paste0('interaction(', group_by_this,', study_name)'), x = 'bin', col = group_by_this)) +
#      theme_classic() +
#      geom_line(stat="smooth", span = 0.05, alpha = 0.5) +
#      geom_vline(xintercept = 0, linetype = 'dashed') +
#      facet_wrap(~sites, scales = 'free')

    
  } else {
#    toplot_orig = toplot
    auc_score = toplot %>% 
      group_by(study_name, sites) %>%
#      subset(bin < -80 & bin > -1500 | bin > 80 & bin < 1500) %>% #tmp - may change this
      summarize(auc = sum(counts_norm))

    toplot = merge(toplot, auc_score) %>% 
    #  dplyr::select(-c(bin,frag_counts)) %>%
      subset(bin == 0) %>%
      distinct()

    # get p-vals for each comparison
    grps = unique(toplot[,group_by_this])
    stopifnot(length(grps) == 3)

    sites_labels = NULL
    group1_ids = toplot$study_name[toplot[,group_by_this] == type1] %>% unique()
    group2_ids = toplot$study_name[toplot[,group_by_this] == type2] %>% unique()

    for(sites_to_plot in unique(auc_score$sites)) {

#if(FALSE){ # test at the sample level, rather than peak level
        x = auc_score$auc[auc_score$study_name %in% group1_ids & auc_score$sites == sites_to_plot]     
        y = auc_score$auc[auc_score$study_name %in% group2_ids & auc_score$sites == sites_to_plot]
        tmp = wilcox.test(x, y)   


#} else { # test at the peak level

#      x = toplot_orig %>% group_by(sites, study_name) %>% 
#        subset(study_name %in% group1_ids & sites == sites_to_plot)
#      y = toplot_orig %>% group_by(sites, study_name) %>% 
#        subset(study_name %in% group2_ids & sites == sites_to_plot)
#      tmp = wilcox.test(x$counts_norm, y$counts_norm)              
#      browser()
#}
      


#      tmp = ks.test(x, y)       
        #TODO: fix this - sometimes the site labels are swapped (need to be assigned based on type1/type2)
      assign('sites', sites_to_plot)
      toadd = c(paste0(sites, '\np = ', signif(tmp$p.value,3)))
    sites_labels = c(sites_labels, toadd)  
      
    }
    names(sites_labels) = unique(auc_score$sites)

    # plot on log scale if there are major outliers
#    if(max(toplot$auc) > (median(toplot$auc) + 10 * min(toplot$auc))) {
      p = ggplot(toplot, aes_string(y = 'log(auc + 1, 2)', x = group_by_this, col = group_by_this, label = 'study_name'))
#    } else {
#      p = ggplot(toplot, aes_string(y = 'auc', x = group_by_this, col = group_by_this, label = 'study_name'))      
#    }
    p = p +
      theme_classic() +
      geom_boxplot(width = 0.8, outlier.shape = NA) +
      geom_jitter(width=0.3, aes(size = ctDNA_estimate)) + #aes(size = ctDNA_estimate)
      #geom_label(size = 2) +
      theme_classic() + 
      scale_color_manual(values = cols) +
      facet_wrap(~sites, scales = 'free', labeller = labeller(sites = sites_labels))
    
    if(!is.na(out_file)) {
      write.table(toplot, file = out_file, quote = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE)
    }

            
  }
  p
}
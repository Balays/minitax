
wrap.fun.complete <- function(taxa_or_bamfile, pattern='.bam', sample=NULL,
                              input='bamfile', ## possible inputs: c('DT', 'list', 'bamfile', '.rds'), 
                              #DT_or_List=NULL,
                              config=config, 
                              steps=c('1.Import_Alns', '2.Find_taxidentity', '3.Refine', '4.Summarise'),
                              keep.max.cigar=T, CIGAR_points=NULL, 
                              methods=NULL, 
                              saveRDS=F, cropnumbers=F,
                              outputs=unlist(strsplit(config$outputs, ', ')),
                              best.mapq=T, mapq.filt=NA,
                              db='proGcontigs', db.data=prog.db,  db.uni.data=prog.db.uni,
                              ranks = c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
                              
                              ) {
  bam <- data.table(NULL)
  bam.filt <- data.table(NULL)
  minimap2 <- data.table(NULL)
  minimap2.filt <- data.table(NULL)
  taxa <- data.table(NULL)
  taxident.all <- data.table(NULL)
  taxident.all.sum <- data.table(NULL)
  
  #### Inputs ####  --->>> This part is not final
          if (input == 'bamfile') {
    sample   <- gsub(pattern, '', gsub('.*/', '', taxa_or_bamfile))
    bamfile  <- taxa_or_bamfile
    ##
  } else if (input == 'taxa.DT') {
    taxa  <- fread(taxa_or_bamfile, na.strings = '')
    taxa[taxa == ''] <- NA
    sample <- gsub(paste0('_best_alignments_w_taxa.tsv'), '', taxa_or_bamfile)
    sample <- gsub('.*\\/', '', sample)
    #taxa$sample <- sample
    steps <- c('4.Summarise')
  }
  
  message('Starting minitax WF on: ', sample, ', at: ', Sys.time())
  message('The following steps will be performed: ', paste(steps, collapse = '; '))
  
  #### methods
  ## avilable methods: 'LCA', 'BestAln', 'SpeciesEstimate', 'RandAln'
  
  ## Calculate and keep alignments with highest CIGAR scores?
  if (keep.max.cigar) {
    calc.cigar.score <- T } else { calc.cigar.score <- F }
  
  ## try start
  try({ 
  
  #### 1. Import and possibly filter alignments based on flags and MAPQ ####
  #### wrap.fun0 ####
  if(any(steps %in% '1.Import_Alns')) {
    
    message('Importing: \n', bamfile, '...')
    
    #system.time({
    bam <- ov.from.bam2(
      bamfile =  bamfile,
      #bamnames,
      #write.filtered = F, force.create=force.create, filtering=F
      rm.gaps.in.aln = F, mapq.filt = NA,
      flag.tokeep    = NA,    flag.tocrop = NA, seqnames.tofilt = NA,
      is.lortia = F, add.primes = F
    )
    #})
    bam$sample <- sample
    
    setDT(bam)
    
    ## Filtering alignments for MAPQ
    if(!all(is.na(mapq.filt))) {
      message('Removing alignments with MAPQ of : ', paste0(mapq.filt, collapse = ';'))
      bam.filt <- bam[!is.element(mapq, mapq.filt), ]
    } else { bam.filt <- bam }
    
    ## Selecting the best alignment for multi-alignment reads based on MAPQ
    if(best.mapq) {
      message('Selecting the best alignment for multi-alignment reads based on MAPQ ...')
    } else {
      message('No filtering of alignments for multi-alignment reads based on MAPQ ...')
    }
    minimap2 <- get.best.aln(bam.filt, keep.chim = T, 
                                  rm.supplementary = F, 
                                  rm.secondary = F, 
                                  best.mapq = best.mapq, 
                                  give.qname.ID.to.secondary=F)
    setDT(minimap2)
    
    message('The number of filtered alignments is: ',
            nrow(bam.filt) - nrow(minimap2), ' (',
            round(nrow(minimap2) / nrow(bam.filt), 6) * 100, '%)')
    
    
    ## summarise
    message('Summarising alignments statistics ...')
    bam.dist <- bam.filt %>% distinct(sample, qname, mapq) # rname
    bam.dist$mapq[is.na(bam.dist$mapq)] <- 'NA-CIGAR'
    bam.dist <- bam.dist %>% group_by(sample, qname) %>% reframe(mapq, n_unique_aln=seq_along(mapq))
    bam.dist.sp <- bam.dist %>% spread(n_unique_aln, mapq)
    bam.dist.sp$is.all.NA <- apply(bam.dist.sp %>% dplyr::select(-c(qname)),            1, function(x) all(is.na(x) | x == 'NA-CIGAR')) 
    bam.dist.sp$max.mapq  <- apply(bam.dist.sp %>% dplyr::select(-c(qname, is.all.NA)), 1, function(x) max(x, na.rm = T)) 
    
    bam.sum  <- bam.dist.sp %>% group_by(sample, is.all.NA) %>% reframe(is.mapped=n())
    #bam.sum$sample <- sample
    
    ### write output
    if (any(outputs %in% 'alignments_w_taxa')) { 
      fwrite(bam.filt, paste0(outdir, '/', 'alignments_w_taxa/', sample,  '_alignments_w_taxa.tsv'), sep = '\t') }
    
    if (any(outputs %in% 'bam.sum')) { 
      fwrite(bam.sum, paste0(outdir, '/', 'bam.sum/', sample,  '_bamstats.tsv'), sep = '\t') }
   
  } else {
    #minimap2 <- DT
  }
  ####  ####
  
  #### 2. Refine alignments based on CIGAR scores ####
  #if (!exists('taxa')) {taxa <- DT}
  
  if(any(steps %in% '2.Refine')) {
    message('\n Refining alignments based on CIGAR scores... \n')
    
    ## CIGAR SCORES
    if (calc.cigar.score) {
      
      if(is.null(CIGAR_points)) {
        CIGAR_points <- data.frame( 
                 match_score           =  1, 
                 mismatch_score        = -3,
                 insertion_score       = -2, 
                 deletion_score        = -2, 
                 gap_opening_penalty   = -2, 
                 gap_extension_penalty = -1)
      }
      
      message('Calculating CIGAR scores...')
      
      # Apply the function using future
      ## -->> data.table is much faster!
      # minimap2[, cigar_score := unlist(future_lapply(cigar, cigar_score, ...))]
      
      # Apply the function using data.table
      #system.time({
      minimap2[, cigar_score := cigar_score(cigar, 
                                            match_score = CIGAR_points$match_score, 
                                            mismatch_score = CIGAR_points$mismatch_score,
                                            insertion_score = CIGAR_points$insertion_score, 
                                            deletion_score = CIGAR_points$deletion_score, 
                                            gap_opening_penalty = CIGAR_points$gap_opening_penalty, 
                                            gap_extension_penalty = CIGAR_points$gap_extension_penalty)]
      #})
      
      # Add a new column 'max_cigar_score'
      columns_to_group_by <- c('sample', 'qname')
      minimap2[, max_cigar := (cigar_score == max(cigar_score)), by = columns_to_group_by]
    }
   

    ## Select the rows with the highest 'cigar_score' for each 'qname'
    if (keep.max.cigar) {
      message('Keeping alignments with the highest CIGAR score for each read...')
      
      message('The number of filtered alignments is: ',
              nrow(minimap2[max_cigar == F, ]), ' (',
              round(nrow(minimap2[max_cigar == F, ]) / nrow(minimap2[, ]), 6) * 100, '%)')
      
      minimap2.filt <- minimap2[max_cigar == T, ]
    }
    
    print(head(minimap2.filt))
    
  } else {
    try({ minimap2.filt <- minimap2 })
    #if (exists('DT')) {taxa <- DT}
  }
    
  ####  ####
    
    
  #### 3. Find tax identity (minitax) ####
  ## running minitax -->> finding tax.identity and tax.identity level for each read
  if(any(steps %in% '3.Find_taxidentity')) {
    message('Finding tax.identity and tax.identity level for each read...')
    
    # Identify the last entry of the ranks vector
    #last_rank <- tail(ranks, n = 1)
    
    try({
    taxa <- minitax2(minimap2.filt, db=db, db.uni.data=db.uni.data, 
                     ranks = ranks )
    })
       
    taxa$sample <- sample
    taxa[,aln_nr := seq_along(cigar), by=.(sample, qname)]
    
    ## Where the tax.identity is below the lowest rank (superkingdom)
    taxa[is.na(tax.identity),       tax.identity       := 'unclassified'] 
    taxa[is.na(tax.identity.level), tax.identity.level := 'unclassified'] 
    
    ### write output
    if (any(outputs %in% 'best_alignments_w_taxa')) { 
      fwrite(taxa, paste0(outdir, '/', 'best_alignments_w_taxa/', sample,  '_best_alignments_w_taxa.tsv')) }

  } else {
    ##
  }
   
  ####  ####
  
    
  
  #### 4.Refine hits and Summarise  #### 
  if(any(steps %in% '4.Summarise')) {
    message('Summarising taxa on the follwing taxon levels: ', paste0(ranks, collapse='; '))
    
   
    ##
    try({
      
      if(methods == 'BestAln' | methods == 'RandAln') {
      
      dupqnames  <- unique(dup(taxa$qname))
      message('Trying to find most likely taxa for this number of reads: \n', 
              '\t', length(dupqnames), ' (', round((length(dupqnames) / length(unique(taxa$qname))) * 100, 3), '% of the reads) \n',
              '- as in these cases, both the MAPQ and the CIGAR scores were equal')
      
      taxa.filt <- taxa[!qname %in% dupqnames ]
      
      taxa.dup  <- taxa[ qname %in% dupqnames ]
      
      ##
      print(list('Distribution of tax.identity.levels (that is the taxon level, at which every aligment of a read share the same taxon):', 
                 taxa[,.N,tax.identity.level]))
    
      
      #### BestAln method  #### 
      ## 
      
      if (methods == 'BestAln') {
        message('Method: BestAln -->> Selecting the most probable taxa,
        by calculating the probability (ratio of alignments) for each taxon at each taxon level,
        and selecting that if it is above acertain thrshold and is better than the assigned tax.identity,
        for each read, in those cases where both the CIGAR and the MAPQ was the same...')
        ### NEW CODE
        # try to find most likely taxa - that is not a random alignment should be selected, but a probability assigned based on taxonomic lineage
        # and the most probable taxa be selected
        # this can be done with or without prior filtering for best alignment based on CIGAR and/or MAPQ as well
        ###
     
        ### THRESSHOLDS !!!
        thresholds  <- NULL; try({ 
          thresholds  <- config$BestAln_thresholds 
          stopifnot(length(thresholds) != length(ranks) | length(thresholds) != 1)
          thresholds  <- data.frame(rank=ranks, threshold = thresholds)
        })
        if (is.null(thresholds)) { thresholds  <- data.frame(rank=ranks, threshold = 0.6) }
        
        print(list('The used threshold to update the tax.identity for each level was: ', thresholds))
        
        taxa.nodup     <- BestAln(taxa.dup, ranks, thresholds=thresholds)
        
        
        ## if there are still duplicates remove them, randomly
        taxa.dup  <- unique(taxa.nodup[aln_nr > 1, qname])
        
        #taxa.nodup[qname == 'M06102:89:000000000-DDBG4:1:1101:13184:14802']
        taxa.nodup  <- taxa.nodup[aln_nr == 1, ]
        
      #### #####  
        
      ####  Random Alignment   #### 
      ##
      } else if (methods == 'RandAln') {
        
        message('Method: RandAln -->> Selecting a random alignment 
        for each read, in those cases where both the CIGAR and the MAPQ was the same...')
        taxa.nodup  <- taxa.dup[aln_nr == 1, ]
        
        ## overwriting tax identity to the randomly selected alignment's taxon
        taxa.nodup[,tax.identity       := NULL]
        taxa.nodup[,tax.identity.level := NULL]
        
        add_rightmost_non_na(taxa.nodup, ranks, 'tax.identity')
        add_rightmost_non_na_colname(taxa.nodup, ranks, 'tax.identity.level')
        
        
      }
      #### #####
      
      
      ## Add taxonomic lineage from database based on tax.identity and tax.identity level
      cols_by        <- c('qname', 'aln_nr', 'sample', ranks, 'tax.identity', 'tax.identity.level')
      taxa.nodupuni  <- unique(taxa.nodup[, ..cols_by] )  
      taxa.nodupuni  <- add_lineage(data = taxa.nodupuni, ranks, db.uni.data)
      
      ## Merge back the updated tax.identity with the read info
      cols_by        <- colnames(taxa)[!colnames(taxa) %in% colnames(taxa.nodupuni)]
      cols_by        <- c(cols_by, 'qname', 'aln_nr', 'sample')
      taxa.nodupuni  <- merge(taxa[,..cols_by], taxa.nodupuni[,], by=c('qname', 'aln_nr', 'sample'), all.y=T)
      
      cols_by        <- colnames(taxa.dup)
      taxa.nodupuni  <- unique(taxa.nodupuni[,..cols_by] )  
      
      ##
      taxa.nodup     <- taxa.nodupuni
      
      ## merge with non-multi aligned reads
      taxident.all  <- rbind(taxa.filt, taxa.nodup)
      
      ## Dereplicate reads, without alignments, keeping tax.identity and ranks only
      cols_by       <- c('qname', 'sample', ranks, 'tax.identity', 'tax.identity.level')
      taxident.all  <- unique(taxident.all[,..cols_by] )  
      taxident.all  <- unique(taxident.all)
      
      ##
      print(list('Distribution of tax.identity.levels, after refining:', 
                 taxident.all[,.N,tax.identity.level]))
    
    
      ## Adding lineage to differentiate taxa with the same rank name
      cols_by        <- c(ranks, 'tax.identity', 'tax.identity.level')
      taxident.uni   <- unique(taxident.all[, ..cols_by])
      paste_columns_dt(taxident.uni, ranks, '::', 'lineage', 1)
      
      cols_by        <- c('qname', 'sample', 'tax.identity', 'tax.identity.level')
      taxident.all   <- unique(taxident.all[, ..cols_by])
      taxident.all   <- merge(taxident.all, taxident.uni, by=c('tax.identity', 'tax.identity.level'))
      
      ## Summarise on ranks
      message('Summarising alignments on: ', paste(ranks, collapse = ', '), '...')
      columns_to_group_by <- c('sample', 'lineage', ranks, 'tax.identity', 'tax.identity.level')
      taxident.all.sum    <- taxident.all[   , .(count = .N), by = columns_to_group_by]
     
    #### #### 
    ##
      
      
      
      
    ##
    #### LCA method #### 
      
    # 1.) Az adott besorolási szinthez hozzá kell rendelni magát a tax.identity.level-t is ! (mert pl ha tax.identity.level == 'Class', akkor a Class rank az NA !!!)
    # 2.) species mindig van a progenomes-ban, ezzel kell valamit kezdeni?
    
    ## This will provide the tax identity for each read (LCA)
    } else if (methods == 'LCA') {
      
      message(
        'Method: LCA -->> This will provide the Lowest Common Ancestor (tax.identity) for each read that were kept after MAPQ filtering.')
      
      taxa.nodup <- taxa
      
      cols_by        <- c('qname', 'sample', 'tax.identity', 'tax.identity.level')
      taxa.nodupuni  <- unique(taxa.nodup[, ..cols_by] )  
      taxa.nodupuni  <- add_lineage(data = taxa.nodupuni, ranks, db.uni.data, only.tax.identity = T)
      taxa.nodupuni  <- unique(taxa.nodupuni)
      
      taxident.all   <- taxa.nodupuni
      
      ##
      print(list('Distribution of tax.identity.levels, after refining:', 
                 taxident.all[,.N,tax.identity.level]))
    
    
      ## Adding lineage to differentiate taxa with the same rank name
      cols_by        <- c(ranks, 'tax.identity', 'tax.identity.level')
      taxident.uni   <- unique(taxident.all[, ..cols_by])
      paste_columns_dt(taxident.uni, ranks, '::', 'lineage', 1)
      cols_by        <- c('qname', 'sample', 'tax.identity', 'tax.identity.level')
      taxident.all   <- unique(taxident.all[, ..cols_by])
      taxident.all   <- merge(taxident.all, taxident.uni, by=c('tax.identity', 'tax.identity.level'))
      
      ## Summarise on ranks
      message('Summarising alignments on: ', paste(ranks, collapse = ', '), '...')
      columns_to_group_by <- c('sample', 'lineage', ranks, 'tax.identity', 'tax.identity.level')
      taxident.all.sum    <- taxident.all[   , .(count = .N), by = columns_to_group_by]
      
      
    #### #### 
    ##
      
      
      
    ##
    #### Species Estimation Method ####
    
    ## This will keep every alignment and divide 
    ## their count back with the number of reads
    } else if (methods == 'SpeciesEstimate') {
      message(
      'Method: SpeciesEstimate -->> This will use every alignment of reads that were kept after MAPQ filtering and
      normalize the counts based on the number alignments.')
      
      ## calculate normalization factor
      readcount   <- n_distinct(taxa[,qname])
      norm_factor <- readcount / nrow(taxa)
      ## so the norm factor for each alignment is the ratio of the number of alignments and the reads 
      
     
      taxident.all  <- taxa
      
      ## Dereplicate reads, without alignments, keeping tax.identity only
      cols_by       <- c('qname', 'sample', ranks, 'tax.identity', 'tax.identity.level')
      taxident.all  <- unique(taxident.all[,..cols_by] )  
      taxident.all  <- unique(taxident.all)
      
      ##
      print(list('Distribution of tax.identity.levels, after refining:', 
                 taxident.all[,.N,tax.identity.level]))
      
      ## Adding lineage to differentiate taxa with the same rank name
      paste_columns_dt(taxident.all, ranks, '::', 'lineage', 3)
      
      
      ## Summarise on ranks and with tax identity
      message('Summarising alignments on: ', paste(ranks, collapse = ', '), '...')
      columns_to_group_by <- c('sample', 'lineage', ranks, 'tax.identity', 'tax.identity.level')
      taxident.all    <- taxident.all[   , .(count = .N), by = columns_to_group_by]
      
      
      ## Summarise taxon (rank) results, without tax identity, and normalize
      cols_by <- c('sample', ranks)
      taxident.all.sum <- taxa[   , .(norm_count = .N * norm_factor), by = cols_by]
      taxident.all.sum[taxident.sum == ''] <- NA
      
      ### 
    }
    ####
  
    })
    
    #### crop numbers for progenomes
    if (grepl('proGcontigs', db) & cropnumbers & exists('taxident.sum')) {
      taxident.sum <- as.data.frame(taxident.sum)
      taxident.sum[,ranks]       <- as.data.frame(apply(taxident.sum[,ranks], 2, function(x) str_replace(as.character(unlist(x)),  pattern = "^\\d+\\s", replacement = "")))
      taxident.sum$tax.identity  <- str_replace(as.character(unlist(taxident.sum$tax.identity)),  pattern = "^\\d+\\s", replacement = "")
      taxident.sum$count  <- as.integer(taxident.sum$count)
    }  
    
    #### Write outputs
    if (any(outputs %in% 'sum_taxa' )) { 
      fwrite(taxident.all.sum, paste0(outdir, '/', 'sum_taxa/', sample,  '_', methods, '_tax.abund.tsv'),    sep = '\t') 
    }
    
    if (any(outputs %in% 'sum_taxa' )) { 
      fwrite(taxident.all,     paste0(outdir, '/', 'sum_taxa/', sample,  '_', methods, '_tax.identity.tsv'), sep = '\t')
    }
    
    #plyr::count(taxident.sum$tax.identity.level)
    
  }
  ####  ####
  
  try({
    
    aln_summary <- data.frame(
           raw_bam       = nrow(bam),
           filt_bam      = nrow(bam.filt),
           mapq_filt     = nrow(minimap2),
           cigar_filt    = nrow(minimap2.filt),
           taxident_filt = nrow(taxa),
           final         = nrow(taxident.all))
    
   
    print(list('Number of alignments in each step:', aln_summary))
    fwrite(aln_summary,     paste0(outdir, '/', sample,  '_', methods, '_aln_summary.tsv'), sep = '\t')
    
  })
    
  
  })
  ## try end
  message('Finsished workflow for: ', sample, ', at: ', Sys.time())
  
  #### Outputs ####
         
         if(any(steps %in% '4.Summarise') & nrow(taxident.all.sum) > 0) { 
    message('Returning result of Step 4. for ', sample)   
    return(taxident.all.sum)
            
  } else if(any(steps %in% '3.Find_taxidentity') & nrow(taxa) > 0) {
    message('Returning result of Step 3. for ', sample)   
    return(taxa)
    
  } else if(any(steps %in% '2.Refine') & nrow(minimap2.filt) > 0) { 
    message('Returning result of Step 2. for ', sample)   
    return(minimap2.filt)
    
  } else if(any(steps %in% '1.Import_Alns') & nrow(bam.filt) > 0) { 
    message('Returning result of Step 1. for ', sample)   
    return(bam.filt)
    
  } else {
    message('minitax yielded no ouptuts for: ', sample, '!')
    return(NULL)
  }
  ####  ####
    
  
}



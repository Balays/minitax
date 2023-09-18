
wrap.fun.complete <- function(taxa_or_bamfile, pattern='.bam', sample=NULL,
                              input='bamfile', ## possible inputs: c('DT', 'list', 'bamfile', '.rds'), 
                              #DT_or_List=NULL,
                              config=config, 
                              steps=c('1.Import_Alns', '2.Find_taxidentity', '3.Refine', '4.Summarise'),
                              CIGAR_points=NULL, 
                              methods=NULL, 
                              saveRDS=F, cropnumbers=F,
                              outputs=unlist(strsplit(config$outputs, ', ')),
                              best.mapq=T, mapq.filt=NA,
                              db='proGcontigs', db.data=prog.db,  db.uni.data=prog.db.uni,
                              ranks = c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
                              
                              ) {
  
  #### Inputs ####  --->>> This part is not final
          if (input == 'bamfile') {
    sample   <- gsub(pattern, '', gsub('.*/', '', taxa_or_bamfile))
    bamfile  <- taxa_or_bamfile
    ##
  } else if (input == 'DT') {
    #sample  <- gsub(pattern, '', gsub('.*/', '', taxa_or_bamfile))
    DT      <- data.table(taxa_or_bamfile[taxa_or_bamfile$sample == sample, ])
    ## which step ???
  } else if (input == '.rds') {
    #sample  <- gsub(pattern, '', gsub('.*/', '', taxa_or_bamfile))
    DT      <- readRDS(paste0(outdir, '/', sample, ".minitax.v2.rds"))
    ## which step ???
  }
  
  message('Starting minitax WF on: ', sample, ', at: ', Sys.time())
  message('The following steps will be performed: ', paste(steps, collapse = '; '))
  
  #### methods
  ## avilable methods: 'LCA', 'BestAln', 'SpeciesEstimate'
  ## methods <- 'SpeciesEstimate'
  if (methods %in% c('LCA', 'SpeciesEstimate')) {
    calc.cigar.score <- F
    keep.max.cigar   <- F
  } else if (methods %in% 'BestAln') {
    calc.cigar.score <- T
    keep.max.cigar   <- T
  } else if (is.na(methods)) {
    calc.cigar.score <- T
    keep.max.cigar   <- F
  } 
  
  ## try start
  try({ 
  #### wrap.fun0 ####
  if(any(steps %in% '1.Import_Alns')) {
    message('Importing ', bamfile, '...')
    bam <- ov.from.bam2(
      bamfile =  bamfile,
      #bamnames,
      #write.filtered = F, force.create=force.create, filtering=F
      rm.gaps.in.aln = F, mapq.filt = NA,
      flag.tokeep    = NA,    flag.tocrop = NA, seqnames.tofilt = NA,
      is.lortia = F, add.primes = F
    )
    bam$sample <- sample
    
    ## Filtering alignments for MAPQ
    if(!all(is.na(mapq.filt))) {
      message('Removing alignments with MAPQ of : ', paste0(mapq.filt, collapse = ';'))
      bam.filt <- bam[!is.element(bam$mapq, mapq.filt), ]
    } else { bam.filt <- bam }
    
    ## Selecting the best alignment for multi-alignment reads based on MAPQ
    if(best.mapq) {
      message('Selecting the best alignment for multi-alignment reads based on MAPQ ...')
    } else {
      message('Every alignment will be kept ...')
    }
    minimap2.filt <- get.best.aln(bam.filt, keep.chim = T, 
                                  rm.supplementary = F, 
                                  rm.secondary = F, 
                                  best.mapq = best.mapq, 
                                  give.qname.ID.to.secondary=F)
    
    
    ## summarise
    message('Summariseing alignments statistics ...')
    bam.dist <- bam.filt %>% distinct(sample, qname, mapq) # rname
    bam.dist$mapq[is.na(bam.dist$mapq)] <- 'NA-CIGAR'
    bam.dist <- bam.dist %>% group_by(sample, qname) %>% reframe(mapq, n_unique_aln=seq_along(mapq))
    bam.dist.sp <- bam.dist %>% spread(n_unique_aln, mapq)
    bam.dist.sp$is.all.NA <- apply(bam.dist.sp %>% dplyr::select(-c(qname)),            1, function(x) all(is.na(x) | x == 'NA-CIGAR')) 
    bam.dist.sp$max.mapq  <- apply(bam.dist.sp %>% dplyr::select(-c(qname, is.all.NA)), 1, function(x) max(x, na.rm = T)) 
    
    bam.sum  <- bam.dist.sp %>% group_by(sample, is.all.NA) %>% reframe(is.mapped=n())
    #bam.sum$sample <- sample
    
    ##
    if (any(outputs %in% 'bam.sum')) {write_tsv(bam.sum, paste0(outdir, '/', sample, '_bamstats.tsv'))}
   
  } else {
    #minimap2.filt <- DT
  }
  ####  ####
  
  #### 2. Find tax identity (minitax) ####
  ## running minitax -->> finding tax.identity and tax.identity level for each read
  if(any(steps %in% '2.Find_taxidentity')) {
    message('Finding tax.identity and tax.identity level for each read...')
    taxa <- data.frame(NULL); try({
    taxa <- minitax2(minimap2.filt, db=db, db.data=db.data, 
                     ranks = ranks, db.uni.data=db.uni.data,
                     )
    })
    taxa$sample <- sample
    
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
      
      columns_to_group_by <- c('sample', 'qname')
      message('Calculating CIGAR scores...')
      taxa$cigar_score <- unlist(future_lapply(taxa$cigar, cigar_score, 
                                               match_score           = CIGAR_points$match_score, 
                                               mismatch_score        = CIGAR_points$mismatch_score,
                                               insertion_score       = CIGAR_points$insertion_score, 
                                               deletion_score        = CIGAR_points$deletion_score, 
                                               gap_opening_penalty   = CIGAR_points$gap_opening_penalty, 
                                               gap_extension_penalty = CIGAR_points$gap_extension_penalty))
      # Add a new column 'max_cigar_score'
      taxa[, max_cigar := (cigar_score == max(cigar_score)), by = columns_to_group_by]
    }
   
    ### output
    if (any(outputs %in% 'alignments_w_taxa')) {write_tsv(taxa, paste0(outdir, '/', sample, '_alignments_w_taxa.tsv'))}
    
  } else {
    #if (exists('DT')) {taxa <- DT}
  }
  
  ####  ####
  
  #### 3. Refine alignments ####
  #if (!exists('taxa')) {taxa <- DT}
  
  ### wrap.fun1
  if(any(steps %in% '3.Refine')) {
    message('Refining alignments...')
     
    ### wrap.fun2
    ## Select the rows with the highest 'cigar_score' for each 'qname'
    if (keep.max.cigar) {
      message('Keeping alignments with the highest CIGAR score for each read...')
      #max_cigar_score <- taxa_not_species[, .SD[which.max(cigar_score)], by = qname]
      taxa.filt <- taxa[taxa$max_cigar == T, ]
    } else {taxa.filt <- taxa}
    if (any(outputs %in% 'best_alignments_w_taxa')) {write_tsv(taxa.filt, paste0(outdir, '/', sample, '_best_alignments_w_taxa.tsv'))}
    
  } else {
    #if (exists('taxa')) {taxa <- taxa}
  }
  ####  ####
  
  
  #### 4.Summarise  #### 
  if(any(steps %in% '4.Summarise')) {
    message('Summarising taxa on the follwing taxon levels: ', paste0(ranks, collapse='; '))
    
    #
    
    ##
    #### BestAln method
    ## This will keep the best matching alignment for each read, based on cigar_score
    if (methods == 'BestAln') {
      message(
      'Method: BestAln -->> This will calculate CIGAR scores and use the best alignent for each read based on them.
      In the case when both the MAPQ and the CIGAR scores are equal, a random alignement will be selected.')
      
      columns_to_group_by <- c('sample', 'tax.identity')
      dupqnames  <- dup(taxa.filt$qname)
      #dupqnames  <- taxa.filt[is.element(taxa.filt$qname, dupqnames), ]
      message('Selecting a random alignement for ', 
              length(dupqnames), ' reads ',
              (length(dupqnames) / length(unique(taxa.filt$qname))) * 100, '%', 
              ', as both the MAPQ and the CIGAR scores were equal...')
      taxa.filt <- taxa.filt[!duplicated(taxa.filt$qname), ]
      ## replacing tax identity with species !!
      message('Replacing tax identity with species in these cases...')
      taxa.filt$tax.identity <- taxa.filt$species
      taxident.all.sum <- taxa.filt[   , .(count = .N), by = columns_to_group_by]
      
      #db.uni.rank      <- data.table(unique.data.frame(data.frame(db.uni.data)[,ranks]))
      db.uni.rank      <- unique(db.uni.data[,c("superkingdom", "phylum", "class", "order", "family", "genus", "species")])
      taxident.all.sum <- merge(taxident.all.sum, db.uni.rank, by.x='tax.identity', by.y='species')
      taxident.all.sum$species <- taxident.all.sum$tax.identity
      taxident.all.sum[taxident.all.sum == ''] <- NA
      taxident.sum  <- taxident.all.sum
      taxident.all  <- taxa.filt
      
    ####
    ##
    
    ##
    #### LCA method
      
      # 1.) Az adott besorolási szinthez hozzá kell rendelni magát a tax.identity.level-t is ! (mert pl ha tax.identity.level == 'Class', akkor a Class rank az NA !!!)
      # 2.) species mindig van a progenomes-ban, ezzel kell valamit kezdeni?
      
      
      
    ## This will provide the tax identity for each read (LCA)
    } else if (methods == 'LCA') {
      message(
        'Method: LCA -->> This will provide the Lowest Common Ancestor (tax.identity) for each read that were kept after MAPQ filtering.')
      
      columns_to_group_by <- c('sample', ranks, 'tax.identity', 'tax.identity.level')
      
      ## Merge taxon identity for each read with the database
      taxident.all     <- unique(taxa.filt [,c('sample', 'qname', 'tax.identity.level', 'tax.identity')])
      ## this contains the best assignment for each read
      
      ## Merge tax identity with the database, on each levels, to get rid of the best assignment and get the LCA for each read
      taxident.lca <- data.frame(NULL)
      for (i in seq_along(ranks)) {
        rank      <- ranks[i]
        rank_cols <- ranks[1:i]
        
        taxident       <- taxident.all[taxident.all$tax.identity.level == rank, ]
        db.uni.rank    <- unique(db.uni.data[,..rank_cols])
        
        taxident.rank  <- merge(taxident, db.uni.rank, by.x='tax.identity', by.y=rank)
        
        #if(i == 1) { 
          taxident.rank <- data.table(taxident.rank, RANK=taxident.rank$tax.identity)
          setnames(taxident.rank, old='RANK', new=rank)
        #} 
        
        if (i == length(ranks)) {
          taxident.rank[,ranks[length(ranks)]] <- taxident.rank$tax.identity
        }
        taxident.lca     <- plyr::rbind.fill(taxident.lca, taxident.rank)
      }
      taxident.lca <- data.table(taxident.lca)
      ## so in this table, for each taxonomic rank, we have the LCA
      
      ## Summarise tax identity table
      taxident.lca.sum <- taxident.lca[, .(count = .N), by = columns_to_group_by]
      ####
      ## 
      taxident.lca.sum[taxident.lca.sum == ''] <- NA
      taxident.sum <- taxident.lca.sum
    
      taxident.all  <- taxident.lca
      
    ####
    ##
    
    #### Species Estimation Method 
    ## This will keep every alignment and divide 
    ## their count back with the number of reads
    } else if (methods == 'SpeciesEstimate') {
      message(
      'Method: SpeciesEstimate -->> This will use every alignment of reads that were kept after MAPQ filtering and
      normalize the counts based on the number alignments.')
      
      ## 
      readcount   <- n_distinct(taxa.filt$qname)
      norm_factor <- readcount / nrow(taxa.filt)
      ## so the norm factor for each alignment is the ratio of the number of alignments and the reads 
      columns_to_group_by <- c('sample', ranks)
      taxident.spec.est.sum <- taxa.filt[   , .(norm_count = .N * norm_factor), by = columns_to_group_by]
      
      ## Summarise tax identity table
      taxident.spec.est.sum[taxident.spec.est.sum == ''] <- NA
      ### 
      taxident.sum <- taxident.spec.est.sum
      
      taxident.all <- taxa.filt 
    }
    ####
    
    
    #### crop numbers for progenomes
    if (grepl('proGcontigs', db) & cropnumbers & exists('taxident.sum')) {
      taxident.sum <- as.data.frame(taxident.sum)
      taxident.sum[,ranks]       <- as.data.frame(apply(taxident.sum[,ranks], 2, function(x) str_replace(as.character(unlist(x)),  pattern = "^\\d+\\s", replacement = "")))
      taxident.sum$tax.identity  <- str_replace(as.character(unlist(taxident.sum$tax.identity)),  pattern = "^\\d+\\s", replacement = "")
      taxident.sum$count  <- as.integer(taxident.sum$count)
    }  
    
    #### Write summary table
    if (any(outputs %in% 'sum_taxa' )) {write_tsv(taxident.sum, paste0(outdir, '/', methods, '_', sample, '_tax.abund.tsv'))}
    #plyr::count(taxident.sum$tax.identity.level)
  }
  ####  ####
  })
  ## try end
    
  #### Outputs ####
         if (exists('taxident.sum')) {
    return(taxident.sum)
    message('returning result of Step 4.')
  } else if (exists('taxident.all')) {
    return(taxident.all)
    message('returning result of Step 3.')
  } else if (exists('taxa')) {
    return(taxa)
    message('returning result of Step 2.')
  } else if (exists('minimap2.filt')) {
    return(minimap2.filt)
    message('returning result of Step 1.')
  } else if (exists('bam')) {
    return(bam)
    message('returning result of Step 0.')
  } else {
    message('mintax yielded no ouptuts for: ', sample, '!')
  }
  ####  ####
    
  message('Finsished workflow for: ', sample, '\n')
}


















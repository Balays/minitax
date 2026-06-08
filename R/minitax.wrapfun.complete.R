minitax_is_config_na <- function(value) {
  length(value) == 0 || is.na(value) || !nzchar(trimws(value)) || toupper(trimws(value)) %in% c("NA", "NULL")
}

minitax_split_config_values <- function(value, default = character()) {
  if (is.null(value) || minitax_is_config_na(value)) return(default)
  values <- trimws(unlist(strsplit(as.character(value), "[;,]")))
  values[nzchar(values)]
}

minitax_normalize_outputs <- function(outputs = NULL, config = NULL) {
  if (!is.null(outputs) && length(outputs) > 0) {
    values <- unlist(lapply(outputs, minitax_split_config_values))
    values <- values[nzchar(values)]
    if (length(values) > 0) return(unique(values))
  }
  if (!is.null(config) && 'outputs' %in% names(config)) {
    return(minitax_split_config_values(config$outputs, default = 'bam.sum'))
  }
  'bam.sum'
}

minitax_outdir <- function() {
  if (exists('outdir', inherits = TRUE)) return(get('outdir', inherits = TRUE))
  stop('Global outdir is not defined.', call. = FALSE)
}

minitax_sample_from_bam <- function(bamfile) {
  sub('\\.bam$', '', basename(bamfile))
}

minitax_cache_file <- function(sample, outdir = minitax_outdir()) {
  file.path(outdir, 'best_alignments_w_taxa', paste0(sample, '_best_alignments_w_taxa.tsv'))
}

minitax_cache_meta_file <- function(sample, outdir = minitax_outdir()) {
  file.path(outdir, 'best_alignments_w_taxa', paste0(sample, '_best_alignments_w_taxa.meta.tsv'))
}

minitax_scalar <- function(value) {
  paste(as.character(unlist(value)), collapse = ';')
}

minitax_cache_signature <- function(bamfile, config, db, mapq.filt, best.mapq,
                                    keep.max.cigar, CIGAR_points, ranks) {
  bam.info <- file.info(bamfile)
  db.dir <- if (!is.null(config) && 'db.dir' %in% names(config)) config$db.dir else ''
  signature <- data.table(
    cache_key = c('bamfile', 'bam_size', 'bam_mtime', 'db', 'db.dir', 'mapq.filt',
                  'best.mapq', 'keep.max.cigar', 'CIGAR_points', 'ranks'),
    value = c(
      normalizePath(bamfile, mustWork = FALSE),
      minitax_scalar(bam.info$size),
      if (is.na(bam.info$mtime)) '' else format(bam.info$mtime, tz = 'UTC', usetz = TRUE),
      minitax_scalar(db),
      normalizePath(db.dir, mustWork = FALSE),
      minitax_scalar(mapq.filt),
      minitax_scalar(best.mapq),
      minitax_scalar(keep.max.cigar),
      minitax_scalar(CIGAR_points),
      minitax_scalar(ranks)
    )
  )
  setnames(signature, 'cache_key', 'key')
  signature
}

minitax_cache_signature_matches <- function(meta.file, signature) {
  if (!file.exists(meta.file)) return(FALSE)
  old <- tryCatch(fread(meta.file, sep = '\t', colClasses = 'character'),
                  error = function(e) data.table())
  if (!all(c('key', 'value') %in% colnames(old))) return(FALSE)
  old <- old[match(signature$key, key)]
  all(!is.na(old$key)) && identical(old$value, signature$value)
}

minitax_write_cache_signature <- function(meta.file, signature) {
  dir.create(dirname(meta.file), recursive = TRUE, showWarnings = FALSE)
  fwrite(signature, meta.file, sep = '\t', na = 'NA')
}

minitax_score_cigars <- function(cigars, CIGAR_points) {
  unique.cigars <- unique(cigars)
  scores <- vapply(
    unique.cigars,
    cigar_score,
    numeric(1),
    match_score = CIGAR_points$match_score,
    mismatch_score = CIGAR_points$mismatch_score,
    insertion_score = CIGAR_points$insertion_score,
    deletion_score = CIGAR_points$deletion_score,
    gap_opening_penalty = CIGAR_points$gap_opening_penalty,
    gap_extension_penalty = CIGAR_points$gap_extension_penalty
  )
  unname(scores[match(cigars, unique.cigars)])
}

minitax_with_dt_threads <- function(dt.threads = 1L) {
  dt.threads <- as.integer(dt.threads)
  if (is.na(dt.threads) || dt.threads < 1L) dt.threads <- 1L
  old <- data.table::getDTthreads()
  data.table::setDTthreads(dt.threads)
  old
}

prepare_minitax_taxa_cache <- function(taxa_or_bamfile, config = config,
                                       keep.max.cigar = TRUE, CIGAR_points = NULL,
                                       methods = 'TaxaCache', outputs = NULL,
                                       best.mapq = TRUE, mapq.filt = NA,
                                       db = 'proGcontigs', db.data = prog.db,
                                       db.uni.data = prog.db.uni,
                                       ranks = c("superkingdom", "phylum", "class", "order", "family", "genus", "species"),
                                       outdir = minitax_outdir(),
                                       reuse.taxa.cache = TRUE, dt.threads = 1L) {
  old.threads <- minitax_with_dt_threads(dt.threads)
  on.exit(data.table::setDTthreads(old.threads), add = TRUE)

  started <- Sys.time()
  sample <- minitax_sample_from_bam(taxa_or_bamfile)
  cache.file <- minitax_cache_file(sample, outdir = outdir)
  meta.file <- minitax_cache_meta_file(sample, outdir = outdir)
  outputs <- unique(c(minitax_normalize_outputs(outputs, config), 'best_alignments_w_taxa'))
  signature <- minitax_cache_signature(taxa_or_bamfile, config, db, mapq.filt,
                                       best.mapq, keep.max.cigar, CIGAR_points, ranks)

  if (isTRUE(reuse.taxa.cache) &&
      file.exists(cache.file) &&
      minitax_cache_signature_matches(meta.file, signature)) {
    message('Reusing valid best_alignments_w_taxa cache for ', sample, ': ', cache.file)
    finished <- Sys.time()
    return(list(
      cache_file = cache.file,
      timing = data.table(
        sample = sample,
        stage = 'cache',
        method = NA_character_,
        status = 'reused',
        elapsed_sec = as.numeric(difftime(finished, started, units = 'secs')),
        rows_in = NA_integer_,
        rows_out = NA_integer_,
        file = cache.file
      )
    ))
  }

  message('Building best_alignments_w_taxa cache for ', sample, ': ', cache.file)
  taxa <- wrap.fun.complete(
    taxa_or_bamfile,
    config = config,
    input = 'bamfile',
    methods = methods,
    keep.max.cigar = keep.max.cigar,
    CIGAR_points = CIGAR_points,
    mapq.filt = mapq.filt,
    outputs = outputs,
    saveRDS = FALSE,
    steps = c('1.Import_Alns', '2.Refine', '3.Find_taxidentity'),
    best.mapq = best.mapq,
    db = db,
    db.data = db.data,
    db.uni.data = db.uni.data,
    ranks = ranks,
    outdir = outdir
  )
  if (is.null(taxa) || nrow(taxa) == 0) {
    stop('No tax-annotated alignments were produced for ', sample, call. = FALSE)
  }
  if (!file.exists(cache.file)) {
    dir.create(dirname(cache.file), recursive = TRUE, showWarnings = FALSE)
    fwrite(taxa, cache.file, sep = '\t', na = 'NA')
  }
  minitax_write_cache_signature(meta.file, signature)
  finished <- Sys.time()
  list(
    cache_file = cache.file,
    timing = data.table(
      sample = sample,
      stage = 'cache',
      method = NA_character_,
      status = 'rebuilt',
      elapsed_sec = as.numeric(difftime(finished, started, units = 'secs')),
      rows_in = NA_integer_,
      rows_out = nrow(taxa),
      file = cache.file
    )
  )
}

summarise_minitax_taxa_file <- function(taxa_or_bamfile, config = config,
                                        methods = NULL,
                                        keep.max.cigar = TRUE, CIGAR_points = NULL,
                                        outputs = NULL, best.mapq = TRUE,
                                        mapq.filt = NA,
                                        db = 'proGcontigs', db.data = prog.db,
                                        db.uni.data = prog.db.uni,
                                        ranks = c("superkingdom", "phylum", "class", "order", "family", "genus", "species"),
                                        outdir = minitax_outdir(),
                                        dt.threads = 1L) {
  old.threads <- minitax_with_dt_threads(dt.threads)
  on.exit(data.table::setDTthreads(old.threads), add = TRUE)
  started <- Sys.time()
  sample <- gsub('_best_alignments_w_taxa.tsv$', '', basename(taxa_or_bamfile))
  taxa.sum <- wrap.fun.complete(
    taxa_or_bamfile,
    config = config,
    input = 'taxa.DT',
    methods = methods,
    keep.max.cigar = keep.max.cigar,
    CIGAR_points = CIGAR_points,
    mapq.filt = mapq.filt,
    outputs = minitax_normalize_outputs(outputs, config),
    saveRDS = FALSE,
    steps = c('4.Summarise'),
    best.mapq = best.mapq,
    db = db,
    db.data = db.data,
    db.uni.data = db.uni.data,
    ranks = ranks,
    outdir = outdir
  )
  finished <- Sys.time()
  list(
    taxa.sum = taxa.sum,
    timing = data.table(
      sample = sample,
      stage = 'summarise',
      method = methods,
      status = 'completed',
      elapsed_sec = as.numeric(difftime(finished, started, units = 'secs')),
      rows_in = NA_integer_,
      rows_out = if (is.null(taxa.sum)) 0L else nrow(taxa.sum),
      file = taxa_or_bamfile
    )
  )
}

wrap.fun.complete <- function(taxa_or_bamfile, pattern='.bam', sample=NULL,
                              input='bamfile', ## possible inputs: c('DT', 'list', 'bamfile', '.rds'),
                              #DT_or_List=NULL,
                              config=config,
                              steps=c('1.Import_Alns', '2.Find_taxidentity', '3.Refine', '4.Summarise'),
                              keep.max.cigar=T, CIGAR_points=NULL,
                              methods=NULL,
                              saveRDS=F, cropnumbers=F,
                              outputs=NULL,
                              best.mapq=T, mapq.filt=NA,
                              db='proGcontigs', db.data=prog.db,  db.uni.data=prog.db.uni,
                              ranks = c("superkingdom", "phylum", "class", "order", "family", "genus", "species"),
                              outdir = minitax_outdir()

                              ) {
  outputs <- minitax_normalize_outputs(outputs, config)
  bam <- data.table(NULL)
  bam.filt <- data.table(NULL)
  minimap2 <- data.table(NULL)
  minimap2.filt <- data.table(NULL)
  taxa <- data.table(NULL)
  taxident.all <- data.table(NULL)
  taxident.all.sum <- data.table(NULL)

  #### Inputs ####  --->>> This part is not final
          if (input == 'bamfile') {
    sample   <- sub('\\.bam$', '', basename(taxa_or_bamfile))
    bamfile  <- taxa_or_bamfile
    ##
  } else if (input == 'taxa.DT') {
    taxa  <- fread(taxa_or_bamfile, na.strings = '')
    taxa[taxa == ''] <- NA
    sample <- gsub('_best_alignments_w_taxa.tsv$', '', basename(taxa_or_bamfile))
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
  tryCatch({

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
    bam.dist <- unique(bam.filt[, .(sample, qname, mapq)])
    bam.sum <- bam.dist[, .(is.all.NA = all(is.na(mapq)),
                            max.mapq = suppressWarnings(max(mapq, na.rm = TRUE))),
                        by = .(sample, qname)]
    bam.sum[is.infinite(max.mapq), max.mapq := NA_real_]
    bam.sum <- bam.sum[, .(is.mapped = .N), by = .(sample, is.all.NA)]
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
      minimap2[, cigar_score := minitax_score_cigars(cigar, CIGAR_points)]
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
    } else {
      minimap2.filt <- minimap2
    }

    print(head(minimap2.filt))

  } else {
    minimap2.filt <- minimap2
    #if (exists('DT')) {taxa <- DT}
  }

  ####  ####


  #### 3. Find tax identity (minitax) ####
  ## running minitax -->> finding tax.identity and tax.identity level for each read
  if(any(steps %in% '3.Find_taxidentity')) {
    message('Finding tax.identity and tax.identity level for each read...')

    # Identify the last entry of the ranks vector
    #last_rank <- tail(ranks, n = 1)

    {
    taxa <- minitax2(minimap2.filt, db=db, db.uni.data=db.uni.data,
                     ranks = ranks )
    }

    taxa$sample <- sample
    taxa[,aln_nr := seq_along(cigar), by=.(sample, qname)]

    ## Where the tax.identity is below the lowest rank (superkingdom)
    taxa[is.na(tax.identity),       tax.identity       := 'unclassified']
    taxa[is.na(tax.identity.level), tax.identity.level := 'unclassified']

    ### write output
    if (any(outputs %in% 'best_alignments_w_taxa')) {
      dir.create(file.path(outdir, 'best_alignments_w_taxa'), recursive = TRUE, showWarnings = FALSE)
      fwrite(taxa, paste0(outdir, '/', 'best_alignments_w_taxa/', sample,  '_best_alignments_w_taxa.tsv'),
             sep = '\t', na = 'NA') }

  } else {
    ##
  }

  ####  ####



  #### 4.Refine hits and Summarise  ####
  if(any(steps %in% '4.Summarise')) {
    message('Summarising taxa on the follwing taxon levels: ', paste0(ranks, collapse='; '))


    ##
    {

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

        ### THRESHOLDS !!!
        threshold.values <- NULL
        if ('BestAln_thresholds' %in% names(config) &&
            !is.na(config$BestAln_thresholds) &&
            nzchar(trimws(config$BestAln_thresholds)) &&
            toupper(trimws(config$BestAln_thresholds)) != 'NA') {
          threshold.values <- as.numeric(trimws(unlist(strsplit(config$BestAln_thresholds, '[;,]'))))
          if (any(is.na(threshold.values))) {
            stop('BestAln_thresholds must contain numeric values.')
          }
        }
        if (is.null(threshold.values)) {
          threshold.values <- 0.6
        }
        if (length(threshold.values) == 1) {
          threshold.values <- rep(threshold.values, length(ranks))
        }
        if (length(threshold.values) != length(ranks)) {
          stop('BestAln_thresholds must contain either one value or one value per rank.')
        }
        thresholds  <- data.frame(rank=ranks, threshold = threshold.values)

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

    # 1.) Az adott besorolĂ„â€šĂ‹â€ˇsi szinthez hozzĂ„â€šĂ‹â€ˇ kell rendelni magĂ„â€šĂ‹â€ˇt a tax.identity.level-t is ! (mert pl ha tax.identity.level == 'Class', akkor a Class rank az NA !!!)
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

      ## Each read contributes a total weight of 1 split across its retained alignments.
      taxident.all  <- copy(taxa)
      taxident.all[, aln_weight := 1 / .N, by = .(sample, qname)]

      ## Dereplicate reads, without alignments, keeping tax.identity only
      cols_by       <- c('qname', 'sample', ranks, 'tax.identity', 'tax.identity.level', 'aln_weight')
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
      taxident.all    <- taxident.all[   , .(count = sum(aln_weight)), by = columns_to_group_by]


      ## Summarise taxon (rank) results, without tax identity, and normalize
      cols_by <- c('sample', ranks)
      taxident.all.sum <- taxident.all[   , .(norm_count = sum(count)), by = cols_by]
      taxident.all.sum[taxident.all.sum == ''] <- NA

      ###
    }
    ####

    }

    #### crop numbers for progenomes
    if (grepl('proGcontigs', db) & cropnumbers & exists('taxident.all.sum')) {
      taxident.all.sum <- as.data.frame(taxident.all.sum)
      taxident.all.sum[,ranks] <- as.data.frame(apply(taxident.all.sum[,ranks], 2, function(x) str_replace(as.character(unlist(x)),  pattern = "^\\d+\\s", replacement = "")))
      if ('tax.identity' %in% colnames(taxident.all.sum)) {
        taxident.all.sum$tax.identity <- str_replace(as.character(unlist(taxident.all.sum$tax.identity)),  pattern = "^\\d+\\s", replacement = "")
      }
      if ('count' %in% colnames(taxident.all.sum)) {
        taxident.all.sum$count <- as.integer(taxident.all.sum$count)
      }
      setDT(taxident.all.sum)
    }

    #### Write outputs
    if (any(outputs %in% 'sum_taxa' )) {
      fwrite(taxident.all.sum, paste0(outdir, '/', 'sum_taxa/', sample,  '_', methods, '_tax.abund.tsv'),    sep = '\t')
    }

    if (any(outputs %in% 'sum_taxa' )) {
      fwrite(taxident.all,     paste0(outdir, '/', 'sum_taxa/', sample,  '_', methods, '_tax.identity.tsv'), sep = '\t')
    }

  }
  ####  ####

  try({

    method.label <- if (is.null(methods) || length(methods) == 0 || is.na(methods)) 'workflow' else methods
    aln_summary <- data.frame(
           raw_bam       = nrow(bam),
           filt_bam      = nrow(bam.filt),
           mapq_filt     = nrow(minimap2),
           cigar_filt    = nrow(minimap2.filt),
           taxident_filt = nrow(taxa),
           final         = nrow(taxident.all))


    print(list('Number of alignments in each step:', aln_summary))
    fwrite(aln_summary,     paste0(outdir, '/', sample,  '_', method.label, '_aln_summary.tsv'), sep = '\t')

  })


  }, error = function(e) {
    stop('minitax workflow failed for ', sample, ': ', conditionMessage(e), call. = FALSE)
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

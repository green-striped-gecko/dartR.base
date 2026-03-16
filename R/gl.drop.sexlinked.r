#'@name gl.drop.sexlinked
#'@title Removes loci that are sex linked
#'@description
#' This function identifies sex-linked and autosomal loci present in a SNP
#' dataset (genlight object) using individuals with known sex. It identifies
#' five types of loci: w-linked or y-linked, sex-biased, z-linked or
#' x-linked, gametologous and autosomal.
#'
#' This function produces as output a genlight object with autosomal loci only.
#'
#' @param x Name of the genlight object containing the SNP data. This genlight
#' object needs to contain the sex of the individuals. See explanation in
#' details [required].
#' @param system String that declares the sex-determination system of the
#' species: 'zw' or 'xy' [required].
#' @param ncores Number of processes to be used in parallel operation. If ncores
#' > 1 parallel operation is activated, see "Details" section [default 1].
#' @param plot.display Creates four output plots. See explanation in details
#' [default TRUE].
#' @param plot.theme Theme for the plot. See Details for options
#' [default theme_dartR()].[not yet implemented]
#' @param plot.colors [not implemented yet]
#' @param plot.dir Directory to save the plot RDS files [default as specified
#' by the global working directory or tempdir()].
#' @param plot.file Name for the RDS binary file to save (base name only,
#' exclude extension) [default NULL].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#'  [default NULL, unless specified using gl.set.verbosity].
#' @details
#' The genlight object must contain in \code{gl@other$ind.metrics} a column
#' named "id", and a column named "sex" in which individuals with known-sex are
#' assigned 'M' for male, or 'F' for female. The function ignores individuals
#' that are assigned anything else or nothing at all (unknown-sex).
#'
#' The creation of plots can be turned-off (\code{plot.display = FALSE}) in 
#' order to save a little bit of running time for very large datasets
#'  (>50,000 SNPs). However, we strongly encourage you to always inspect the 
#'  output plots at least once to make sure everything is working properly.
#'
#'\strong{ Function's output }
#'
#' This function returns as output a genlight object that contains only 
#' autosomal loci (i.e. sex-linked loci have been dropped.)
#'
#' And four plots:\itemize{
#' \item {A BEFORE plot based on loci call rate by sex, with w/y-linked loci
#'  colored in yellow and sex-biased loci in blue}
#' \item {An AFTER plot based on loci call rate by sex, with sex-linked loci 
#' removed}
#' \item {A BEFORE plot based on loci heterozygosity by sex, with z/x-linked 
#' loci colored in orange and gametologs in green}
#' \item {An AFTER plot based on loci heterozygosity by sex, with sex-linked 
#' loci removed}
#' }
#'
#' @return A genlight object and 4 plots.
#'
#' @author Custodian: Diana Robledo-Ruiz -- Post to
#'   \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' LBP_noSexLinked <- gl.drop.sexlinked(x = LBP, system = "xy", 
#' plot.display = TRUE, ncores = 1)
#' LBP_noSexLinked
#'
#' @references
#' \itemize{
#' \item Robledo-Ruiz, D. A., Austin, L., Amos, J. N., Castrejon-Figueroa, J.,
#'  Harley, D. K., Magrath, M. J., Sunnucks, P., & Pavlova, A. (2023).
#'  Easy to use R functions to separate reduced representation genomic datasets
#'   into sex linked and autosomal loci, and conduct sex assignment. Molecular
#'   Ecology Resources, 00, 1-21.
#'  }
#'
#' @importFrom stats chisq.test
#' @importFrom stats fisher.test
#' @importFrom stats p.adjust
#' @importFrom foreach foreach "%dopar%"
#' @import patchwork
#'
#' @export

gl.drop.sexlinked <- function(x,
                              system = NULL,
                              ncores = 1,
                              plot.display = TRUE,
                              plot.theme = theme_dartR(),
                              plot.colors = NULL,
                              plot.file = NULL,
                              plot.dir = NULL,
                              verbose = NULL) {
  # PRELIMINARIES -- checking ----------------
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  if (verbose == 0) {
    plot.display <- FALSE
  }
  
  # SET WORKING DIRECTORY
  plot.dir <- gl.check.wd(plot.dir, verbose = 0)
  
  # SET COLOURS #not yet implemented...
  if (is.null(plot.colors)) {
    plot.colors <- c("#2171B5", "#6BAED6")
  } else {
    if (length(plot.colors) > 2) {
      if (verbose >= 2) {
        cat(warn(
          "  More than 2 colors specified, only the first 2 are used\n"
        ))
      }
      plot.colors <- plot.colors[1:2]
    }
  }
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "v.2023.3",
                   verbose = verbose)
  
  # CHECK DATATYPE
  datatype <- utils.check.datatype(x,
                                   accept = c("dartR",
                                              "genlight",
                                              "SNP", 
                                              "SilicoDArT"),
                                   verbose = verbose)
  
  if (ncores > 1) {
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
  }
  
  if (is.null(system)) {
    stop(error(
      "You must specify the sex-determination system with the parameter 'system' ('zw' or 'xy')."
    ))
  } else {
    if (!(system == 'zw' | system == 'xy')) {
      stop(error(
        "Parameter 'system' must be 'zw' or 'xy'."
      ))
    }
  }
  
  # Transform genotypes to matrix and transpose
  gen <- as.data.frame(t(as.matrix(x)))
  
  # Extract IDs per sex
  metrics <- x@other$ind.metrics
  
  # Locate the "sex" column regardless of case
  sex.col <- grep("^sex$",
                  names(metrics),
                  ignore.case = TRUE,
                  value = TRUE)
  if (length(sex.col) != 1) {
    stop(error(
      "Could not find a single column named 'sex' (case-insensitive)."
    ))
  }
  
  # Pull values, force upper case, and drop NAs
  sex.values <- toupper(metrics[[sex.col]])
  valid      <- !is.na(sex.values)
  
  # Check that at least one F or M exists
  if (!any(sex.values[valid] %in% c("F", "M"))) {
    stop(error(
      "Females and males in the 'sex' column must be 'F' or 'M' (any case)."
    ))
  }
  
  # subset IDs
  ids.F <- metrics$id[sex.values == "F" & valid]
  ids.M <- metrics$id[sex.values == "M" & valid]
  
  if (verbose > 1){
    message(report(paste(
      "  Detected ",
      length(ids.F),
      " females and ",
      length(ids.M),
      " males.",
      sep = ""
    )))
  }
  
  # Subset genotypes by sex
  gen.F <- gen[, (colnames(gen) %in% ids.F)]
  gen.M <- gen[, (colnames(gen) %in% ids.M)]
  
  ##################### 1. Sex-linked loci by scoring rate
  
  # Create a results table with loci as row names and index number per locus
  table <- data.frame(index = c(1:nrow(gen)), row.names = row.names(gen))
  
  # Count missing (NA) and add as column to table
  table$count.F.miss <- rowSums(is.na(gen.F))
  table$count.M.miss <- rowSums(is.na(gen.M))
  
  # Count scored ("0", "1" or "2") and add as column to table
  table$count.F.scored <- rowSums(!is.na(gen.F))
  table$count.M.scored <- rowSums(!is.na(gen.M))
  
  if (ncores > 1) {
    if (verbose > 1){
      message(report(
        "  Starting phase 1. Working in parallel..."
      ))
    }
  } else {
    if (verbose > 1){
      message(report(
        "  Starting phase 1. May take a while..."
      ))
    }
  }
  
  # Apply Fisher's exact test (because there are observations with less than 5)
  if (ncores > 1) {
    xfisher <- foreach::foreach(i = 1:nrow(table), .combine = rbind) %dopar% {
      # Make vector of observed values
      obs <- matrix(
        c(table[i, "count.F.miss"], table[i, "count.M.miss"],
          table[i, "count.F.scored"], table[i, "count.M.scored"]),
        nrow = 2,
        ncol = 2,
        dimnames = list(c("F", "M"), c("miss", "scored"))
      )
      
      # See if it's possible to use chisq test
      if (sum(obs) >= 1000) {
        # Convert zeros to 1 to not obtain an error with chisq-test
        obs[obs == 0] <- 1
        
        # Chisq-test
        chisq.res <- chisq.test(obs, correct = FALSE)
        
        # Add to results table
        return(data.frame(
          ratio = chisq.res$statistic,
          p.value = chisq.res$p.value
        ))
        
      } else {
        # Convert zeros to 1 to not obtain an error with Fisher's test
        obs[obs == 0] <- 1
        
        # Run Fisher's exact test
        F.test <- fisher.test(obs)
        
        # Add to results table
        return(data.frame(
          ratio = F.test$estimate,
          p.value = F.test$p.value
        ))
      }
    }
    table <- cbind(table, xfisher)
  } else {
    # Add empty columns for chisq-statistic and corresponding p-value
    table$ratio   <- NA
    table$p.value <- NA
    
    # Test for independece of sex and missingness
    for (i in 1:nrow(table)) {
      # Make matrix of observed values
      obs <- matrix(
        c(table[i, "count.F.miss"], table[i, "count.M.miss"], 
          table[i, "count.F.scored"], table[i, "count.M.scored"]),
        nrow = 2,
        ncol = 2,
        dimnames = list(c("F", "M"), c("miss", "scored"))
      )
      
      # See if it's possible to use chisq test
      if (sum(obs) >= 1000) {
        # Convert zeros to 1 to not obtain an error with chisq-test
        obs[obs == 0] <- 1
        
        # Chisq-test
        chisq.res <- chisq.test(obs, correct = FALSE)
        
        # Add to results table
        table[i, "ratio"]   <- chisq.res$statistic
        table[i, "p.value"] <- chisq.res$p.value
        
      } else {
        # Convert zeros to 1 to not obtain an error with Fisher's test
        obs[obs == 0] <- 1
        
        # Run Fisher's exact test
        F.test <- fisher.test(obs)
        
        # Add to results table
        table[i, "ratio"]   <- F.test$estimate
        table[i, "p.value"] <- F.test$p.value
      }
    }
  }
  
  scoringRate.F <- scoringRate.M <- heterozygosity.F <- heterozygosity.M <- NA
  # Adjust p-values for multiple comparisons (False discovery rate)
  table$p.adjusted <- p.adjust(table$p.value, method = "fdr")
  
  # Calculate scoring rate for females and males and add to results table
  table$scoringRate.F <- table$count.F.scored / (table$count.F.scored +
                                                   table$count.F.miss)
  
  table$scoringRate.M <- table$count.M.scored / (table$count.M.scored +
                                                   table$count.M.miss)
  
  ##### 1.1 W-linked or Y-linked loci
  # For zw sex-determination system
  if (system == "zw") {
    table$w.linked <- NA
    
    for (i in 1:nrow(table)) {
      if (table[i, "scoringRate.M"] <= 0.1 &&
          table[i, "p.adjusted"] <= 0.01) {
        table[i, "w.linked"] <- TRUE
      } else {
        table[i, "w.linked"] <- FALSE
      }
    }
    table.wlinked <- table[table$w.linked == TRUE, ]
  }
  
  # For xy sex-determination system
  if (system == "xy") {
    table$y.linked <- NA
    
    for (i in 1:nrow(table)) {
      if (table[i, "scoringRate.F"] <= 0.1 &&
          table[i, "p.adjusted"] <= 0.01) {
        table[i, "y.linked"] <- TRUE
      } else {
        table[i, "y.linked"] <- FALSE
      }
    }
    table.ylinked <- table[table$y.linked == TRUE, ]
  }
  
  ##### 1.2 Loci with sex-biased scoring rate
  table$sex.biased <- NA
  
  for (i in 1:nrow(table)) {
    if (table[i, "p.adjusted"] <= 0.01 && table[i, 11] == FALSE) {
      table[i, "sex.biased"] <- TRUE
    } else {
      table[i, "sex.biased"] <- FALSE
    }
  }
  table.sexbiased <- table[table$sex.biased == TRUE, ]
  
  ##### 1.3 Plot BEFORE vs AFTER
  
  if (verbose > 1){
    message(report(
      "  Building call rate plots."
    ))
  }
  
  # For zw sex-determination system
  if (system == "zw") {
    table.autosomal <- table[table$w.linked == FALSE &
                               table$sex.biased == FALSE, ]
    
    BEF.mis <- ggplot2::ggplot(table.autosomal,
                               aes(x = scoringRate.F, y = scoringRate.M)) +
      geom_point(color = 'grey33') +
      geom_point(data = table.sexbiased, color = 'dodgerblue3') +
      geom_point(data = table.wlinked, color = 'gold') +
      ggtitle("BEFORE dropping sex-linked loci") +
      xlab("Call rate Females") +
      ylab("Call rate Males") +
      xlim(0, 1) + ylim(0, 1)
    
    AFT.mis <- ggplot2::ggplot(table.autosomal, 
                               aes(x = scoringRate.F, y = scoringRate.M)) +
      geom_point(color = 'grey33') +
      ggtitle("AFTER dropping sex-linked loci") +
      xlab("Call rate Females") +
      ylab("Call rate Males") +
      xlim(0, 1) + ylim(0, 1)
  }
  
  # For xy sex-determination system
  if (system == "xy") {
    table.autosomal <- table[table$y.linked == FALSE &
                               table$sex.biased == FALSE, ]
    
    BEF.mis <- ggplot2::ggplot(table.autosomal, 
                               aes(x = scoringRate.F, y = scoringRate.M)) +
      geom_point(color = 'grey33') +
      geom_point(data = table.sexbiased, color = 'dodgerblue3') +
      geom_point(data = table.ylinked, color = 'gold') +
      ggtitle("BEFORE dropping sex-linked loci") +
      xlab("Call rate Females") +
      ylab("Call rate Males") +
      xlim(0, 1) + ylim(0, 1)
    
    AFT.mis <- ggplot2::ggplot(table.autosomal, 
                               aes(x = scoringRate.F, y = scoringRate.M)) +
      geom_point(color = 'grey33') +
      ggtitle("AFTER dropping sex-linked loci") +
      xlab("Call rate Females") +
      ylab("Call rate Males") +
      xlim(0, 1) + ylim(0, 1)
  }
  
  #################### 2. Sex-linked loci by heterozygosity
  # Count heterozygotes ("1") and add as column to results table
  table$count.F.het <- rowSums(gen.F == 1, na.rm = TRUE)
  table$count.M.het <- rowSums(gen.M == 1, na.rm = TRUE)
  
  # Count homozygotes ("0" or "2") and add as column to results table
  table$count.F.hom <- rowSums(
    gen.F != 1,
    # Ignores NAs
    na.rm = TRUE
  )
  table$count.M.hom <- rowSums(
    gen.M != 1,
    # Ignores NAs
    na.rm = TRUE
  )
  
  if (verbose > 1){
    message(report(
      "  Starting phase 2. May take a while..."
    ))
  }
  
  if (ncores > 1) {
    xstat <- foreach::foreach(i = 1:nrow(table), .combine = rbind) %dopar% {
      # Exclude w.y-linked loci and loci with sex-biased score
      if (table[i, 11] == TRUE | table[i, "sex.biased"] == TRUE) {
        stat.value <- NA
        stat.p.value <- NA
        
      } else {
        # Make contingency table
        contingency <- matrix(
          c(table[i, "count.F.het"], table[i, "count.M.het"], 
            table[i, "count.F.hom"], table[i, "count.M.hom"]),
          nrow = 2,
          ncol = 2,
          dimnames = list(c("F", "M"), c("het", "hom"))
        )
        
        # Check if Yate's correction is necessary (n =< 20, Sokhal & Rohlf 1995)
        if (sum(contingency) >= 1000) {
          # Convert zeros to 1 to not obtain an error with chisq-test
          contingency[contingency == 0] <- 1
          
          # Chisq-test
          chisq.res <- chisq.test(contingency, correct = FALSE)
          
          # Add to results table
          stat.value   <- chisq.res$statistic
          stat.p.value <- chisq.res$p.value
          
        } else {
          # Convert zeros to 1 to not obtain an error with chisq-test
          contingency[contingency == 0] <- 1
          
          # Run Fisher's exact test
          F.test <- fisher.test(contingency)
          
          # Add to results table
          stat.value   <- F.test$estimate
          stat.p.value <- F.test$p.value
        }
      }
      return(data.frame(stat = stat.value, stat.p.value = stat.p.value))
    }
    table <- cbind(table, xstat)
  } else {
    # Add empty columns for statistic and corresponding p-value
    table$stat         <- NA
    table$stat.p.value <- NA
    
    # Apply test for independence of sex and heterozygosity
    for (i in 1:nrow(table)) {
      # Exclude w.y-linked loci and loci with sex-biased score
      if (table[i, 11] == TRUE | table[i, "sex.biased"] == TRUE) {
        table[i, "stat"]         <- NA
        table[i, "stat.p.value"] <- NA
        
      } else {
        # Make contingency table
        contingency <- matrix(
          c(table[i, "count.F.het"], table[i, "count.M.het"], 
            table[i, "count.F.hom"], table[i, "count.M.hom"]),
          nrow = 2,
          ncol = 2,
          dimnames = list(c("F", "M"), c("het", "hom"))
        )
        
        # See if it's possible to use chisq test
        if (sum(contingency) >= 1000) {
          # Convert zeros to 1 to not obtain an error with chisq-test
          contingency[contingency == 0] <- 1
          
          # Chisq-test
          chisq.res <- chisq.test(contingency, correct = FALSE)
          
          # Add to results table
          table[i, "stat"]         <- chisq.res$statistic
          table[i, "stat.p.value"] <- chisq.res$p.value
          
        } else {
          # Convert zeros to 1 to not obtain an error with Fisher's test
          contingency[contingency == 0] <- 1
          
          # Run Fisher's exact test
          F.test <- fisher.test(contingency)
          
          # Add to results table
          table[i, "stat"]         <- F.test$estimate
          table[i, "stat.p.value"] <- F.test$p.value
        }
      }
    }
  }
  
  # Adjust p-values for multiple comparisons (False discovery rate, least 
  # conservative)
  table$stat.p.adjusted <- p.adjust(table$stat.p.value, method = "fdr")
  
  # Calculate for heterozygosity per sex and add to results table
  table$heterozygosity.F <- table$count.F.het / 
    (table$count.F.het + table$count.F.hom)
  
  table$heterozygosity.M <- table$count.M.het /
    (table$count.M.het + table$count.M.hom)
  
  ##### 2.1 Z-linked or X-linked loci AND gametologs
  # For zw sex-determination system
  if (system == "zw") {
    table$z.linked  <- FALSE
    table$gametolog <- FALSE
    
    for (i in 1:nrow(table)) {
      # Exclude w-linked loci and loci with sex-biased score
      if (!is.na(table[i, "stat.p.adjusted"])) {
        # Exclude autosomal
        if (table[i, "stat.p.adjusted"] <= 0.01) {
          # Identify if heterozygosity is larger in males
          if (table[i, "heterozygosity.M"] > table[i, "heterozygosity.F"]) {
            table[i, "z.linked"] <- TRUE
          } else {
            table[i, "gametolog"] <- TRUE
          }
        }
      }
    }
    table.zlinked <- table[table$z.linked  == TRUE, ]
    table.gametol <- table[table$gametolog == TRUE, ]
  }
  
  # For xy sex-determination system
  if (system == "xy") {
    table$x.linked  <- FALSE
    table$gametolog <- FALSE
    
    for (i in 1:nrow(table)) {
      # Exclude y-linked loci and loci with sex-biased score
      if (!is.na(table[i, "stat.p.adjusted"])) {
        # Exclude autosomal
        if (table[i, "stat.p.adjusted"] <= 0.01) {
          # Identify if heterozygosity is larger in females
          if (table[i, "heterozygosity.F"] > table[i, "heterozygosity.M"]) {
            table[i, "x.linked"] <- TRUE
          } else {
            table[i, "gametolog"] <- TRUE
          }
        }
      }
    }
    table.xlinked <- table[table$x.linked  == TRUE, ]
    table.gametol <- table[table$gametolog == TRUE, ]
  }
  
  ##### 2.2 Plot BEFORE vs AFTER
  
  if (verbose > 1){
    message(report(
      "  Building heterozygosity plots."
    ))
  }
  
  # For zw sex-determination system
  if (system == "zw") {
    table.autosomal <- table[table$w.linked   == FALSE &
                               table$sex.biased == FALSE &
                               table$z.linked   == FALSE &
                               table$gametolog  == FALSE, ]
    
    BEF.het <- ggplot2::ggplot(table.autosomal,
                               aes(x = heterozygosity.F, y = heterozygosity.M)) +
      geom_point(color = 'grey33') +
      geom_point(data = table.gametol, color = 'chartreuse3') +
      geom_point(data = table.zlinked, color = 'darkorange1') +
      ggtitle("BEFORE dropping sex-linked loci") +
      xlab("% Heterozygous Females") +
      ylab("% Heterozygous Males") +
      xlim(0, 1) + ylim(0, 1)
    
    AFT.het <- ggplot2::ggplot(table.autosomal,
                               aes(x = heterozygosity.F, y = heterozygosity.M)) +
      geom_point(color = 'grey33') +
      ggtitle("AFTER dropping sex-linked loci") +
      xlab("% Heterozygous Females") +
      ylab("% Heterozygous Males") +
      xlim(0, 1) + ylim(0, 1)
  }
  
  # For xy sex-determination system
  if (system == "xy") {
    table.autosomal <- table[table$y.linked   == FALSE &
                               table$sex.biased == FALSE &
                               table$x.linked   == FALSE &
                               table$gametolog  == FALSE, ]
    
    BEF.het <- ggplot2::ggplot(table.autosomal,
                               aes(x = heterozygosity.F, y = heterozygosity.M)) +
      geom_point(color = 'grey33') +
      geom_point(data = table.gametol, color = 'chartreuse3') +
      geom_point(data = table.xlinked, color = 'darkorange1') +
      ggtitle("BEFORE dropping sex-linked loci") +
      xlab("% Heterozygous Females") +
      ylab("% Heterozygous Males") +
      xlim(0, 1) + ylim(0, 1)
    
    AFT.het <- ggplot2::ggplot(table.autosomal,
                               aes(x = heterozygosity.F, y = heterozygosity.M)) +
      geom_point(color = 'grey33') +
      ggtitle("AFTER dropping sex-linked loci") +
      xlab("% Heterozygous Females") +
      ylab("% Heterozygous Males") +
      xlim(0, 1) + ylim(0, 1)
  }
  
  if (verbose > 1){
    message(report(
      "  Done building heterozygosity plots."
    ))
  }
  
  #################### 3. Create output of function
  ##### 3.1 Save the indices of each category of loci to later subset x
  # For zw sex-determination system
  if (system == "zw") {
    a <- table[table$w.linked   == TRUE, "index"]
    b <- table[table$sex.biased == TRUE, "index"]
    c <- table[table$z.linked   == TRUE, "index"]
    d <- table[table$gametolog  == TRUE, "index"]
    
    autosomal <- table[table$w.linked   == FALSE &
                         table$sex.biased == FALSE &
                         table$z.linked   == FALSE &
                         table$gametolog  == FALSE, "index"]
    
    if (verbose>1) message("**FINISHED** \nTotal of analysed loci: ", nrow(table), ".\n",
            "Dropped ", length(a)+length(b)+length(c)+length(d), " sex-linked loci:\n",
            "   ",    length(a), " W-linked loci (yellow)\n",
            "   ",    length(b), " sex-biased loci (blue)\n",
            "   ",    length(c), " Z-linked loci (orange)\n",
            "   ",    length(d), " gametologs (green).\n",
            "And kept ",   length(autosomal), " autosomal loci (grey).")
  }
  
  if (system == "xy") {
    a <- table[table$y.linked   == TRUE, "index"]
    b <- table[table$sex.biased == TRUE, "index"]
    c <- table[table$x.linked   == TRUE, "index"]
    d <- table[table$gametolog  == TRUE, "index"]
    
    autosomal <- table[table$y.linked   == FALSE &
                         table$sex.biased == FALSE &
                         table$x.linked   == FALSE &
                         table$gametolog  == FALSE, "index"]
    
    if (verbose>1) message("**FINISHED** \nTotal of analyzed loci: ", nrow(table), ".\n",
                           "Dropped ", length(a)+length(b)+length(c)+length(d), " sex-linked loci:\n",
                           "   ",    length(a), " Y-linked loci (yellow)\n",
                           "   ",    length(b), " sex-biased loci (blue)\n",
                           "   ",    length(c), " X-linked loci (orange)\n",
                           "   ",    length(d), " gametologs (green).\n",
                           "And kept ",   length(autosomal), " autosomal loci (grey).")
  }
  
  ##### 3.2 Subset x object
  gl.autosomal <- x[, autosomal]
  
  #################### 4. Output
  if (ncores > 1) {
    parallel::stopCluster(cl)
  }
  
  p4 <- BEF.mis + AFT.mis + BEF.het + AFT.het
  
  if (plot.display) {
    print(p4)
  }
  
  # Optionally save the plot ---------------------
  
  if (!is.null(plot.file)) {
    tmp <- utils.plot.save(p4,
                           dir = plot.dir,
                           file = plot.file,
                           verbose = verbose)
  }
  
  # FLAG SCRIPT END ---------------
  
  if (verbose >= 1) {
    cat(report(
      "Completed:", funname, "\n"
    ))
  }
  # ----------------------
  
  # RETURN
  
  return(gl.autosomal)
}
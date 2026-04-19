#' @name sq.get.table
#' @title pull a table from a MariaDB database
#' @family utility

#' @description
#' This is a utility function that pulls a table from a MariaDB database
#'
#' @param con -- name of connection to the database [default con]
#' @param table -- name of the table to import [required]
#' @param verbose -- TRUE, progress displayed; FALSE, progress lines suppressed [default FALSE] 
#' @return returns a dataframe with the table contents
# @import RMariaDB
#' @export
#' @author Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}

sq.get.table <- function(con = con, 
                         table = NULL, 
                         verbose=FALSE) {
  
  
  if(is.null(con)){
    stop("FATAL ERROR: Database connection not specified, use sq.open.db to open a database\n")
  }
  if(is.null(table)){
    stop("FATAL ERROR: Database table not specified, no default\n")
  }
  
  sqlstr <- paste0("SELECT * FROM ",table)
  df <- RMariaDB::dbGetQuery(con,sqlstr)
  #df[df==''] <- "Null"
  df[apply(df,MARGIN=2, FUN=as.character) == ''] <- "Null"
  df[apply(df,MARGIN=2, FUN=as.character) == 'null'] <- "Null"
  df[apply(df,MARGIN=2, FUN=as.character) == 'NULL'] <- "Null"
  #df[apply(df,MARGIN=2, FUN=as.character) == "-001-11-30"] <- "Null"
  df[apply(df,MARGIN=2, FUN=as.character) == 'NA'] <- "Null"
  df[apply(df,MARGIN=2, FUN=as.character) == ' '] <- "Null"
  df[apply(df,MARGIN=2, FUN=as.character) == 'Null'] <- NA
  if(table=="C_clutches"){
    df <- df[!is.na(df$ClutchID),]
  } else if(table=="Z_cagelist" | table=="Z_freezerlist" | table=="Z_valuelist" | table=="Z_synonyms"){
    df <- df
  } else {
    df <- df[!is.na(df$SpecimenID),]
  }
  
  if(verbose){
    cat("Fields in table",table,"\n")
    cat(paste("   ",names(df),"\n"))
  }
  
  return(df)
  
}


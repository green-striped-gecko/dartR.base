#' @name sq.open.db
#' @title open a MariaDB database
#' @family utility
#'
#' @description
#' This is a utility function that opens a MariaDB connection for a specified database
#'
#' @param db -- name of the database to open [required]
#' @param table -- name of the table to update [required]
#' @param host -- name of the host for the database [required]
#' @param user -- name of the user to use to access the database [required]
#' @param passwd -- password to use to access the database [required]
#' @param verbose -- TRUE, progress displayed; FALSE, progress lines suppressed [default FALSE] 
#' @return returns a connection to the database (con)
# @import RMariaDB
#' @export
#' @author Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' 

sq.open.db <- function(host,
                       db,
                       table,
                       user,
                       passwd,
                       verbose=FALSE) {
  
  con <- RMariaDB::dbConnect(RMariaDB::MariaDB(), 
                   user=user, 
                   password=passwd, 
                   dbname=db, 
                   host=host)
  if(db!="mysql"){
    if(verbose){
      cat("Tables in database",db,"\n")
      cat(paste("  ",RMariaDB::dbListTables(con),"\n"))
    }
  }
  
  return(con)
}


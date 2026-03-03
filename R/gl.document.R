#' @name gl.document
#' @title Generate a roxygen2 Documentation Template for a Function
#' @description
#' Creates a skeleton \code{roxygen2} documentation file for a specified 
#' function. The generated file contains standard documentation fields 
#' including title, description, parameters, details, return value, examples, 
#' and references. The output is written as a new \code{.R} file in the 
#' specified directory.
#'
#' The function inspects the formal arguments of the target function and 
#' automatically generates \code{@param} entries for each argument. This 
#' provides a structured starting point for developing consistent 
#' documentation across a package.
#'
#' @param func_name Name of the function to be documented (unquoted).
#' @param author_name Character string specifying the author of the function.
#' @param example_dataset Name of an example dataset to include in the 
#' documentation examples (currently placeholder, not implemented).
#' @param outputDir Character string specifying the directory where the 
#' documentation file will be written. The file will be named 
#' \code{<func_name>.R}.
#'
#' @details 
#' The function constructs a list of standard documentation fields and writes 
#' them in roxygen2 format to a new file connection. Parameter names are 
#' extracted using \code{formals()} and written as placeholder entries for 
#' subsequent manual editing.
#' 
#' The generated template includes the following sections:
#' \itemize{
#'   \item \code{@name}
#'   \item \code{@title}
#'   \item \code{@description}
#'   \item \code{@param}
#'   \item \code{@details}
#'   \item \code{@return}
#'   \item \code{@author}
#'   \item \code{@examples}
#'   \item \code{@references}
#' }
#'
#' @return 
#' A new \code{.R} file containing a roxygen2 documentation template is 
#' written to the specified directory. The function returns \code{NULL} 
#' invisibly.
#'
#' @author 
#' Author name supplied via \code{author_name}.
#'
#' @examples
#' # Example usage:
#' # gl.document(
#' #   func_name = myFunction,
#' #   author_name = "Your Name",
#' #   example_dataset = "gl.example",
#' #   outputDir = "R/"
#' # )
#'
#'
#' @export

gl.document <- function(func_name, 
                        author_name,
                        example_dataset = NULL,
                        outputDir){
  
  # Convert function to character name
  funcName <- as.character(substitute(func_name))
  
  # Extract parameter names
  parameters <- names(formals(func_name))
  
  # ---- Construct example call ----
  example_call <- paste0(funcName, "(", 
                         paste(parameters, collapse = ", "),
                         ")")
  
  # Build example block
  example_lines <- c(
    "#' @examples",
    "#' # Example usage:",
    if(!is.null(example_dataset)) 
      paste0("#' data(", example_dataset, ")"),
    paste0("#' ", example_call)
  )
  
  # Remove NULL if no dataset
  example_lines <- example_lines[!is.na(example_lines)]
  
  # ---- Create file connection ----
  fileConn <- file(paste0(outputDir, funcName, ".R"), "wt")
  
  # ---- Write header fields ----
  writeLines(paste0("#' @name ", funcName), fileConn)
  writeLines(paste0("#' @title Title for ", funcName), fileConn)
  writeLines(paste0("#' @description ", funcName, " does:"), fileConn)
  writeLines("#'", fileConn)
  
  # ---- Write parameters ----
  for(params in parameters){
    writeLines(paste0("#' @param ", params, " Insert description."), fileConn)
  }
  
  writeLines("#'", fileConn)
  
  # ---- Write details ----
  writeLines(paste0("#' @details"), fileConn)
  writeLines("#' Detailed description goes here.", fileConn)
  writeLines("#'", fileConn)
  
  # ---- Write return ----
  writeLines(paste0("#' @return ", funcName, " returns:"), fileConn)
  writeLines("#'", fileConn)
  
  # ---- Write author ----
  writeLines(paste0("#' @author ", author_name), fileConn)
  writeLines("#'", fileConn)
  
  # ---- Write examples ----
  writeLines(example_lines, fileConn)
  writeLines("#'", fileConn)
  
  # ---- Write references ----
  writeLines("#' @references", fileConn)
  writeLines("#' Patterson, J. (2005). Maximum ride. New York: Little, Brown.", fileConn)
  writeLines("#'", fileConn)
  
  writeLines("#' @export", fileConn)
  
  close(fileConn)
  
  invisible(NULL)
}

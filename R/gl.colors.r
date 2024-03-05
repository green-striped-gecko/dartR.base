#' @name gl.colors
#' @title Returns a list of colors for use in plots
#' @family graphics
#' 
#' @description
#' Creates a vector of colors in hex (e.g. "#3B9AB2" "#78B7C5") based on user selected category 
#' (parameter type).
#'  \itemize{
#'  \item "2" [two colors]
#'  \item "2c" [two colors contrast]
#'  \item "3" [three colors]
#'  \item  4" [four colors]
#'  \item "pal" [need to be specify the palette type and the number of colors]
#'  }
#'  
#' A palette of colors 
#' can be specified via "div" [divergent], "dis" [discrete], "con" [convergent], "vir" [viridis]. 
#' Be aware a palette needs the number of colors specified as well. It returns a function 
#' and therefore the number of colors needs to be a part of the function call. 
#' Check the examples to see how this works. 
#' 
#' @param type Type of color (2, 3 or 4 colors, or palette, see description) [default 2]. 
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#' 
#' @author Custodian: Bernd Gruber -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' 
#' @examples
#' gl.colors(2)
#' gl.colors("2")
#' gl.colors("2c")
#' #five discrete colors
#' gl.colors(type="dis")(5)
#' #seven divergent colors
#' gl.colors("div")(7)
#' @export
#' @return returns colors as a vector 

#' 
gl.colors <- function(
          type=2,
          verbose=NULL) {
  
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)

  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "v.2023.3",
                   verbose = verbose)
 
  # FUNCTION SPECIFIC CHECKS
  if (is.numeric(type))
    if (type<2 | type>4)
    {
      cat(error("No valid color option specified. Check ?gl.colors\n"))
      stop(-1)
    }
  
  check <- match(type, c("2","2c","3","4","structure","dis", "div", "vir","con"))
  
  if (is.na(check)) {
    cat(error("No valid color option specified. Check ?gl.colors\n"))
    stop(-1)
  }
  
  # DO THE JOB
  
  if(verbose >=2){cat(report("Selected color type",type,"\n"))}
  
  type <- as.character(type)
  cols <- NA
  # SET PLOTS COLORS
  if (type=="2")  cols <- c("#3B9AB2", "#78B7C5")
  if (type=="2c")  cols <- c("deeppink", "chartreuse3")
  if (type=="3")  cols <-  c("#3B9AB2", "deeppink", "lemonchiffon")
  if (type=="4")  cols <-  c("lemonchiffon", "deeppink", "dodgerblue", "chartreuse3")
  if (type=="structure")  cols<-  c( "#E4E1E3", "#F6222E", "#FE00FA", "#16FF32","#3283FE", "#FEAF16", "#B00068", "#1CFFCE", "#90AD1C","#2ED9FF", "#DEA0FD", "#AA0DFE", "#F8A19F", "#325A9B","#C4451C", "#1C8356", "#85660D", "#B10DA1", "#FBE426","#1CBE4F", "#FA0087", "#FC1CBF", "#F7E1A0", "#C075A6","#782AB6", "#AAF400", "#BDCDFF", "#822E1C", "#B5EFB5","#7ED7D1", "#1C7F93", "#D85FF7", "#683B79", "#66B0FF","#3B00FB")
  
    if (type=="div") cols <- colorRampPalette(c("#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00"))
    
    if (type=="dis") cols <- function(n) {
      hues <-seq(15, 375, length = n + 1)
      return(hcl(h = hues, l = 65, c = 100)[1:n])
    }
    if (type=="con") {
      cool <-
        rainbow(50, start = rgb2hsv(col2rgb("cyan"))[1], end = rgb2hsv(col2rgb("blue"))[1])
      warm <-
        rainbow(50, start = rgb2hsv(col2rgb("red"))[1], end = rgb2hsv(col2rgb("yellow"))[1])
      dummy <- c(rev(cool), rev(warm))
      cols <- colorRampPalette(dummy)
    }
    if (type=="vir")  {
      cols <- colorRampPalette(
        c(
          "#440154FF",
          "#482173FF",
          "#433E85FF",
          "#38598CFF",
          "#2D708EFF",
          "#25858EFF",
          "#1E9B8AFF",
          "#2BB07FFF",
          "#51C56AFF",
          "#85D54AFF",
          "#C2DF23FF",
          "#FDE725FF"
        ))
    }
  
  # FLAG SCRIPT END ---------------
  
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  
return(invisible(cols))
}
#' @title Compute Derivative function
#' @description This function finds partial derivatives of the regression equation with respect a given variable in the model.
#' @param varcov the variance-covariance matrix from any model of interest
#' @param diffTerm  the variable in the model with respect to which to partial derivatives of the regression equation are computed
#' @export
#' @import foreign
#' @import graphics
#' @import stats
#' @return partial derivative of the regression equation
#' @examples
#'
#' dat <- read.csv(system.file("example.csv", package="Interactions"))
#' reg <- lm(polity2 ~ Exclusivity_Index * time_since_initial_const, data = dat)
#' compute.derivative(vcov(reg), "Exclusivity_Index")

compute.derivative<-function(varcov, diffTerm){

  eq_terms <- colnames(varcov)
  eq_terms_bs <- paste("b",(1:length(eq_terms))-1, sep="")
  colnames(varcov) <- rownames(varcov) <- eq_terms_bs

  if(grepl("Intercept", eq_terms[1])){
    eq_m <- paste(paste("b",1:(length(eq_terms)-1),"*",eq_terms[-1], sep=""), collapse="+") #reconstruct equation
  }else{
    eq_m <- paste(paste("b",0:(length(eq_terms)-1),"*",eq_terms, sep=""), collapse="+")
  }

  eq_m <- gsub(":", "*", eq_m); eq_m <- gsub("\\(", "__", eq_m); eq_m <- gsub("\\)", "___", eq_m)
  differ_eq <- deparse(D(parse(text=eq_m), diffTerm))
  differ_eq <- gsub("__", "(", gsub("___", ")", paste(differ_eq, collapse="")))
  differ_eq <- gsub("\\s+", " ", paste(differ_eq, collapse=""))

  return(differ_eq)}

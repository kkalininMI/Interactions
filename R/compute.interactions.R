#' @title Compute Derivative function
#' @description This function helps to plot marginal effects of the linear model.
#' @param formula the model formula
#' @param summarytable summary table with beta coefficients
#' @param varcov the variance-covariance matrix from any model of interest
#' @param data data frame
#' @param diffTerm  the variable in the model with respect to which to partial derivatives of the regression equation a computed
#' @param Scalor the list of values for the variables used in the derivative of the regression model
#' @param ci 95 percent confidence interval (by default)
#' @param df degrees of freedom
#' @param graph TRUE/FALSE R plot is drawn, FALSE, by default
#' @param drawgg TRUE/FALSE ggplot is drawn, FALSE, by default
#' @param xlab_text x-axis label
#' @param ylab_text y-axis label
#' @param main_text add title to a plot
#' @param mean_lty line type for the average line
#' @param mean_lwd line width for the average line
#' @param mean_col line color for the average line
#' @param ci_lty line type for the confidence interval
#' @param ci_lwd line width for the average line
#' @param ci_col line color for the average line
#' @param yint_lty line type for the intercept line
#' @param yint_lwd line width for the intercept line
#' @param yint_col line color for the intercept line
#' @export
#' @import foreign
#' @import ggplot2
#' @import graphics
#' @import stats
#' @return partial derivative of the regression equation
#' @examples
#'
#' library(Amelia)
#' library(mitools)
#' library(plm)
#' library(ggplot2)
#'
#' compute_cluster<-function(dat, formula, model, effect){
#'            fit <- lapply(dat, function(x){
#'            pdta <- pdata.frame(x, index=c("COUNTRY", "YEAR"));
#'            if(model!='between'|model!='random'){
#'            plm(formula, data = pdta, model = model, effect = effect)
#'            }else{
#'            plm(formula, data = pdta, model = model)
#'            }
#'            })
#'
#'            vcovHC_matrix <- lapply(fit, function(y) vcovBK(y, cluster="group"))
#'            betas <- MIextract(fit, fun = coef)
#'            r.sq <- round(mean(unlist(MIextract(fit, fun = r.squared))), 2)
#'            sum.m  <- summary(MIcombine(betas, vcovHC_matrix))
#'            vars <- vcovHC_matrix
#'            multi.dat<- MIextract(fit, fun=model.frame)
#'            combine.dat <-Reduce(`+`, multi.dat) / length(multi.dat)
#'            sum.m <- summary(MIcombine(betas, vars))
#'            vcovm <-MIcombine(betas, vcovHC_matrix)$variance
#'            res <- list (betas=betas, vars=vars, sum.m=sum.m,
#'            vcovm=vcovm, combine.dat=combine.dat, r.sq=r.sq)
#'
#'            return(res)}
#'
#' dat<-read.csv(system.file("example.csv", package="Interactions"))
#'
#' #Imputation of the panel data
#' a.out <- amelia(dat, m = 3, ts = "YEAR", cs = "COUNTRY",
#'                ords = "polity2", polytime=1, lags="economics1")
#'
#' mydata <- imputationList(list(a.out[[1]]$imp1, a.out[[1]]$imp2, a.out[[1]]$imp3))
#'
#' #Run plm on imputed data
#' mformula <- as.formula("polity2 ~ Exclusivity_Index * time_since_initial_const +
#'                                  Exclusivity_Index * lag(economics1) + al_ethnic + al_language + al_religion")
#' run_mod <- compute_cluster(dat = mydata$imputations, formula = mformula, model = "random", effect=NULL)
#'
#' #Find partial derivative of regression equation
#' compute.derivative(run_mod$vcovm, "Exclusivity_Index")
#' # "b1 + b7 * time_since_initial_const + b8 * lag(economics1)"
#' # For every variable that appears in the derivative, the user must specify the range of values for that variable.
#' # Only one variable can have a range of values more than 1 (i.e. running variable); other variables need to be fixed at a single value, such as the mean.
#' # If a partial derivative contains interaction terms that involve variables with previously defined scales, those terms can be ignored.
#' # "b1 + b7 * time_since_initial_const + b8 * lag(economics1)"
#' #           ^      seq(1,80,1)      ^        ^     mean    ^
#' #Subsequently, all the defined scales or rulers can be merged into a list.
#' mod.scalor = list(seq(0,80,.01), 8.64)
#'
#' #Compute and draw marginal effects
#'
#' compute.interactions(formula = mformula,
#'                     summarytable = run_mod$sum.m,
#'                     varcov = run_mod$vcovm,
#'                     data = run_mod$combine.dat,
#'                     diffTerm = "Exclusivity_Index",
#'                     xlab_text = "Time",
#'                     ylab_text = expression(paste(partialdiff,"Polity"," / ",partialdiff, EI[general])),
#'                     Scalor = mod.scalor,
#'                     ci=0.95,
#'                     df = 1000,
#'                     graph=FALSE,
#'                     drawgg=TRUE,
#'                     ci_col = "blue")

compute.interactions<-function(formula, summarytable, varcov,
                               data = NULL, diffTerm, Scalor,
                               ci = 0.95, df = 1000, graph = FALSE, drawgg = FALSE,
                               xlab_text = NULL, ylab_text = NULL, main_text = NULL,
                               mean_lty = 1, mean_lwd = 1, mean_col = "black",
                               ci_lty = 1, ci_lwd = 0.5, ci_col = "black",
                               yint_lty = 2, yint_lwd = 1, yint_col = "black"){

  #Build equation out of variance-covariance matrix
  eq_terms <- colnames(varcov)
  eq_terms <- gsub("\\(", "__", eq_terms); eq_terms <- gsub("\\)", "___", eq_terms)

  eq_terms_bs <- paste("b",(1:length(eq_terms))-1, sep="")
  colnames(varcov) <- rownames(varcov) <- eq_terms_bs

  if(grepl("Intercept", eq_terms[1])){
    eq_m <- paste(paste("b",1:(length(eq_terms)-1),"*",eq_terms[-1], sep=""), collapse="+") #reconstruct equation
  }else{
    eq_m <- paste(paste("b",0:(length(eq_terms)-1),"*",eq_terms, sep=""), collapse="+")
  }

  eq_m <- gsub(":", "*", eq_m)
  differ_eq <- deparse(D(parse(text=eq_m), diffTerm))
  differ_eq <- gsub("\\s+", " ", paste(differ_eq, collapse=""))

  #Difference section
  betas <- unlist(regmatches(
    differ_eq, gregexpr("b\\d+((?=\\s+\\*)|(?=\\s+\\+))", differ_eq, perl=TRUE)))

  diffvar <- gsub("(^\\* )|( \\+$)", "", regmatches(differ_eq, gregexpr("(\\*.*?\\+)|(\\*.*?$)", differ_eq))[[1]])
  v_coef <- summarytable[,1]; names(v_coef) <- eq_terms
  betasEQ <- paste("b",seq(0, length(v_coef)-1), sep="")

  #Create new variables
  for(i in 1:length(betasEQ)) assign(betasEQ[i],v_coef[i])

  for(i in 1:length(Scalor)){
    assign(paste("Scalor_Variable",i, sep=""), Scalor[[i]])
    differ_eq<-gsub(diffvar[i], paste("Scalor_Variable",i, sep=""), differ_eq)
  }

  #Estimate mean and SD
  mean.est <- eval(parse(text = differ_eq))

  VarCovBetas <- expand.grid(betas,betas)
  VarCovBetas$Var1F <- VarCovBetas$Var1T <-
    eq_terms[as.numeric(gsub("b", "", VarCovBetas$Var1))+1]
  VarCovBetas$Var2F <- VarCovBetas$Var2T <-
    eq_terms[as.numeric(gsub("b", "", VarCovBetas$Var2))+1]

  for(i in 1:length(diffvar)){
    VarCovBetas$Var1F <-
      sub(diffvar[i],paste("Scalor_Variable",i, sep=""), VarCovBetas$Var1F)
    VarCovBetas$Var2F <-
      sub(diffvar[i],paste("Scalor_Variable",i, sep=""), VarCovBetas$Var2F)
    VarCovBetas$Var1F <-
      sub(diffTerm,"Intercept",VarCovBetas$Var1F)
    VarCovBetas$Var2F <-
      sub(diffTerm,"Intercept",VarCovBetas$Var2F)
  }

  VarCovBetas$VarEst <- c(varcov[rownames(varcov) %in% VarCovBetas$Var1,
                                 colnames(varcov) %in% VarCovBetas$Var2])

  d1 <- data.frame(do.call("rbind", strsplit(VarCovBetas$Var1F, "\\:")))
  d2 <- data.frame(do.call("rbind", strsplit(VarCovBetas$Var2F, "\\:")))

  VarCovBetas[c(paste("Var1S", 1:length(d1), sep=""))] <-
    data.frame(do.call("rbind", strsplit(VarCovBetas$Var1F, "\\:")))
  VarCovBetas[paste("Var2S", 1:length(d1), sep="")] <-
    data.frame(do.call("rbind", strsplit(VarCovBetas$Var2F, "\\:")))

  VarCovBetas<-VarCovBetas[!duplicated(VarCovBetas$VarEst),]
  VarCovBetas[is.na(VarCovBetas)] <- "Intercept"

  acc_var <- rep(0, max(unlist(lapply(Scalor,FUN=length))))
  addit <- 0; Intercept <- 1

  for(i in 1:nrow(VarCovBetas)){
    multiplier <- paste(
      paste(lapply(1:length(d1), function(x)
        paste("get(VarCovBetas$Var1S", x, "[i])", sep="")), collapse=" * "),
      paste(lapply(1:length(d2), function(x)
        paste("get(VarCovBetas$Var2S", x , "[i])", sep="")), collapse=" * "), sep = " * ")

    if(VarCovBetas$Var1[i] == VarCovBetas$Var2[i]){
      addit = VarCovBetas$VarEst[i] * eval(parse(text = multiplier))}
    if(VarCovBetas$Var1[i] != VarCovBetas$Var2[i]){
      addit = 2 * VarCovBetas$VarEst[i] * sqrt(eval(parse(text = multiplier)))
    }

    acc_var <- acc_var + addit; addit=0
  }

  upperbound <- mean.est + abs(qt((1-ci)/2, df = df)) * sqrt(acc_var)
  lowerbound <- mean.est - abs(qt((1-ci)/2, df = df)) * sqrt(acc_var)

  max_id_scalor <- which.max(unlist(lapply(Scalor,FUN=length)))
  min_id_scalor <- which.min(unlist(lapply(Scalor,FUN=length)))

  if(graph==TRUE){

    if(is.null(ylab_text)) ylab_text = paste("d", as.character(formula)[2],"/","d",diffTerm, sep="")
    if(is.null(xlab_text)) xlab_text = diffvar[max_id_scalor]
    if(is.null(main_text)) main_text = "Marginal Effect Plot"

    if (length(diffvar) > 1){

      title_part <- paste(apply(cbind(
        diffvar[-c(max_id_scalor, which(grepl("\\*", diffvar)))],  #remove interaction term
        unlist(Scalor[-max_id_scalor])), 1, function(x) paste (x, collapse="=")), collapse=", ")


      plot(Scalor[[max_id_scalor]], mean.est,
           ylim = c(min(c(upperbound,lowerbound)), max(c(upperbound,lowerbound))),
           xlim = c(min(Scalor[[max_id_scalor]]), max(Scalor[[max_id_scalor]])),
           ylab = ylab_text,
           xlab = xlab_text,
           lty = 1, type = "l",
           main = paste("Marginal Effect (",
                        gsub("__", "(", gsub("___", ")", paste(title_part, collapse=""))),")", sep=""))

      abline(h=0, lty = yint_lty, lwd = yint_lwd, col = yint_col)

    }else{

      plot(Scalor[[max_id_scalor]], mean.est,
           ylim=c(min(c(upperbound,lowerbound)),max(c(upperbound,lowerbound))),
           xlim = c(min(Scalor[[max_id_scalor]]), max(Scalor[[max_id_scalor]])),
           ylab = ylab_text,
           xlab = xlab_text,
           lty = mean_lty,
           lwd = mean_lwd,
           col = mean_col,
           type = "l",
           main = main_text)

      abline(h=0, lty = yint_lty, lwd = yint_lwd, col = yint_col)
    }

    par(new=TRUE)

    plot(Scalor[[max_id_scalor]],upperbound,
         ylim = c(min(c(upperbound,lowerbound)), max(c(upperbound,lowerbound))),
         xlim = c(min(Scalor[[max_id_scalor]]), max(Scalor[[max_id_scalor]])),
         ylab = "", xlab = "", lty = ci_lty, lwd = ci_lwd, type = "l", col = ci_col)

    par(new=TRUE)

    plot(Scalor[[max_id_scalor]],lowerbound,
         ylim = c(min(c(upperbound,lowerbound)),max(c(upperbound,lowerbound))),
         xlim = c(min(Scalor[[max_id_scalor]]), max(Scalor[[max_id_scalor]])),
         ylab = "", xlab = "", lty = ci_lty, lwd = ci_lwd, type = "l", col = ci_col)

    if (!is.null(data)) rug(data[, colnames(data) %in% diffvar[max_id_scalor]])
  }

  if (drawgg==TRUE){

    if (!is.null(data)){
      data$rug_var <- data[, colnames(data) %in% diffvar[max_id_scalor]]
      data$rug_var <- ifelse(data$rug_var > max(Scalor[[max_id_scalor]]),
                             max(Scalor[[max_id_scalor]]), data$rug_var)
    }

    if(length(Scalor)==1 & length(Scalor[[1]])==2){
      ggraph <- ggplot(data=NULL,aes(x = Scalor[[max_id_scalor]], group = Scalor[[max_id_scalor]]))+
        geom_boxplot(aes(lower = lowerbound, upper = upperbound,
                         middle = mean.est,
                         ymin = mean.est - abs(qt((1 - ci)/2, df = df)) * sqrt(acc_var),
                         ymax = mean.est + abs(qt((1 - ci)/2, df = df)) * sqrt(acc_var)), stat="identity") +
        labs(x = xlab_text, y = ylab_text, title = "")

    }else{

      ggraph <- ggplot() +
        geom_line(aes(Scalor[[max_id_scalor]], mean.est), size = mean_lwd,
                  linetype = mean_lty, color = mean_col) +
        geom_line(aes(Scalor[[max_id_scalor]], lowerbound), size = ci_lwd,
                  linetype = ci_lty, color = ci_col) +
        geom_line(aes(Scalor[[max_id_scalor]], upperbound), size = ci_lwd,
                  linetype = ci_lty, color = ci_col) +
        geom_hline(yintercept=0, linetype=yint_lty,
                   size = yint_lwd, color = yint_col) +
        geom_ribbon(aes(x = Scalor[[max_id_scalor]],
                        ymin = lowerbound, ymax = upperbound), alpha=0.2, color = ci_col, fill = ci_col) +
        labs(x = xlab_text, y = ylab_text, title = "") +
        geom_rug(data = data, aes(x = rug_var), color = "black", sides = "b")

    }

  }else{

    ggraph <- NULL

  }

  results = list (mean.est = mean.est,
                  upperbound = upperbound,
                  lowerbound = lowerbound,
                  ggraph = ggraph)

  return(results)
}

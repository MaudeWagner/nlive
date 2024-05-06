#' Automated Estimation of Sigmoidal and Piecewise Linear Mixed Models
#'
#'
#' The nlive() function allows to fit a Sigmoidal Mixed Model with 4 parameters, a Piecewise Linear Mixed Model
#' with abrupt change, or a Piecewise Linear Mixed Model with a smooth polynomial transition in the context of
#' longitudinal Gaussian outcomes. This function was designed to be
#' intuitive enough to the less sophisticated users, while using recent developments such as the
#' stochastic approximation expectation-maximization (SAEM) algorithm for efficient estimation. It was
#' designed to optimize the initial values of the main parameters and help
#' interpretation of the output by providing different features such as annotated outputs and graphs.
#'
#'
#' CAUTIONS REGARDING THE USE OF THE FUNCTION
#'
#' traj.marg: if "TRUE", this argument automatically plots the estimated marginal trajectories of the longitudinal outcome
#' for the most common profile of covariates, if any (i.e., ref "1" for binary variables and mean values for continuous variables).
#' Thus, users must ensure that continuous variables are centered on the mean.
#'
#'
#' @param model indicator of the model to fit (1=Sigmoidal Mixed Model, 2=Piecewise Mixed Model with abrupt change, 3=Piecewise Mixed Model with smooth transition)
#' @param dataset data frame containing the variables ID, outcome, time, predictor.all, and predictor.par1 to predictor.par4.
#' @param ID name of the variable representing the grouping structure specified with " (e.g., "ID" representing the unique identifier of participants).
#' @param outcome name of the time-varying variable representing the longitudinal outcome specified with " (e.g., "outcome").
#' @param time name of the variable representing the timescale specified with " (e.g., "time"), which can be negative or positive.
#' @param predictor.all optional vector indicating the name of the variable(s) that the four main parameters of the model of interest will be adjusted to (e.g. predictor.all=c("X1","X2")). Default to NULL.
#' @param predictor.par1 optional vector indicating the name of the variable(s) that the first main parameter of the model of interest will be adjusted to (e.g. predictor.all=c("X1","X2")). For model 1, the first parameter = last level. For models 2 and 3, first parameter = intercept. Default to NULL.
#' @param predictor.par2 optional vector indicating the name of the variable(s) that the second main parameter of the model of interest will be adjusted to (e.g. predictor.all=c("X1","X2")). For model 1, the second parameter = initial level. For models 2 and 3, second parameter = slope before the change-point. Default to NULL.
#' @param predictor.par3 optional vector indicating the name of the variable(s) that the third main parameter of the model of interest will be adjusted to (e.g. predictor.all=c("X1","X2"). For model 1, the third parameter = midpoint. For models 2 and 3, third parameter = slope after the change-point. Default to NULL.
#' @param predictor.par4 optional vector indicating the name of the variable(s) that the fourth main parameter of the model of interest will be adjusted to (e.g. predictor.all=c("X1","X2")). For model 1, the fourth parameter is the Hill slope. For models 2 and 3, the fourth parameter is the changepoint. Default to NULL.
#' @param start optional vector to override the specification of the four initial values for the main parameters. For model 1, the values must be included in the following order: last level, initial level, midpoint, Hill slope. For models 2 and 3, the values must be included in the following order: intercept, slope before the changepoint, slope after the changepoint, changepoint. Default to NULL.
#' @param traj.marg optional logical indicating if the marginal estimated trajectory should be plotted for the most common profile of covariates, if any. Default to FALSE.
#' @param traj.marg.group  optional name of the grouping variable listed in one of the predictor arguments to plot and contrast the estimated marginal trajectories between two specific groups, specified with " (e.g., traj.marg.group="X1"). If the variable is binary, the trajectories are contrasted between the two groups of interest. If the variable is continuous, the 10th and 90th percentile values will automatically be considered. The default value is NULL.
#' @param plot.xlabel optional text for the title of the x-axis of all plots
#' @param plot.ylabel optional text for the title of the y-axis of all plots
#' @param spag.plot.title optional text for the title of the spaghetti plot
#' @param traj.marg.title optional text for the title of the marginal estimated trajectory
#' @param traj.marg.group.title optional text for the title of the marginal estimated trajectories contrasted between groups
#' @param traj.marg.group.val optional vector that can be used when _traj.marg.group_ receives a quantitative variable and that allows to manually specify two percentile values to be considered for contrasting the traj.marg.group. The two values must be between 0 and 1 (e.g., traj.marg.group.val=c(0.2,0.8); for percentiles 20th and 80th). Default to 10th and 90th percentiles (i.e., traj.marg.group.val=c(0.1,0.9)).
#'
#' @return An object of class SaemixObject (from the existing _saemix_ R package) containing the results of
#' the fit of the data by the non-linear mixed _model_ of interest. The _nlive_ function automatically provides
#' (i) a spaghetti plot of the observed outcome for 70 randomly selected statistical units in the dataset and
#' (ii) the standard _saemix_ output, including the fixed effects estimates, the variance of random effects,
#' and Likelihood of the fitted _model_. The outputs are printed on the terminal and the numerical and
#' graphical outputs are stored in a directory.
#'
#' @author Maude Wagner, Ana W. Capuano, Emmanuelle Comets
#'
#' \email{maude_wagner@@rush.edu}
#'
#' @seealso
#'
#' \code{\link{nlive.smm}}, \code{\link{nlive.pmma}}, \code{\link{nlive.pmms}}
#'
#' @references
#'
#' Capuano AW, Wagner M. nlive: an R package to facilitate the application of the sigmoidal and random changepoint mixed models. BMC Medical Research Methodology. 2023;23(1):257.
#' van den Hout A, Muniz-Terrera G, Matthews F. Smooth random change point models. Statistics in Medicine. 2011;30(6):599-610.
#' Comets E, Lavenu A, Lavielle MM. Parameter estimation in nonlinear mixed effect models using saemix, an R implementation of the SAEM algorithm. Journal of Statistical Software. 2017;80(3):1-41.
#'
#' @examples
#'
#' #### Fitting a sigmoidal mixed model - with no covariate
#' \dontrun{
#' head(dataCog)
#' requireNamespace('nlraa')
#' smm.fit = nlive(model=1, dataset=dataCog, ID="ID", outcome="cognition", time="time")
#' }
#'
#' #### plot(smm.fit): diagnostic plots to assess the goodness-of-fit of smm.fit
#' #### psi(smm.fit): estimates of individual parameters
#'
#' #### Fitting a piecewise mixed model with abrupt change - with no covariate
#' \dontrun{
#' pmm.abrupt.fit = nlive(model=2, dataset=dataCog, ID="ID", outcome="cognition", time="time")
#' }
#' #### plot(pmm.abrupt.fit): diagnostic plots to assess the goodness-of-fit of pmm.abrupt.fit
#' #### psi(pmm.abrupt.fit): estimates of individual parameters
#'
#'
#' #### Fitting a piecewise mixed model with smooth change - with all parameters
#' #### adjusted for ageDeath90. Here, the nlive() function will also provide a
#' #### plot of the estimated marginal trajectory in the whole study sample.
#' \dontrun{
#' pmm.smooth.fit = nlive(model=3, dataset=dataCog, ID="ID", outcome="cognition", time="time",
#' predictor.all = c("ageDeath90"), traj.marg = TRUE)
#' }
#' #### plot(pmm.smooth.fit): diagnostic plots to assess the goodness-of-fit of the model
#' #### psi(pmm.smooth.fit): estimates of individual parameters
#'
#'
#' @import Rmisc
#' @import sitar
#' @import nlraa
#' @import knitr
#' @importFrom nlraa SSlogis5
#' @importFrom dplyr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr row_number
#' @importFrom dplyr sample_n
#' @importFrom dplyr filter
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 guides
#' @importFrom graphics points
#' @importFrom graphics legend
#' @importFrom graphics plot
#' @importFrom lcmm hlme
#' @importFrom saemix saemix
#' @importFrom saemix saemixModel
#' @importFrom saemix saemixData
#' @importFrom saemix saemixControl
#' @importFrom saemix normcdf
#' @importFrom saemix psi
#' @importFrom saemix plot
#' @importFrom stats quantile
#' @importFrom stats na.omit
#' @importFrom sqldf sqldf
#' @importFrom fastDummies dummy_cols
#'
#' @export

nlive <- function(model, dataset, ID, outcome, time, predictor.all = NULL,
                  predictor.par1 = NULL,
                  predictor.par2 = NULL,
                  predictor.par3 = NULL,
                  predictor.par4 = NULL,
                  start = NULL,
                  plot.xlabel = NULL,
                  plot.ylabel = NULL,
                  traj.marg = FALSE,
                  traj.marg.group = NULL,
                  spag.plot.title = NULL,
                  traj.marg.title = NULL,
                  traj.marg.group.title = NULL,
                  traj.marg.group.val = NULL){


  if (missing(model))
    stop("The argument model must be specified in any model")
  if (missing(dataset))
    stop("The argument dataset should be specified and defined as a data.frame")
  if (nrow(dataset) == 0)
    stop("Data should not be empty")
  if (missing(ID))
    stop("The argument ID must be specified in any model")
  if (!is.numeric(dataset[[ID]]))
    stop("The argument ID must be numeric")
  if (missing(outcome))
    stop("The argument outcome must be specified in any model")
  if (!is.numeric(dataset[[outcome]]))
    stop("The argument outcome must be numeric")
  if (missing(time))
    stop("The argument time must be specified in any model")
  if (!is.numeric(dataset[[time]]))
    stop("The argument time must be numeric")
  if (!is.null(predictor.all) & !is.vector(predictor.all))
    stop("The argument predictor.all should receive a vector indicating the name of the variable(s) that the four main parameters of the model will be adjusted to.")
  if (!is.null(predictor.par1) & !is.vector(predictor.par1))
    stop("The argument var.first.level should receive a vector indicating the name of the variable(s) that the first main parameters of the model will be adjusted to.")
  if (!is.null(predictor.par2) & !is.vector(predictor.par2))
    stop("The argument var.last.level should receive a vector indicating the name of the variable(s) that the second main parameters of the model will be adjusted to.")
  if (!is.null(predictor.par3) & !is.vector(predictor.par3))
    stop("The argument var.midpoint should receive a vector indicating the name of the variable(s) that the third main parameters of the model will be adjusted to.")
  if (!is.null(predictor.par4) & !is.vector(predictor.par4))
    stop("The argument var.Hslope should receive a vector indicating the name of the variable(s) that the fourth main parameters of the model will be adjusted to.")
  if (!is.null(start) & !is.vector(start))
    stop("The argument start must be a vector specifying the four initial values for the four main parameters in the following order: c(last level, initial level, midpoint, Hill slope).")

  if (!is.null(traj.marg.group) & is.null(predictor.all) & is.null(predictor.par1) & is.null(predictor.par2) & is.null(predictor.par3) & is.null(predictor.par4))
    stop("The covariate in traj.marg.group should also be in at least one of these arguments: predictor.all, predictor.par1, predictor.par2, predictor.par3, predictor.par4.")
  if (!is.null(traj.marg.group.val) & is.null(predictor.all) & is.null(predictor.par1) & is.null(predictor.par2) & is.null(predictor.par3) & is.null(predictor.par4))
    stop("The values in traj.marg.group.val are not used because no covariate specified in the argument traj.marg.group.")



  ## dataset
  dataset$ID      = na.omit(dataset[,ID])
  dataset$outcome = na.omit(dataset[,outcome])
  dataset$time    = na.omit(dataset[,time])


  ##########################################
  ## accounting for categorical variables ##
  ##########################################

  ## predictor.all
  data_subset     = subset(dataset, select = unique(unlist(as.list(c(predictor.all)))))
  factor_cat_sup2 = as.list(sapply(data_subset, function(col) is.factor(col) && nlevels(col) > 2))
  names_col       = names(which(factor_cat_sup2 == T))
  ##  Dummy coding when nlevels(var.factor) > 2
  if (length(names_col) == 0) {predictor.all2 = predictor.all } else if (length(names_col) > 0) {
    data_subset2    = dummy_cols(data_subset, select_columns = names(which(factor_cat_sup2 == T)), remove_first_dummy = T, ignore_na = T, remove_selected_columns = T)
    predictor.all2  = as.vector(colnames(data_subset2))
  }

  ## predictor.par1
  data_subset     = subset(dataset, select = unique(unlist(as.list(c(predictor.par1)))))
  factor_cat_sup2 = as.list(sapply(data_subset, function(col) is.factor(col) && nlevels(col) > 2))
  names_col       = names(which(factor_cat_sup2 == T))
  ##  Dummy coding when nlevels(var.factor) > 2
  if (length(names_col) == 0) { predictor.par1_2 = predictor.par1 } else if (length(names_col) > 0) {
    data_subset2    = dummy_cols(data_subset, select_columns = names(which(factor_cat_sup2 == T)), remove_first_dummy = T, ignore_na = T, remove_selected_columns = T)
    predictor.par1_2 = as.vector(colnames(data_subset2))
  }

  ## predictor.par2
  data_subset     = subset(dataset, select = unique(unlist(as.list(c(predictor.par2)))))
  factor_cat_sup2 = as.list(sapply(data_subset, function(col) is.factor(col) && nlevels(col) > 2))
  names_col       = names(which(factor_cat_sup2 == T))
  ##  Dummy coding when nlevels(var.factor) > 2
  if (length(names_col) == 0) { predictor.par2_2 = predictor.par2 } else if (length(names_col) > 0) {
    data_subset2  = dummy_cols(data_subset, select_columns = names(which(factor_cat_sup2 == T)), remove_first_dummy = T, ignore_na = T, remove_selected_columns = T)
    predictor.par2_2  = as.vector(colnames(data_subset2))
  }

  ## predictor.par3
  data_subset     = subset(dataset, select = unique(unlist(as.list(c(predictor.par3)))))
  factor_cat_sup2 = as.list(sapply(data_subset, function(col) is.factor(col) && nlevels(col) > 2))
  names_col       = names(which(factor_cat_sup2 == T))
  ##  Dummy coding when nlevels(var.factor) > 2
  if (length(names_col) == 0) { predictor.par3_2 = predictor.par3 } else if (length(names_col) > 0) {
    data_subset2  = dummy_cols(data_subset, select_columns = names(which(factor_cat_sup2 == T)), remove_first_dummy = T, ignore_na = T, remove_selected_columns = T)
    predictor.par3_2  = as.vector(colnames(data_subset2))
  }

  ## predictor.par4
  data_subset     = subset(dataset, select = unique(unlist(as.list(c(predictor.par4)))))
  factor_cat_sup2 = as.list(sapply(data_subset, function(col) is.factor(col) && nlevels(col) > 2))
  names_col       = names(which(factor_cat_sup2 == T))
  ##  Dummy coding when nlevels(var.factor) > 2
  if (length(names_col) == 0) { predictor.par4_2 = predictor.par4 } else if (length(names_col) > 0) {
    data_subset2     = dummy_cols(data_subset, select_columns = names(which(factor_cat_sup2 == T)), remove_first_dummy = T, ignore_na = T, remove_selected_columns = T)
    predictor.par4_2 = as.vector(colnames(data_subset2))
  }

  ## definition of nb_cov & predictors & new data_subset with dummy variables (if any)
  x          = unique(unlist(as.list(c(predictor.all2, predictor.par1_2, predictor.par2_2, predictor.par3_2, predictor.par4_2))))
  nb_cov     = length(x)
  predictors = unique(x)
  ##  New dataset with dummy var (if any)
  data_subset     = subset(dataset, select = unique(unlist(as.list(c(predictor.all, predictor.par1, predictor.par2, predictor.par3, predictor.par4)))))
  factor_cols0 = as.list(sapply(data_subset, function(col) is.factor(col) && nlevels(col) > 2))
  names_col    = names(which(factor_cols0 == T))

  if (length(names_col) == 0) { dataset2 = dataset } else if (length(names_col) > 0) {
    dataset2 = dummy_cols(dataset, select_columns = names_col, remove_first_dummy = T, ignore_na = T, remove_selected_columns = T)
    ## remove NA from predictors (if any)
    if (length(predictors) > 0){
      for (i in 1:length(predictors)){
        dataset2 = subset(dataset2, is.na(dataset2[,predictors[i]]) == F)
        i = i + 1
      }
    }
  }

  ## Values used throughout the code
  first    = as.data.frame(dataset %>% group_by(ID) %>% filter(row_number(ID) == 1))
  nb_indiv = dim(first)[1]
  nb_line  = dim(dataset)[1]
  nb_col   = dim(dataset)[2]

  ## Options
  # common to all plots #
  if (is.null(plot.xlabel) == T){x.lab = time} else if (is.null(plot.xlabel) == F){x.lab = plot.xlabel}
  if (is.null(plot.ylabel) == T){y.lab = outcome} else if (is.null(plot.ylabel) == F){y.lab = plot.ylabel}

  # specific #
  if (is.null(spag.plot.title)  == T){main.spag.lab ="Spaghetti plot, random sample of 70 individuals"} else if (is.null(spag.plot.title) == F){main.spag.lab = spag.plot.title}
  if (is.null(traj.marg.title) == T){main.traj.marg = "Marginal estimated trajectories"} else if (is.null(traj.marg.title) == F){main.traj.marg = traj.marg.title}
  if (is.null(traj.marg.group.title) == T){main.traj.marg.group = "Marginal estimated trajectories between groups"} else if (is.null(traj.marg.group.title) == F){main.traj.marg.group = traj.marg.group.title}
  if (is.null(traj.marg.group.val) == T){percentiles = c(0.1,0.9)} else if (is.null(traj.marg.group.val) == F){percentiles = traj.marg.group.val}

  #############################################
  ######          SPAGHETTI PLOT         ######
  ###### 70 indivudals randomly selected ######
  ######     (automatically generated)   ######
  #############################################
  if (nb_indiv >= 70){

    ## Random selection of 70 individuals
    tempo1 = sample_n(first, 70)
    tempo2 = sqldf('SELECT * FROM dataset WHERE ID IN (SELECT ID FROM tempo1)')
    min    = min(tempo2[,time])
    max    = max(tempo2[,time])
    plot   = ggplot(tempo2, aes(time, outcome, group = ID, colour = factor(ID))) +
      geom_line(size=0.8) +
      scale_x_continuous(limits=c(min, max)) +
      guides(colour="none") +
      labs(x=x.lab, y=y.lab) +
      ggtitle(main.spag.lab)
    print(plot)

  } else if (nb_indiv < 70){

    ## Random selection for the whole study sample
    tempo1 = sample_n(first, dim(first)[1])
    tempo2 = sqldf('SELECT * FROM dataset WHERE ID IN (SELECT ID FROM tempo1)')
    min    = min(tempo2[,time])
    max    = max(tempo2[,time])
    plot   = ggplot(tempo2, aes(time, outcome, group = ID, colour = factor(ID))) +
      scale_x_continuous(limits=c(min, max)) +
      guides(colour="none") +
      geom_line(size=0.8) +
      labs(x=x.lab, y=y.lab) +
      ggtitle("Spaghetti plot, whole sample")
    print(plot)
  }


  #####################################
  ######  SIGMOIDAL MIXED MODEL  ######
  #####################################
  if (model == 1){

    dataset$time_pos = abs(dataset[,time])

    ## STARTING VALUES  ##
    if (is.null(start) == T) {

      frag   = quantile(dataset$time, probs = c(0.2,0.4,0.6,0.8))
      # interval 1
      tab_int = subset(dataset, time <= frag[1])
      lmm    = hlme(outcome ~ 1 + time, random = ~1 + time, subject=ID, data=tab_int, verbose = F)
      slope1 = mean(coef(lmm)["time"])
      # interval 2
      tab_int = subset(dataset, time >= frag[1] & time <= frag[2])
      lmm  = hlme(outcome ~ 1 + time, random = ~1 + time, subject=ID, data=tab_int, verbose = F)
      slope2 = mean(coef(lmm)["time"])
      # interval 3
      tab_int = subset(dataset, time >= frag[2] & time <= frag[3])
      lmm  = hlme(outcome ~ 1 + time, random = ~1 + time, subject=ID, data=tab_int, verbose = F)
      slope3 = mean(coef(lmm)["time"])
      # interval 4
      tab_int = subset(dataset, time >= frag[3] & time <= frag[4])
      lmm  = hlme(outcome ~ 1 + time, random = ~1 + time, subject=ID, data=tab_int, verbose = F)
      slope4 = mean(coef(lmm)["time"])
      # interval 5
      tab_int = subset(dataset, time >= frag[4])
      lmm  = hlme(outcome ~ 1 + time, random = ~1 + time, subject=ID, data=tab_int, verbose = F)
      slope5 = mean(coef(lmm)["time"])
      ###
      y = c(slope1,slope2,slope3,slope4,slope5)
      val = min(y)
      pos = frag[which(y == val)-1]

      # last level (intercept)
      tempo   = subset(dataset, time > quantile(dataset$time, probs=c(0.97)))
      last.level  = mean(tempo[, outcome])
      # first level
      tempo   = subset(dataset, time < quantile(dataset$time, probs=c(0.05)))
      first.level  = mean(tempo[, outcome])
      ##
      tempo = subset(dataset, time_pos < -pos)
      lmm   = hlme(outcome ~ 1+time_pos, random =~1+time_pos, subject=ID, data=tempo, verbose = F)
      low   = as.numeric(coef(lmm)["time_pos"])
      ##
      tempo = subset(dataset, time_pos > -pos)
      lmm   = hlme(outcome ~ 1+time_pos, random =~1+time_pos, subject=ID, data=tempo, verbose = F)
      high  = as.numeric(coef(lmm)["time_pos"])
      c(low,high)
      ##
      midpoint   = ifelse(abs(low/high)> 0.5 & abs(low/high) < 1.5, 300, 2)
      hill.slope = ifelse(abs(low/high)> 0.5 & abs(low/high) < 1.5, 1.05, 0.5)

    } else if (length(start) == 4) {
      last.level  = start[1]
      first.level = start[2]
      midpoint    = start[3]
      hill.slope  = start[4]
    }
    c(last.level, first.level, midpoint, hill.slope)

    ## CORRELATION BETWEEN n1/n2 only ##
    varCov      = diag(0, ncol=4, nrow=4)
    varCov[1,1] = varCov[2,2] = varCov[1,2] = varCov[2,1] = 1

    ## SPECIFICATION OF THE MODEL ##
    model_fct = function(psi, ID, xidep){

      t           = xidep[, 1]

      last.level  = psi[ID, 1]
      first.level = psi[ID, 2]
      midpoint    = psi[ID, 3]
      hill.slope  = psi[ID, 4]

      t2 = abs(t)
      Y_pred = SSlogis5(t2, last.level, first.level, midpoint, hill.slope, theta = 1)

      return(Y_pred)
    }


    ## SPECIFICATION OF THE MATRIX OF COVARIATES ##

    if (nb_cov == 0) {

      mat_predictor = matrix(1, ncol=4, nrow=nb_cov)

    } else if (is.null(predictor.all)==F & nb_cov == length(predictor.all)){  ##!!!!!!!

      mat_predictor = matrix(1, ncol=4, nrow=nb_cov)
      y = unique(unlist(as.list(predictor.all)))
      y = as.list(y)
      rownames(mat_predictor) = y

    } else if (is.null(predictor.all)==F & nb_cov != length(predictor.all)){
      # y related to all parameters
      mat.all = matrix(1, ncol=4, nrow=length(predictor.all))
      y = unique(unlist(as.list(predictor.all)))
      y = as.list(y)
      rownames(mat.all) = y
      # x related to other than all parameters
      mat.tempo = matrix(0, ncol=4, nrow=nb_cov-length(predictor.all))
      x = unique(unlist(as.list(c(predictor.par1, predictor.par2, predictor.par3, predictor.par4))))
      x = as.list(x)
      rownames(mat.tempo) = x

      for (i in 1:length(x)){
        if (isTRUE(x[i] %in% predictor.par1)){mat.tempo[i,1] = 1}
      }
      for (i in 1:length(x)){
        if (isTRUE(x[i] %in% predictor.par2)){mat.tempo[i,2] = 1}
      }
      for (i in 1:length(x)){
        if (isTRUE(x[i] %in% predictor.par3)){mat.tempo[i,3] = 1}
      }
      for (i in 1:length(x)){
        if (isTRUE(x[i] %in% predictor.par4)){mat.tempo[i,4] = 1}
      }
      mat_predictor = rbind(mat.all,mat.tempo)
    }
    mat_predictor

    ##
    saemix.model.cog = saemixModel(model = model_fct,
                                   psi0  = rbind(c(last.level  = last.level,
                                                   first.level = first.level,
                                                   midpoint    = midpoint,
                                                   hill.slope  = hill.slope),
                                                 matrix(0, ncol=4, nrow=1)),
                                   covariate.model  = mat_predictor,
                                   covariance.model = varCov,
                                   verbose = F)

    ##
    write.table(dataset2, file = "dataset.txt")
    dataset <- read.table('dataset.txt', header = TRUE, sep = "",dec=".")
    path = paste0(as.character(getwd()),"/dataset.txt")
    ##
    saemix.cog = saemixData(name.data       = dataset,
                            name.group      = ID,
                            name.predictors = time,
                            name.response   = outcome,
                            name.covariates = predictors,
                            verbose = F)
    ##
    saemix.options = saemixControl(map = F,
                                   nbiter.saemix = c(300,100), # nb of iterations for the exploration and smoothing phase
                                   nbiter.burn = 20,           # nb of iterations for burning
                                   nbiter.mcmc = c(5,5,5,0),   # nb of iterations in each kernel during the MCMC step
                                   ll.is = T,                  # default = T (TRUE) to estimate the log-likelihood
                                   seed = 123,
                                   displayProgress = F,        # default = T = to output the convergence plots
                                   save = F,                   # default = T = the results of the fit should be saved to a file
                                   save.graphs = F)            # default = T = to save the diagnostic and individual graphs


    ## fit
    ptm<-proc.time()
    model_SMM = saemix(saemix.model.cog, saemix.cog, saemix.options)
    model.fit = model_SMM
    cost<-proc.time()-ptm
    coef.SMM  = model_SMM@results@fixed.effects

    # calculation of p-values for the 4 parameters
    tab = cbind(c(model_SMM@results@name.fixed, "residual standard error"),
                c(round(model_SMM@results@fixed.effects,3), round(model_SMM@results@respar[model_SMM@results@indx.res],3)),
                c(round(model_SMM@results@se.fixed,3),round(model_SMM@results@se.respar[model_SMM@results@indx.res],3)))
    colnames(tab) = c("Parameter","Estimate","  SE")
    wstat = as.double(tab[,2])/as.double(tab[,3])
    pval  = rep(0, length(wstat))
    pval2 = round(1 - normcdf(abs(wstat[1:length(pval)])), 4)
    tab   = cbind(tab,"p-value" = pval2)
    tab = as.data.frame(tab)
    tab$`p-value` = ifelse(tab$`p-value` == "0", "P<.0001", tab$`p-value`)
    tab
    # output
    print(tab)
    cat("----------------------------------------------------\n The program took", round(cost[3],2), "seconds \n")


    ### ESTIMATED MARGINAL TRAJECTORIES IN THE WHOLE SAMPLE (most common profile of covariates) ###
    if (traj.marg == TRUE){
      ##
      A     = coef.SMM[1]
      B     = coef.SMM[1*length(predictor.all)+length(predictor.par1) + 2]
      C     = coef.SMM[2*length(predictor.all)+length(predictor.par1)+length(predictor.par2) + 3]
      D     = coef.SMM[3*length(predictor.all)+length(predictor.par1)+length(predictor.par2)+length(predictor.par3) + 4]
      # marginal predictions
      if (mean(dataset[,time]) > 0){
        min_plot = round(quantile(dataset[,time], probs=c(0.01)),0)
        max_plot = round(quantile(dataset[,time], probs=c(0.90)),0)
        tab_SMM  = data.frame(time_pos = seq(min_plot, max_plot,0.1))
        tab_SMM$traj_SMM = SSlogis5(tab_SMM$time_pos, A, B, C, D, theta=1)

        plot(tab_SMM$time_pos, tab_SMM$traj_SMM,  las=1,
             lwd  = 3, type = "l",
             xlab = x.lab,  ylab = y.lab,
             main = main.traj.marg,
             col  = "cyan4",
             legend = NULL)

      } else if (mean(dataset[,time]) < 0) {
        min_plot = round(quantile(dataset[,time], probs=c(0.1)),0)
        max_plot = round(quantile(dataset[,time], probs=c(0.99)),0)
        tab_SMM  = data.frame(time_neg = seq(min_plot, max_plot,0.1), time_pos = seq(-min_plot, max_plot,-0.1))
        tab_SMM$traj_SMM = SSlogis5(tab_SMM$time_pos, A, B, C, D, theta=1)

        plot(tab_SMM$time_neg, tab_SMM$traj_SMM,  las=1,
             lwd  = 4, type = "l",
             xlab = x.lab,  ylab = y.lab,
             main = main.traj.marg,
             col  = "cyan4",
             legend = NULL)
      }
    }



    ######################################################
    ### ESTIMATED MARGINAL TRAJECTORIES BETWEEN GROUPS ###
    ######################################################
    ################################
    #   BINARY GROUPING VARIABLE   #
    ################################
    if (is.null(traj.marg.group)==F & is.factor(dataset[,traj.marg.group]) == TRUE){

      ## GROUP OF REFENCE (val == 0)
      A_group0 = coef.SMM[1]
      B_group0 = coef.SMM[1*length(predictor.all)+length(predictor.par1) + 2]
      C_group0 = coef.SMM[2*length(predictor.all)+length(predictor.par1)+length(predictor.par2) + 3]
      D_group0 = coef.SMM[3*length(predictor.all)+length(predictor.par1)+length(predictor.par2)+length(predictor.par3) + 4]
      ## SECOND GROUP (val == 1)
      pos_raw = which(traj.marg.group == rownames(mat_predictor))

      A_group1 = coef.SMM[1]+
        coef.SMM[1+pos_raw]*mat_predictor[pos_raw,1]
      B_group1 = coef.SMM[1*length(predictor.all)+length(predictor.par1) + 2]+
        coef.SMM[1*length(predictor.all)+length(predictor.par1) + 2 + pos_raw]*mat_predictor[pos_raw,2]
      C_group1 = coef.SMM[2*length(predictor.all)+length(predictor.par1)+length(predictor.par2) + 3]+
        coef.SMM[2*length(predictor.all)+length(predictor.par1)+length(predictor.par2) + 3 + pos_raw]*mat_predictor[pos_raw,3]

      if (is.na(coef.SMM[3*length(predictor.all)+length(predictor.par1)+length(predictor.par2)+length(predictor.par3) + 4 + pos_raw])==F){

        D_group1 = coef.SMM[3*length(predictor.all)+length(predictor.par1)+length(predictor.par2)+length(predictor.par3) + 4]+
          coef.SMM[3*length(predictor.all)+length(predictor.par1)+length(predictor.par2)+length(predictor.par3) + 4 + pos_raw]*mat_predictor[pos_raw,4]

      } else if (is.na(coef.SMM[3*length(predictor.all)+length(predictor.par1)+length(predictor.par2)+length(predictor.par3) + 4 + pos_raw])==T){

        D_group1 = coef.SMM[3*length(predictor.all)+length(predictor.par1)+length(predictor.par2)+length(predictor.par3) + 4]

      }

      # POSITIVE TIMESCALE #
      if (mean(dataset[,time]) > 0){
        min_x = round(quantile(dataset[,time], probs=c(0.01)),0)
        max_x = round(quantile(dataset[,time], probs=c(0.90)),0)
        tab_SMM  = data.frame(time_pos = seq(min_x, max_x,0.1))
        tab_SMM$traj_SMM_group0 = SSlogis5(tab_SMM$time_pos, A_group0, B_group0, C_group0, D_group0, theta=1)
        tab_SMM$traj_SMM_group1 = SSlogis5(tab_SMM$time_pos, A_group1, B_group1, C_group1, D_group1, theta=1)
        ##
        min_y = min(tab_SMM$traj_SMM_group0, tab_SMM$traj_SMM_group1)
        max_y = max(tab_SMM$traj_SMM_group0, tab_SMM$traj_SMM_group1)
        ##
        plot(tab_SMM$time_neg, tab_SMM$traj_SMM_group0,  las=1,
             lwd  = 5, type = "l",
             xlab = x.lab,  ylab = y.lab,
             ylim = c(min_y, max_y),
             main = main.traj.marg.group,
             col  = "cyan4",
             legend = NULL)
        points(tab_SMM$time_neg, tab_SMM$traj_SMM_group1, col="red",
               lty=3, lwd = 3, type = "l")
        # legend
        group0 = as.character(paste(traj.marg.group,"=",unique(dataset[,traj.marg.group])[1], collapse = ""))
        group1 = as.character(paste(traj.marg.group,"=",unique(dataset[,traj.marg.group])[2], collapse = ""))
        legend("bottomleft", bty = "n", lwd=3, lty = c(1,3),
               col=c("cyan4","red"), legend = c(group0,group1))

        # NEGATIVE TIMESCALE #
      } else if (mean(dataset[,time]) < 0) {

        min_x = round(quantile(dataset[,time], probs=c(0.1)),0)
        max_x = round(quantile(dataset[,time], probs=c(0.99)),0)
        tab_SMM  = data.frame(time_neg = seq(min_x, max_x,0.1), time_pos = seq(-min_x, max_x,-0.1))
        tab_SMM$traj_SMM_group0 = SSlogis5(tab_SMM$time_pos, A_group0, B_group0, C_group0, D_group0, theta=1)
        tab_SMM$traj_SMM_group1 = SSlogis5(tab_SMM$time_pos, A_group1, B_group1, C_group1, D_group1, theta=1)
        ##
        min_y = min(tab_SMM$traj_SMM_group0, tab_SMM$traj_SMM_group1)
        max_y = max(tab_SMM$traj_SMM_group0, tab_SMM$traj_SMM_group1)
        ##
        plot(tab_SMM$time_neg, tab_SMM$traj_SMM_group0,  las=1,
             lwd  = 5, type = "l",
             xlab = x.lab,  ylab = y.lab,
             ylim = c(min_y, max_y),
             main = main.traj.marg.group,
             col  = "cyan4",
             legend = NULL)
        points(tab_SMM$time_neg, tab_SMM$traj_SMM_group1, col="red",
               lty=3, lwd = 3, type = "l")
        # legend
        group0 = as.character(paste(traj.marg.group,"=",unique(dataset[,traj.marg.group])[1], collapse = ""))
        group1 = as.character(paste(traj.marg.group,"=",unique(dataset[,traj.marg.group])[2], collapse = ""))
        legend("bottomleft", bty = "n", lwd=3, lty = c(1,3),
               col=c("cyan4","red"), legend = c(group0,group1))
      }
    }



    ################################
    # CONTINUOUS GROUPING VARIABLE #
    ################################
    if (is.null(traj.marg.group)==F & is.numeric(dataset[,traj.marg.group]) == TRUE){

      bound = round(quantile(dataset[,traj.marg.group], probs = percentiles),1)
      pos_raw = which(traj.marg.group == rownames(mat_predictor))

      ## GROUP OF REFENCE (val in the 10th percentile)
      A_group0 = coef.SMM[1]+
        coef.SMM[1+pos_raw]*mat_predictor[pos_raw,1]*bound[1]
      B_group0 = coef.SMM[1*length(predictor.all)+length(predictor.par1) + 2]+
        coef.SMM[1*length(predictor.all)+length(predictor.par1) + 2 + pos_raw]*mat_predictor[pos_raw,2]*bound[1]
      C_group0 = coef.SMM[2*length(predictor.all)+length(predictor.par1)+length(predictor.par2) + 3]+
        coef.SMM[2*length(predictor.all)+length(predictor.par1)+length(predictor.par2) + 3 + pos_raw]*mat_predictor[pos_raw,3]*bound[1]

      if (is.na(coef.SMM[3*length(predictor.all)+length(predictor.par1)+length(predictor.par2)+length(predictor.par3) + 4 + pos_raw])==F){

        D_group0 = coef.SMM[3*length(predictor.all)+length(predictor.par1)+length(predictor.par2)+length(predictor.par3) + 4]+
          coef.SMM[3*length(predictor.all)+length(predictor.par1)+length(predictor.par2)+length(predictor.par3) + 4 + pos_raw]*mat_predictor[pos_raw,4]*bound[1]

      } else if (is.na(coef.SMM[3*length(predictor.all)+length(predictor.par1)+length(predictor.par2)+length(predictor.par3) + 4 + pos_raw])==T){

        D_group0 = coef.SMM[3*length(predictor.all)+length(predictor.par1)+length(predictor.par2)+length(predictor.par3) + 4]

      }

      ## SECOND GROUP (val in the 90th percentile)
      A_group1 = coef.SMM[1]+
        coef.SMM[1+pos_raw]*mat_predictor[pos_raw,1]*bound[2]
      B_group1 = coef.SMM[1*length(predictor.all)+length(predictor.par1) + 2]+
        coef.SMM[1*length(predictor.all)+length(predictor.par1) + 2 + pos_raw]*mat_predictor[pos_raw,2]*bound[2]
      C_group1 = coef.SMM[2*length(predictor.all)+length(predictor.par1)+length(predictor.par2) + 3]+
        coef.SMM[2*length(predictor.all)+length(predictor.par1)+length(predictor.par2) + 3 + pos_raw]*mat_predictor[pos_raw,3]*bound[2]

      if (is.na(coef.SMM[3*length(predictor.all)+length(predictor.par1)+length(predictor.par2)+length(predictor.par3) + 4 + pos_raw])==F){

        D_group1 = coef.SMM[3*length(predictor.all)+length(predictor.par1)+length(predictor.par2)+length(predictor.par3) + 4]+
          coef.SMM[3*length(predictor.all)+length(predictor.par1)+length(predictor.par2)+length(predictor.par3) + 4 + pos_raw]*mat_predictor[pos_raw,4]*bound[2]

      } else if (is.na(coef.SMM[3*length(predictor.all)+length(predictor.par1)+length(predictor.par2)+length(predictor.par3) + 4 + pos_raw])==T){

        D_group1 = coef.SMM[3*length(predictor.all)+length(predictor.par1)+length(predictor.par2)+length(predictor.par3) + 4]

      }

      # POSITIVE TIMESCALE #
      if (mean(dataset[,time]) > 0){
        min_x = round(quantile(dataset[,time], probs=c(0.01)),0)
        max_x = round(quantile(dataset[,time], probs=c(0.90)),0)
        tab_SMM  = data.frame(time_pos = seq(min_x, max_x,0.1))
        tab_SMM$traj_SMM_group0 = SSlogis5(tab_SMM$time_pos, A_group0, B_group0, C_group0, D_group0, theta=1)
        tab_SMM$traj_SMM_group1 = SSlogis5(tab_SMM$time_pos, A_group1, B_group1, C_group1, D_group1, theta=1)
        ##
        min_y = min(tab_SMM$traj_SMM_group0, tab_SMM$traj_SMM_group1)
        max_y = max(tab_SMM$traj_SMM_group0, tab_SMM$traj_SMM_group1)
        ##
        plot(tab_SMM$time_neg, tab_SMM$traj_SMM_group0,  las=1,
             lwd  = 5, type = "l",
             xlab = x.lab,  ylab = y.lab,
             ylim = c(min_y, max_y),
             main = main.traj.marg.group,
             col  = "cyan4",
             legend = NULL)
        points(tab_SMM$time_neg, tab_SMM$traj_SMM_group1, col="red",
               lty=3, lwd = 3, type = "l")
        # legend
        group0 = as.character(paste0(traj.marg.group," (",percentiles[1]*100,"th percentile)", collapse =""))
        group1 = as.character(paste0(traj.marg.group," (",percentiles[2]*100,"th percentile)", collapse =""))
        legend("bottomleft", bty = "n", lwd=3, lty = c(1,3),
               col=c("cyan4","red"), legend = c(group0,group1), cex=0.9)

        # NEGATIVE TIMESCALE #
      } else if (mean(dataset[,time]) < 0) {

        min_x = round(quantile(dataset[,time], probs=c(0.1)),0)
        max_x = round(quantile(dataset[,time], probs=c(0.99)),0)
        tab_SMM  = data.frame(time_neg = seq(min_x, max_x,0.1), time_pos = seq(-min_x, max_x,-0.1))
        tab_SMM$traj_SMM_group0 = SSlogis5(tab_SMM$time_pos, A_group0, B_group0, C_group0, D_group0, theta=1)
        tab_SMM$traj_SMM_group1 = SSlogis5(tab_SMM$time_pos, A_group1, B_group1, C_group1, D_group1, theta=1)
        ##
        min_y = min(tab_SMM$traj_SMM_group0, tab_SMM$traj_SMM_group1)
        max_y = max(tab_SMM$traj_SMM_group0, tab_SMM$traj_SMM_group1)
        ##
        plot(tab_SMM$time_neg, tab_SMM$traj_SMM_group0,  las=1,
             lwd  = 5, type = "l",
             xlab = x.lab,  ylab = y.lab,
             ylim = c(min_y, max_y),
             main = main.traj.marg.group,
             col  = "cyan4",
             legend = NULL)
        points(tab_SMM$time_neg, tab_SMM$traj_SMM_group1, col="red",
               lty=3, lwd = 3, type = "l")
        # legend
        group0 = as.character(paste0(traj.marg.group," (",percentiles[1]*100,"th percentile)", collapse =""))
        group1 = as.character(paste0(traj.marg.group," (",percentiles[2]*100,"th percentile)", collapse =""))
        legend("bottomleft", bty = "n", lwd=3, lty = c(1,3),
               col=c("cyan4","red"), legend = c(group0,group1), cex=0.9)
      }
    }
  }







  ######################################
  ######  PIECEWISE MIXED MODELS  ######
  ######    (models in c(2,3)     ######
  ######################################

  if (model == 2) {
    ## SPECIFICATION - ABRUPT CHANGE ##
    model_fct = function(psi, ID, xidep){
      t           = xidep[, 1]
      last.level  = psi[ID, 1]
      slope1      = psi[ID, 2]
      slope2      = psi[ID, 3]
      changepoint = psi[ID, 4]
      I1 = as.numeric(t <= changepoint)
      I2 = as.numeric(t >  changepoint)
      Y_pred = I1*(last.level + slope2*changepoint + slope1*(t - changepoint)) + I2*(last.level + slope2*t)
      return(Y_pred)
    }
  } else if (model == 3) {
    ## SPECIFICATION - POLYNOMIAL SMOOTH TRANSITION ##
    model_fct = function(psi, ID, xidep){
      t           = xidep[, 1]
      last.level  = psi[ID, 1]
      slope1      = psi[ID, 2]
      slope2      = psi[ID, 3]
      changepoint = psi[ID, 4]
      transition  = 2
      I1 = as.numeric(t <= changepoint + transition/2)
      I2 = as.numeric(t >  changepoint + transition/2)
      Y_pred = I1*(last.level + slope2*(changepoint + transition/2) + slope1*(t - (changepoint + transition/2))) +I2*(last.level + slope2*t)
      return(Y_pred)
    }
  }


  if (model == 2 | model == 3) {

    ## CORRELATION BETWEEN n1/n2 only ##
    varCov      = diag(1, ncol=4, nrow=4)
    varCov[2,3] = varCov[3,2] = 1

    ## STARTING VALUES OF THE 4 PARAMETERS ##
    if (is.null(start) == T){

      frag   = quantile(dataset$time, probs = c(0.2,0.4,0.6,0.8))
      # interval 1
      tab_int = subset(dataset, time <= frag[1])
      lmm    = hlme(outcome ~ 1 + time, random = ~1 + time, subject=ID, data=tab_int, verbose = F)
      slope1 = mean(coef(lmm)["time"])
      # interval 2
      tab_int = subset(dataset, time >= frag[1] & time <= frag[2])
      lmm  = hlme(outcome ~ 1 + time, random = ~1 + time, subject=ID, data=tab_int, verbose = F)
      slope2 = mean(coef(lmm)["time"])
      # interval 3
      tab_int = subset(dataset, time >= frag[2] & time <= frag[3])
      lmm  = hlme(outcome ~ 1 + time, random = ~1 + time, subject=ID, data=tab_int, verbose = F)
      slope3 = mean(coef(lmm)["time"])
      # interval 4
      tab_int = subset(dataset, time >= frag[3] & time <= frag[4])
      lmm  = hlme(outcome ~ 1 + time, random = ~1 + time, subject=ID, data=tab_int, verbose = F)
      slope4 = mean(coef(lmm)["time"])
      # interval 5
      tab_int = subset(dataset, time >= frag[4])
      lmm  = hlme(outcome ~ 1 + time, random = ~1 + time, subject=ID, data=tab_int, verbose = F)
      slope5 = mean(coef(lmm)["time"])
      ###
      y = c(slope1,slope2,slope3,slope4,slope5)
      val = min(y)
      pos = frag[which(y == val)-1]

      # last level (intercept)
      tempo      = subset(dataset, time > quantile(dataset$time, probs=c(0.97))) #KEEP
      last.level = mean(tempo[, outcome])
      # first slope
      tempo  = subset(dataset, time < pos)
      m1     = hlme(outcome ~ 1 + time, random = ~1 + time, subject=ID, data=tempo, verbose = F)
      slope1 = mean(coef(m1)["time"])
      # second slope
      tempo  = subset(dataset, time > pos)
      m2     = hlme(outcome ~ 1 + time, random = ~1+time, subject=ID, data=tempo, verbose = F)
      slope2 = mean(coef(m2)["time"])
      # changepoint
      changepoint = mean(pos)
      c(last.level, slope1, slope2, changepoint)

    } else if (length(start) == 4) {
      last.level  = start[1]
      slope1      = start[2]
      slope2      = start[3]
      changepoint = start[4]
    }
    c(last.level, slope1, slope2, changepoint)

    ## SPECIFICATION OF THE MATRIX OF COVARIATES ##

    if (nb_cov == 0) {

      mat_predictor = matrix(1, ncol=4, nrow=nb_cov)

    } else if (is.null(predictor.all)==F & nb_cov == length(predictor.all)){

      mat_predictor = matrix(1, ncol=4, nrow=nb_cov)
      y = unique(unlist(as.list(predictor.all)))
      y = as.list(y)
      rownames(mat_predictor) = y

    } else if (is.null(predictor.all)==F & nb_cov != length(predictor.all)){

      mat.all = matrix(1, ncol=4, nrow=length(predictor.all))
      y = unique(unlist(as.list(predictor.all)))
      y = as.list(y)
      rownames(mat.all) = y
      mat.tempo = matrix(0, ncol=4, nrow=nb_cov-length(predictor.all))
      x = unique(unlist(as.list(c(predictor.par1, predictor.par2, predictor.par3, predictor.par4))))
      x = as.list(x)
      rownames(mat.tempo) = x

      for (i in 1:length(x)){
        if (isTRUE(x[i] %in% predictor.par1)){mat.tempo[i,1] = 1}
      }
      for (i in 1:length(x)){
        if (isTRUE(x[i] %in% predictor.par2)){mat.tempo[i,2] = 1}
      }
      for (i in 1:length(x)){
        if (isTRUE(x[i] %in% predictor.par3)){mat.tempo[i,3] = 1}
      }
      for (i in 1:length(x)){
        if (isTRUE(x[i] %in% predictor.par4)){mat.tempo[i,4] = 1}
      }
      mat_predictor = rbind(mat.all,mat.tempo)
    }
    mat_predictor

    ##
    saemix.model.cog = saemixModel(model = model_fct,
                                   psi0  = rbind(c(last.level  = last.level,
                                                   slope1      = slope1,
                                                   slope2      = slope2,
                                                   changepoint = changepoint),
                                                 matrix(0, ncol=4, nrow=1)),
                                   covariate.model  = mat_predictor,
                                   covariance.model = varCov,
                                   verbose = F)

    ##
    saemix.cog = saemixData(name.data       = dataset,
                            name.group      = ID,
                            name.predictors = time,
                            name.response   = outcome,
                            name.covariates = predictors,
                            verbose = F)
    ##
    saemix.options = saemixControl(map = F, print = F,
                                   nbiter.saemix = c(300,100), # nb of iterations for the exploration and smoothing phase
                                   nbiter.burn = 20,           # nb of iterations for burning
                                   nbiter.mcmc = c(5,5,5,0),   # nb of iterations in each kernel during the MCMC step
                                   ll.is = T,                  # default = T (TRUE) to estimate the log-likelihood
                                   seed = 123,
                                   displayProgress = F,        # default = T = to output the convergence plots
                                   save = F,                   # default = T = the results of the fit should be saved to a file
                                   save.graphs = F)            # default = T = to save the diagnostic and individual graphs

    ptm<-proc.time()
    ## fit using the estimation function saemix()
    model_PMM = saemix(saemix.model.cog, saemix.cog, saemix.options)
    model.fit = model_PMM
    cost<-proc.time()-ptm
    coef.PMM  = model_PMM@results@fixed.effects

    # calculation of p-values for the 4 parameters
    tab = cbind(c(model_PMM@results@name.fixed, "error"),
                c(round(model_PMM@results@fixed.effects,4), round(model_PMM@results@respar[model_PMM@results@indx.res],4)),
                c(round(model_PMM@results@se.fixed,4),round(model_PMM@results@se.respar[model_PMM@results@indx.res],4)))
    colnames(tab) = c("Parameter","Estimate","  SE")
    wstat = as.double(tab[,2])/as.double(tab[,3])
    pval  = rep(0, length(wstat))
    pval2 = round(1 - normcdf(abs(wstat[1:length(pval)])), 4)
    tab   = cbind(tab,"p-value" = pval2)
    tab = as.data.frame(tab)
    tab$`p-value` = ifelse(tab$`p-value` == "0", "P<.0001", tab$`p-value`)
    tab


    ##########################################################
    ######   MARGINAL ESTIMATED TRAJECTORIES            ######
    ######   FOR THE MOST COMMON PROFILE OF COVARIATES  ######
    ##########################################################
    if(traj.marg == TRUE){
      if (mean(dataset[,time]) > 0){
        min_plot = round(quantile(dataset[,time], probs=c(0.01)),0)
        max_plot = round(quantile(dataset[,time], probs=c(0.90)),0)
      } else if (mean(dataset[,time]) < 0) {
        min_plot = round(quantile(dataset[,time], probs=c(0.1)),0)
        max_plot = round(quantile(dataset[,time], probs=c(0.99)),0)
      }

      tab_traj = data.frame(time = seq(min_plot, max_plot, 0.1),
                            ## parameters
                            last.level  = coef.PMM[1],
                            slope1      = coef.PMM[1*length(predictor.all)+length(predictor.par1) + 2],
                            slope2      = coef.PMM[2*length(predictor.all)+length(predictor.par1)+length(predictor.par2) + 3],
                            changepoint = coef.PMM[3*length(predictor.all)+length(predictor.par1)+length(predictor.par2)+length(predictor.par3) + 4],
                            ## transition parameter
                            transition = 2)

      n0 = tab_traj$last.level[1]
      n1 = tab_traj$slope1[1]
      n2 = tab_traj$slope2[1]
      n3 = tab_traj$changepoint[1]
      c(n0,n1,n2,n3)

      if (model == 2) {
        tab_traj$I1       = as.numeric(tab_traj$time <= tab_traj$changepoint)
        tab_traj$I2       = as.numeric(tab_traj$time >  tab_traj$changepoint)
        tab_traj$traj_PMM = tab_traj$I1*(n0 + n2*n3 + n1*(tab_traj$time - n3)) +
          tab_traj$I2*(n0 + n2*tab_traj$time)

      } else if (model == 3) {
        e  = tab_traj$transition[1]
        lambda = n0 + n2*(n3+e/2) - n1*(n3 + e/2)
        tab_traj$I1 = as.numeric(tab_traj$time <= tab_traj$changepoint)
        tab_traj$I2 = as.numeric(tab_traj$time >  tab_traj$changepoint & tab_traj$time < tab_traj$changepoint + 2)
        tab_traj$I3 = as.numeric(tab_traj$time >= tab_traj$changepoint + 2)

        # system of 4 linear equations with 4 unknown parameters
        A = matrix(c(n3^3, n3^2, n3, 1,
                     (n3+e)^3, (n3+e)^2, (n3+e), 1,
                     2*(n3^2), 2*n3, 1, 0,
                     2*((n3+e)^2), 2*(n3+e), 1, 0), ncol=4, byrow=T)
        B = matrix(c(n0+n2*(n3+e/2)-n1*(n3+e/2)+n1*n3,
                     n0+n2*(n3+e),
                     n1,
                     n2))
        mat = solve(A) %*% B

        tab_traj$traj_PMM =
          tab_traj$I1*(lambda+n1*tab_traj$time)+
          tab_traj$I2*(mat[1]*(tab_traj$time^3) + mat[2]*(tab_traj$time^2) + mat[3]*tab_traj$time + mat[4]) +
          tab_traj$I3*(n0+n2*tab_traj$time)
      }

      ## Plot
      plot(tab_traj$time, tab_traj$traj_PMM,
           lwd = 3, type = "l", las=1,
           xlim=c(min_plot, max_plot),
           xlab = x.lab,  ylab = y.lab,
           main = main.traj.marg,
           col  = "cyan4",
           legend = NULL)
    }


    ######################################################
    ### ESTIMATED MARGINAL TRAJECTORIES BETWEEN GROUPS ###
    ######################################################
    ################################
    #   BINARY GROUPING VARIABLE   #
    ################################
    if (is.null(traj.marg.group)==F & is.factor(dataset[,traj.marg.group]) == TRUE){

      if (mean(dataset[,time]) > 0){
        min_plot = round(quantile(dataset[,time], probs=c(0.01)),0)
        max_plot = round(quantile(dataset[,time], probs=c(0.90)),0)
      } else if (mean(dataset[,time]) < 0) {
        min_plot = round(quantile(dataset[,time], probs=c(0.1)),0)
        max_plot = round(quantile(dataset[,time], probs=c(0.99)),0)
      }

      ## GROUP OF REFENCE (val == 0)
      tab_traj = data.frame(time = seq(min_plot, max_plot, 0.1),
                            ## parameters
                            last.level  = coef.PMM[1],
                            slope1      = coef.PMM[1*length(predictor.all)+length(predictor.par1) + 2],
                            slope2      = coef.PMM[2*length(predictor.all)+length(predictor.par1)+length(predictor.par2) + 3],
                            changepoint = coef.PMM[3*length(predictor.all)+length(predictor.par1)+length(predictor.par2)+length(predictor.par3) + 4],
                            ## transition parameter
                            transition = 2)

      n0_group0 = tab_traj$last.level[1]
      n1_group0 = tab_traj$slope1[1]
      n2_group0 = tab_traj$slope2[1]
      n3_group0 = tab_traj$changepoint[1]

      ## GROUP OF REFENCE (val == 1)
      pos_raw = which(traj.marg.group == rownames(mat_predictor))
      ##
      tab_traj = data.frame(time = seq(min_plot, max_plot, 0.1), transition = 2)
      ##
      tab_traj$last.level  = coef.PMM[1] + coef.PMM[1+pos_raw]*mat_predictor[pos_raw,1]
      tab_traj$slope1      = coef.PMM[1*length(predictor.all)+length(predictor.par1) + 2] + coef.PMM[1*length(predictor.all)+length(predictor.par1) + 2 + pos_raw]*mat_predictor[pos_raw,2]
      tab_traj$slope2      = coef.PMM[2*length(predictor.all)+length(predictor.par1)+length(predictor.par2) + 3] + coef.PMM[2*length(predictor.all)+length(predictor.par1)+length(predictor.par2) + 3 + pos_raw]*mat_predictor[pos_raw,3]

      if (is.na(coef.PMM[3*length(predictor.all)+length(predictor.par1)+length(predictor.par2)+length(predictor.par3) + 4 + pos_raw])==F){

        tab_traj$changepoint = coef.PMM[3*length(predictor.all)+length(predictor.par1)+length(predictor.par2)+length(predictor.par3) + 4]+
          coef.PMM[3*length(predictor.all)+length(predictor.par1)+length(predictor.par2)+length(predictor.par3) + 4 + pos_raw]*mat_predictor[pos_raw,4]

      } else if (is.na(coef.PMM[3*length(predictor.all)+length(predictor.par1)+length(predictor.par2)+length(predictor.par3) + 4 + pos_raw])==T){

        tab_traj$changepoint = coef.PMM[3*length(predictor.all)+length(predictor.par1)+length(predictor.par2)+length(predictor.par3) + 4]

      }

      n0_group1 = tab_traj$last.level[1]
      n1_group1 = tab_traj$slope1[1]
      n2_group1 = tab_traj$slope2[1]
      n3_group1 = tab_traj$changepoint[1]


      if (model == 2) {
        ## GROUP OF REFENCE (val == 0) ##
        tab_traj$I1       = as.numeric(tab_traj$time <= tab_traj$changepoint)
        tab_traj$I2       = as.numeric(tab_traj$time >  tab_traj$changepoint)
        tab_traj$traj_PMM_group0 = tab_traj$I1*(n0_group0 + n2_group0*n3_group0 + n1_group0*(tab_traj$time - n3_group0)) +
          tab_traj$I2*(n0_group0 + n2_group0*tab_traj$time)
        ## SECOND GROUP (val == 1) ##
        tab_traj$traj_PMM_group1 = tab_traj$I1*(n0_group1 + n2_group1*n3_group1 + n1_group1*(tab_traj$time - n3_group1)) +
          tab_traj$I2*(n0_group1 + n2_group1*tab_traj$time)

      } else if (model == 3) {
        e  = tab_traj$transition[1]
        ## GROUP OF REFENCE (val == 0) ##
        lambda = n0_group0 + n2_group0*(n3_group0+e/2) - n1_group0*(n3_group0 + e/2)
        tab_traj$I1 = as.numeric(tab_traj$time <= tab_traj$changepoint)
        tab_traj$I2 = as.numeric(tab_traj$time >  tab_traj$changepoint & tab_traj$time < tab_traj$changepoint + 2)
        tab_traj$I3 = as.numeric(tab_traj$time >= tab_traj$changepoint + 2)

        # system of 4 linear equations with 4 unknown parameters
        A = matrix(c(n3_group0^3, n3_group0^2, n3_group0, 1,
                     (n3_group0+e)^3, (n3_group0+e)^2, (n3_group0+e), 1,
                     2*(n3_group0^2), 2*n3_group0, 1, 0,
                     2*((n3_group0+e)^2), 2*(n3_group0+e), 1, 0), ncol=4, byrow=T)
        B = matrix(c(n0_group0+n2_group0*(n3_group0+e/2)-n1_group0*(n3_group0+e/2)+n1_group0*n3_group0,
                     n0_group0+n2_group0*(n3_group0+e),
                     n1_group0,
                     n2_group0))
        mat = solve(A) %*% B

        tab_traj$traj_PMM_group0 =
          tab_traj$I1*(lambda+n1_group0*tab_traj$time)+
          tab_traj$I2*(mat[1]*(tab_traj$time^3) + mat[2]*(tab_traj$time^2) + mat[3]*tab_traj$time + mat[4]) +
          tab_traj$I3*(n0_group0+n2_group0*tab_traj$time)

        ## SECOND GROUP (val == 1) ##
        lambda = n0_group1 + n2_group1*(n3_group1+e/2) - n1_group1*(n3_group1 + e/2)
        tab_traj$I1 = as.numeric(tab_traj$time <= tab_traj$changepoint)
        tab_traj$I2 = as.numeric(tab_traj$time >  tab_traj$changepoint & tab_traj$time < tab_traj$changepoint + 2)
        tab_traj$I3 = as.numeric(tab_traj$time >= tab_traj$changepoint + 2)

        # system of 4 linear equations with 4 unknown parameters
        A = matrix(c(n3_group1^3, n3_group1^2, n3_group1, 1,
                     (n3_group1+e)^3, (n3_group1+e)^2, (n3_group1+e), 1,
                     2*(n3_group1^2), 2*n3_group1, 1, 0,
                     2*((n3_group1+e)^2), 2*(n3_group1+e), 1, 0), ncol=4, byrow=T)
        B = matrix(c(n0_group1+n2_group1*(n3_group1+e/2)-n1_group1*(n3_group1+e/2)+n1_group1*n3_group1,
                     n0_group1+n2_group1*(n3_group1+e),
                     n1_group1,
                     n2_group1))
        mat = solve(A) %*% B

        tab_traj$traj_PMM_group1 =
          tab_traj$I1*(lambda+n1_group1*tab_traj$time)+
          tab_traj$I2*(mat[1]*(tab_traj$time^3) + mat[2]*(tab_traj$time^2) + mat[3]*tab_traj$time + mat[4]) +
          tab_traj$I3*(n0_group1+n2_group1*tab_traj$time)
      }

      ##
      min_y = min(tab_traj$traj_PMM_group0, tab_traj$traj_PMM_group1)
      max_y = max(tab_traj$traj_PMM_group0, tab_traj$traj_PMM_group1)
      ##
      plot(tab_traj$time, tab_traj$traj_PMM_group0,  las=1,
           lwd  = 5, type = "l",
           xlab = x.lab,  ylab = y.lab,
           ylim = c(min_y, max_y),
           xlim=c(min_plot, max_plot),
           main = main.traj.marg.group,
           col  = "cyan4",
           legend = NULL)
      points(tab_traj$time, tab_traj$traj_PMM_group1, col="red",
             lty=3, lwd = 3, type = "l")
      # legend
      group0 = as.character(paste(traj.marg.group,"=",unique(dataset[,traj.marg.group])[1], collapse = ""))
      group1 = as.character(paste(traj.marg.group,"=",unique(dataset[,traj.marg.group])[2], collapse = ""))
      legend("bottomleft", bty = "n", lwd=3, lty = c(1,3),
             col=c("cyan4","red"), legend = c(group0,group1))

    }


    ####################################
    #   CONTINUOUS GROUPING VARIABLE   #
    ####################################
    if (is.null(traj.marg.group)==F & is.numeric(dataset[,traj.marg.group]) == TRUE){

      bound = round(quantile(dataset[,traj.marg.group], probs = percentiles),1)
      pos_raw = which(traj.marg.group == rownames(mat_predictor))

      if (mean(dataset[,time]) > 0){
        min_plot = round(quantile(dataset[,time], probs=c(0.01)),0)
        max_plot = round(quantile(dataset[,time], probs=c(0.90)),0)
      } else if (mean(dataset[,time]) < 0) {
        min_plot = round(quantile(dataset[,time], probs=c(0.1)),0)
        max_plot = round(quantile(dataset[,time], probs=c(0.99)),0)
      }
      c(min_plot,max_plot)

      ## GROUP OF REFERENCE (val in the 10th percentile)
      tab_traj_group0 = data.frame(time = seq(min_plot, max_plot, 0.1), transition = 2)
      ##
      tab_traj_group0$last.level  = coef.PMM[1] + coef.PMM[1+pos_raw]*mat_predictor[pos_raw,1]*bound[1]
      tab_traj_group0$slope1      = coef.PMM[1*length(predictor.all)+length(predictor.par1) + 2] + coef.PMM[1*length(predictor.all)+length(predictor.par1) + 2 + pos_raw]*mat_predictor[pos_raw,2]*bound[1]
      tab_traj_group0$slope2      = coef.PMM[2*length(predictor.all)+length(predictor.par1)+length(predictor.par2) + 3] + coef.PMM[2*length(predictor.all)+length(predictor.par1)+length(predictor.par2) + 3 + pos_raw]*mat_predictor[pos_raw,3]*bound[1]

      if (is.na(coef.PMM[3*length(predictor.all)+length(predictor.par1)+length(predictor.par2)+length(predictor.par3) + 4 + pos_raw])==F){

        tab_traj_group0$changepoint = coef.PMM[3*length(predictor.all)+length(predictor.par1)+length(predictor.par2)+length(predictor.par3) + 4]+
          coef.PMM[3*length(predictor.all)+length(predictor.par1)+length(predictor.par2)+length(predictor.par3) + 4 + pos_raw]*mat_predictor[pos_raw,4]*bound[1]

      } else if (is.na(coef.PMM[3*length(predictor.all)+length(predictor.par1)+length(predictor.par2)+length(predictor.par3) + 4 + pos_raw])==T){

        tab_traj_group0$changepoint = coef.PMM[3*length(predictor.all)+length(predictor.par1)+length(predictor.par2)+length(predictor.par3) + 4]

      }

      n0_group0 = tab_traj_group0$last.level[1]
      n1_group0 = tab_traj_group0$slope1[1]
      n2_group0 = tab_traj_group0$slope2[1]
      n3_group0 = tab_traj_group0$changepoint[1]

      if (model == 2) {
        ## GROUP OF REFENCE (val == 0) ##
        tab_traj_group0$I1       = as.numeric(tab_traj_group0$time <= tab_traj_group0$changepoint)
        tab_traj_group0$I2       = as.numeric(tab_traj_group0$time >  tab_traj_group0$changepoint)
        tab_traj_group0$traj_PMM_group0 = tab_traj_group0$I1*(n0_group0 + n2_group0*n3_group0 + n1_group0*(tab_traj_group0$time - n3_group0)) +
          tab_traj_group0$I2*(n0_group0 + n2_group0*tab_traj_group0$time)
      } else if (model == 3) {
        e  = tab_traj_group0$transition[1]
        ## GROUP OF REFENCE (val == 0) ##
        lambda = n0_group0 + n2_group0*(n3_group0+e/2) - n1_group0*(n3_group0 + e/2)
        tab_traj_group0$I1 = as.numeric(tab_traj_group0$time <= tab_traj_group0$changepoint)
        tab_traj_group0$I2 = as.numeric(tab_traj_group0$time >  tab_traj_group0$changepoint & tab_traj_group0$time < tab_traj_group0$changepoint + 2)
        tab_traj_group0$I3 = as.numeric(tab_traj_group0$time >= tab_traj_group0$changepoint + 2)

        # system of 4 linear equations with 4 unknown parameters
        A = matrix(c(n3_group0^3, n3_group0^2, n3_group0, 1,
                     (n3_group0+e)^3, (n3_group0+e)^2, (n3_group0+e), 1,
                     2*(n3_group0^2), 2*n3_group0, 1, 0,
                     2*((n3_group0+e)^2), 2*(n3_group0+e), 1, 0), ncol=4, byrow=T)
        B = matrix(c(n0_group0+n2_group0*(n3_group0+e/2)-n1_group0*(n3_group0+e/2)+n1_group0*n3_group0,
                     n0_group0+n2_group0*(n3_group0+e),
                     n1_group0,
                     n2_group0))
        mat = solve(A) %*% B

        tab_traj_group0$traj_PMM_group0 =
          tab_traj_group0$I1*(lambda+n1_group0*tab_traj_group0$time)+
          tab_traj_group0$I2*(mat[1]*(tab_traj_group0$time^3) + mat[2]*(tab_traj_group0$time^2) + mat[3]*tab_traj_group0$time + mat[4]) +
          tab_traj_group0$I3*(n0_group0+n2_group0*tab_traj_group0$time)
      }


      ## SECOND GROUP (val in the 90th percentile)
      tab_traj_group1 = data.frame(time = seq(min_plot, max_plot, 0.1), transition = 2)
      ##
      tab_traj_group1$last.level  = coef.PMM[1] + coef.PMM[1+pos_raw]*mat_predictor[pos_raw,1]*bound[2]
      tab_traj_group1$slope1      = coef.PMM[1*length(predictor.all)+length(predictor.par1) + 2] + coef.PMM[1*length(predictor.all)+length(predictor.par1) + 2 + pos_raw]*mat_predictor[pos_raw,2]*bound[2]
      tab_traj_group1$slope2      = coef.PMM[2*length(predictor.all)+length(predictor.par1)+length(predictor.par2) + 3] + coef.PMM[2*length(predictor.all)+length(predictor.par1)+length(predictor.par2) + 3 + pos_raw]*mat_predictor[pos_raw,3]*bound[2]

      if (is.na(coef.PMM[3*length(predictor.all)+length(predictor.par1)+length(predictor.par2)+length(predictor.par3) + 4 + pos_raw])==F){

        tab_traj_group1$changepoint = coef.PMM[3*length(predictor.all)+length(predictor.par1)+length(predictor.par2)+length(predictor.par3) + 4]+
          coef.PMM[3*length(predictor.all)+length(predictor.par1)+length(predictor.par2)+length(predictor.par3) + 4 + pos_raw]*mat_predictor[pos_raw,4]*bound[2]

      } else if (is.na(coef.PMM[3*length(predictor.all)+length(predictor.par1)+length(predictor.par2)+length(predictor.par3) + 4 + pos_raw])==T){

        tab_traj_group1$changepoint = coef.PMM[3*length(predictor.all)+length(predictor.par1)+length(predictor.par2)+length(predictor.par3) + 4]

      }

      n0_group1 = tab_traj_group1$last.level[1]
      n1_group1 = tab_traj_group1$slope1[1]
      n2_group1 = tab_traj_group1$slope2[1]
      n3_group1 = tab_traj_group1$changepoint[1]


      if (model == 2) {
        ## SECOND GROUP (val == 1) ##
        tab_traj_group1$I1       = as.numeric(tab_traj_group1$time <= tab_traj_group1$changepoint)
        tab_traj_group1$I2       = as.numeric(tab_traj_group1$time >  tab_traj_group1$changepoint)
        tab_traj_group1$traj_PMM_group1 = tab_traj_group1$I1*(n0_group1 + n2_group1*n3_group1 + n1_group1*(tab_traj_group1$time - n3_group1)) +
          tab_traj_group1$I2*(n0_group1 + n2_group1*tab_traj_group1$time)
      } else if (model == 3) {
        e  = tab_traj_group1$transition[1]
        ## SECOND GROUP (val == 1) ##
        lambda = n0_group1 + n2_group1*(n3_group1+e/2) - n1_group1*(n3_group1 + e/2)
        tab_traj_group1$I1 = as.numeric(tab_traj_group1$time <= tab_traj_group1$changepoint)
        tab_traj_group1$I2 = as.numeric(tab_traj_group1$time >  tab_traj_group1$changepoint & tab_traj_group1$time < tab_traj_group1$changepoint + 2)
        tab_traj_group1$I3 = as.numeric(tab_traj_group1$time >= tab_traj_group1$changepoint + 2)

        # system of 4 linear equations with 4 unknown parameters
        A = matrix(c(n3_group1^3, n3_group1^2, n3_group1, 1,
                     (n3_group1+e)^3, (n3_group1+e)^2, (n3_group1+e), 1,
                     2*(n3_group1^2), 2*n3_group1, 1, 0,
                     2*((n3_group1+e)^2), 2*(n3_group1+e), 1, 0), ncol=4, byrow=T)
        B = matrix(c(n0_group1+n2_group1*(n3_group1+e/2)-n1_group1*(n3_group1+e/2)+n1_group1*n3_group1,
                     n0_group1+n2_group1*(n3_group1+e),
                     n1_group1,
                     n2_group1))
        mat = solve(A) %*% B

        tab_traj_group1$traj_PMM_group1 =
          tab_traj_group1$I1*(lambda+n1_group1*tab_traj_group1$time)+
          tab_traj_group1$I2*(mat[1]*(tab_traj_group1$time^3) + mat[2]*(tab_traj_group1$time^2) + mat[3]*tab_traj_group1$time + mat[4]) +
          tab_traj_group1$I3*(n0_group1+n2_group1*tab_traj_group1$time)
      }

      ##
      min_y = min(tab_traj_group0$traj_PMM_group0, tab_traj_group1$traj_PMM_group1)
      max_y = max(tab_traj_group0$traj_PMM_group0, tab_traj_group1$traj_PMM_group1)
      ##
      plot(tab_traj_group0$time, tab_traj_group0$traj_PMM_group0,  las=1,
           lwd  = 5, type = "l",
           xlab = x.lab,  ylab = y.lab,
           ylim = c(min_y, max_y),
           xlim=c(min_plot, max_plot),
           main = main.traj.marg.group,
           col  = "cyan4",
           legend = NULL)
      points(tab_traj_group1$time, tab_traj_group1$traj_PMM_group1, col="red",
             lty=3, lwd = 3, type = "l")
      # legend
      group0 = as.character(paste0(traj.marg.group," (",percentiles[1]*100,"th percentile)", collapse =""))
      group1 = as.character(paste0(traj.marg.group," (",percentiles[2]*100,"th percentile)", collapse =""))
      legend("bottomleft", bty = "n", lwd=3, lty = c(1,3),
             col=c("cyan4","red"), legend = c(group0,group1), cex=0.9)
    }
  }
  # output
  print(list(model.fit, tab))
  cat("----------------------------------------------------\n The program took", round(cost[3],2), "seconds \n")
  return(model.fit)

}

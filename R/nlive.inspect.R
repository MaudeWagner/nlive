
#' Generation of key plots for a longitudinal variable of interest
#'
#'
#' The nlive.inspect() function allows to generate basic graphs to describe the longitudinal observed
#' measures of a variable of interest in the dataset
#'
#' @param dataset data frame containing the ID, variable, and time.
#' @param ID name of the variable representing the grouping structure specified with " (e.g., "ID" representing the unique identifier of participants).
#' @param variable name of the time-varying variable of interest specified with " (e.g., "variable").
#' @param time name of the variable representing the timescale specified with " (e.g., "time"). Can be negative or positive.
#' @param plot.xlabel optional text for the title of the x-axis of all plots.
#' @param plot.ylabel optional text for the title of the y-axis of all plots.
#' @param spag.plot.title optional text for the title of the spaghetti plot.
#'
#'
#' @return The _nlive.inspect_ function automatically provides
#' (i) an histogram of all the repeated measures of the variable available over time,
#' (ii) a spaghetti plot of the longitudinal observed variable for 70 randomly selected statistical units,
#' (iii) repeated boxplots of the longitudinal observed variable for each time unit.
#' The outputs are printed on the terminal and the numerical and graphical outputs are stored in a directory
#'
#'
#' @author Maude Wagner, Ana W. Capuano, Emmanuelle Comets
#'
#' \email{maude_wagner@@rush.edu}
#'
#' @references
#'
#' Capuano AW, Wagner M. nlive: an R package to facilitate the application of the sigmoidal and random changepoint mixed models. BMC Medical Research Methodology. 2023;23(1):257.
#' Hadley Wickham (2016). ggplot2: elegant graphics for data analysis. Springer.
#'
#' @examples
#'
#' \dontrun{
#' nlive.inspect(dataset=dataCog, ID="ID", variable="cognition", time="time")
#' }
#'
#' @import Rmisc
#' @import sitar
#' @import knitr
#' @importFrom dplyr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr row_number
#' @importFrom dplyr sample_n
#' @importFrom dplyr filter
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_histogram
#' @importFrom ggplot2 ggsave
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 geom_boxplot
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 guides
#' @importFrom graphics points
#' @importFrom graphics legend
#' @importFrom graphics plot
#' @importFrom stats quantile
#' @importFrom stats na.omit
#' @importFrom sqldf sqldf
#'
#' @export
########

nlive.inspect <- function(dataset, ID, variable, time,
         plot.xlabel = NULL,
         plot.ylabel = NULL,
         spag.plot.title = NULL){

  if (missing(dataset))
    stop("The argument dataset should be specified and defined as a data.frame")
  if (nrow(dataset) == 0)
    stop("Data should not be empty")
  if (missing(variable))
    stop("The argument variable must be specified")
  if (!is.numeric(dataset[[variable]]))
    stop("The argument variable must be numeric")
  if (missing(time))
    stop("The argument time must be specified")
  if (!is.numeric(dataset[[time]]))
    stop("The argument time must be numeric")


  ## dataset
  dataset$ID      = na.omit(dataset[,ID])
  dataset$variable = na.omit(dataset[,variable])
  dataset$time    = na.omit(dataset[,time])


  ## Values used throughout the code
  first    = as.data.frame(dataset %>% group_by(ID) %>% filter(row_number(ID) == 1))
  nb_indiv = dim(first)[1]

  ## Options
  # common to all plots #
  if (is.null(plot.xlabel) == T){x.lab = time} else if (is.null(plot.xlabel) == F){x.lab = plot.xlabel}
  if (is.null(plot.ylabel) == T){y.lab = variable} else if (is.null(plot.ylabel) == F){y.lab = plot.ylabel}
  # specific #
  if (is.null(spag.plot.title)  == T){main.spag.lab ="Spaghetti plot, random sample of 70 individuals"} else if (is.null(spag.plot.title) == F){main.spag.lab = spag.plot.title}


  #############################################
  ######           HISTOGRAM             ######
  #############################################
  histo  =   ggplot(dataset, aes(x = variable)) +
    geom_histogram(bins = 15, lwd=0.2, colour="white", alpha = 0.7) +
    labs(x = "", y = y.lab, title = "Histogram, whole sample")
  ggsave("inspect.histo.pdf", scale = 0.75)
  #
  print(histo)



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
    plot   = ggplot(tempo2, aes(time, variable, group = ID, colour = factor(ID))) +
      geom_line(size=0.8) +
      scale_x_continuous(limits=c(min, max)) +
      guides(colour="none") +
      labs(x=x.lab, y=y.lab, title = "Spaghetti plot, random selection of 70 individuals")
    ggsave("inspect.spagh.pdf", scale = 0.75)
    print(plot)

  } else if (nb_indiv < 70){

    ## Random selection for the whole study sample
    tempo1 = sample_n(first, dim(first)[1])
    tempo2 = sqldf('SELECT * FROM dataset WHERE ID IN (SELECT ID FROM tempo1)')
    min    = min(tempo2[,time])
    max    = max(tempo2[,time])
    plot   = ggplot(tempo2, aes(time, variable, group = ID, colour = factor(ID))) +
      geom_line(size=0.8) +
      scale_x_continuous(limits=c(min, max)) +
      guides(colour="none") +
      labs(x=x.lab, y=y.lab, title = "Spaghetti plot, whole sample")
    ggsave(path = "/myDirectory",
           device = "pdf", filename = "inspect.spagh", plot = plot)
    #
    print(plot)

  }


  #############################################
  ######            BOXPLOTS             ######
  #############################################
  box_plot  = ggplot(dataset, aes(x =round(time,0), y = variable)) +
    geom_boxplot(aes(fill=as.factor(round(time,0)))) +
    scale_y_continuous(limits=c(min(dataset$variable), max(dataset$variable))) +
    labs(x = x.lab, y = y.lab, title = "Boxplots over time, whole sample")+
    theme(legend.position = "none")
  ggsave("inspect.box.pdf", scale = 0.75)
  #
  print(box_plot)

}



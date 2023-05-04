# packages ----

library(tidyverse)
library(INLA)
library(latex2exp)

# directories ---- 

# extract file location of this script
codePath <- rstudioapi::getActiveDocumentContext()$path
codePathSplitted <- strsplit(codePath, "/")[[1]]

# retrieve directories
homeDir <- paste(codePathSplitted[1: (length(codePathSplitted)-3)], collapse = "/")
codeDir <- paste0(homeDir, '/Code')
dataDir <- paste0(homeDir, '/Data')

# need to make sure results is made

if(!dir.exists(paths = paste0(codeDir, '/Results'))) {
  dir.create(path = paste0(codeDir, '/Results'))
}

resultsDir <- paste0(codeDir, '/Results')

# saving format ----

height <- width <- 10
textSize <- 17.5

# load ----

## functions ----

source(paste0(codeDir, '/functions.R'))

# pc priors ----

n <- 100000
allPrec <- seq(from = 0.05, to = 1000, length.out = n) # precisions
allSD <- 1/sqrt(allPrec)

results <- 
  data.frame(allPrec = rep(allPrec, times = 2),
             allSD = rep(allSD, times = 2),
             prec = c(INLA::inla.pc.dprec(prec = allPrec, u = 1, alpha = 0.01), 
                      INLA::inla.pc.dprec(prec = allPrec, u = 3, alpha = 0.01)),
             sd = c(dexp(allSD, rate = -log(0.01)/1), 
                    dexp(allSD, rate = -log(0.01)/3)),
             U = c(rep(1, times = n),
                       rep(3, times = n)) %>% 
               factor(., levels = c(1, 3), labels = c(1, 3))) %>% 
  tibble::as_tibble()

pcPriorsPlot <- 
  ggplot2::ggplot(data = results, aes(x = allSD, y = sd, group = U, colour = U)) +
  ggplot2::geom_line() +
  ggplot2::labs(x = 'Standard Deviation', y = '') +
  ggplot2::scale_color_manual(labels = c(latex2exp::TeX(r'($\alpha = 0.01 \ and \ U = 1$)'), 
                                           latex2exp::TeX(r'($\alpha = 0.01 \ and \ U = 3$)')),
                                values = c('red3', 'blue3')) +
  my.theme(legend.title = element_blank(),
           legend.position = 'top',
           text = element_text(size = textSize)); pcPriorsPlot

ggplot2::ggsave(filename = paste0(resultsDir, '/pcPriorsPlot.png'),
                plot = pcPriorsPlot,
                height = height, width = width)

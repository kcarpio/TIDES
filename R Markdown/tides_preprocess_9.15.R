# Tides Snyder qPCR Analysis Combined preprocess (9/15/20)
tides <- read.table(pipe("pbpaste"), sep="\t", header = TRUE) 

library(PerformanceAnalytics)
library(viridis)
library(dplyr)
library(ggplot2)

# matrix showing bivariate scatterplots, distributions, and correlation values with significance levels
tt <- tides[which(tides$Tissuetype == 'BP'
                    & tides$Replicatenum == 1), c(4:10)]
chart.Correlation(tt, histogram = TRUE, pch = 19)
mtext("Tides BP rep.1", side = 3, line = 3)

# convert columns to factors 
tides$Replicatenum <- as.factor(tides$Replicatenum)
tides$`Sample Name` <- as.factor(tides$`Sample Name`)

# simple dot correlations
num <- tides[which(tides$Tissuetype == 'BP'
              & tides$Replicatenum == 1), c(1:2, 4:10)]
numvars <- which(sapply(num, is.numeric))
cat('There are', length(numvars), 'numeric variables')
all_num <- num[, numvars]
cor_numvars <- cor(all_num, use="pairwise.complete.obs")
corrplot.mixed(cor_numvars, tl.col="black", tl.pos = "lt", tl.cex = 1, cl.cex = 1, number.cex=1, 
               title="Tides BP rep.1",mar=c(0,0,1,0),tl.offset = 0.5)

# sort on decreasing correlations 
# cor_sorted <- as.matrix(sort(cor_numvar[,'SDHA'], decreasing = TRUE))
# select only high corelations
# corhigh <- names(which(apply(cor_sorted, 1, function(x) abs(x)>0.01)))
# cor_numvars <- cor_numvars[corhigh, corhigh]

# RNA concentration by tissue type
ggplot(tidesgm, aes(x = RNA_concentration, y = Tissuetype, fill = ..x..)) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  scale_fill_viridis(name = "RNA_concentration", option = "C") +
  labs(title = 'RNA concentration by tissue type and replicate') +
  facet_grid(. ~ Replicatenum)

# shapiro-wilk normality test (pre-norm)
myfun <- function(x, group) {
  data.frame(x, group) %>%
    group_by(group) %>%
    summarise(
      statistic = ifelse(sd(x)!=0,shapiro.test(x)$statistic,NA), 
      p.value = ifelse(sd(x)!=0,shapiro.test(x)$p.value,NA)
    )
}
(lst <- lapply(tides[,-1], myfun, group=tides[,1]))
# CS rep.1 18S: 0.0001031; CS rep.2 18S: 0.0006744; CS rep.1 SNORA46: 0.00562

# geometric means for 18S, SDHA, UBP1
tides$hkg.geom <- apply(tides[c(4:6)], 1, 
                        function(x) (prod(x[x!=0]))^(1/sum(x!=0)))
geomeans <- tides$hkg.geom

# normalize by geometric means
size.factors <- geomeans(data.matrix(tides))
tides.norm <- t(apply(tides, 1, function(x){ x / size.factors }))

# compare distributions pre- and post- normalization 
library("getopt")
print_use <- function(file=stderr())
{
  cat("norm-check: 
use: norm-check [options] < expr-matrix.txt
  options:
    -a|--after: FILE     plot distributions of expression values post-
                         norm to this PNG file: 'post-norm.png'
    -b|--before: FILE    plot distributions of expression values pre-
                         norm to this PNG file: 'pre-norm.png'
    -r|--res: INT        resolution (dpi), default is 150
    -t|--height: INT     height (pixels), default is 1200
    -w|--width: INT      width (pixels), default is 1200
    -y|--ylim: REAL      the visible range of the Y axis depends on the first
                         sample distribution plotted; if other distributions are
                         getting cut off, can use this setting to override the
                         default\n\n")
}

spec <- matrix( c("after",  'a', 1, "character",
                  "before", 'b', 1, "character",
                  "res",    'r', 1, "integer",
                  "height", 't', 1, "integer",
                  "width",  'w', 1, "integer",
                  "ylim",   'y', 1, "double"),
                byrow=TRUE, ncol=4)

opt  <- getopt(spec)
if(is.null(opt$after))  { opt$after  <- "post-norm.png" }
if(is.null(opt$before)) { opt$before <- "pre-norm.png"  }
if(is.null(opt$height)) { opt$height <- 1200            }
if(is.null(opt$res))    { opt$res    <- 150             }
if(is.null(opt$width))  { opt$width  <- 1200            }
if(!is.null(opt$ylim))  { opt$ylim   <- c(0, opt$ylim)  }

# determine number of samples
nsamp <- dim(tidesgm)[2] - 1
tides.norm  <- tidesgm[,1:nsamp+1]

# distribution of RNA expression values post-normalization
png(opt$after, height=opt$height, width=opt$width, res=opt$res)
h <- hist(log(tides.norm[,1]), plot=FALSE)
plot(h$mids, h$density, type="l", col=rainbow(nsamp)[1], main="",
     xlab="RNA expression value", ylab="Tissue type",
     ylim=opt$ylim)
for(i in 2:nsamp)
{
  h <- hist(log(tides.norm[,i]), plot=FALSE)
  lines(h$mids, h$density, col=rainbow(nsamp)[i])
}
devnum <- dev.off()

if(!is.null(opt$out))
{
  write.table(data.norm, file=opt$out, sep="\t")
}


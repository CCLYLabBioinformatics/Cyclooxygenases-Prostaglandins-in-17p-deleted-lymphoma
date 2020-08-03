##################################################
#Cor of PTGS genes with 17pGenes in TCGA-LAML ####
#Cor of PTGS genes with 17pGenes in TCGA-LAML ####
##################################################
sel <- read.csv("./TCGA_AML_all_sample_FPKM.csv")
rownames(sel) <- sel$X
sel <- sel[,-1]

P17_gene_all <- c("SCO1","MYH3","MYH2","MYH1","MYH4","MYH8","MYH13","RPS27AP1","GAS7","RCVRN","GLP2R","DHRS7C","USP43","WDR16","STX8","NTN1","PIK3R5","PIK3R6","MFSD6L","SPDYE4","CCDC42","RPS26P53","MYH10","NDEL1","RNF222","RPL26","KRBA2","ODF4","ARHGEF15","RANGRF","PFAS","TRNAI24","TRNAS14","TRNAT9","TRNAP8","TRNAD17","TRNAG18","TRNAW4","AURKB","TRNAI23","TRNAT8","TRNAS23","TRNAW3","SNORD118","TMEM107","VAMP2","PER1","TRNAT22","TRNAS6","TRNAG2","TRNAR1","HES7","TRNAL2","TRNAQ1","TRNAK1","ALOXE3","ALOX12B","ALOX15B","GUCY2D","CNTROB","TRAPPC1","KCNAB3","SCARNA21","CHD3","CYB5D1","LSMD1","TMEM88","KDM6B","RPL29P2","DNAH2","EFNB3","WRAP53","TP53","ATP1B2","SAT2","SHBG","FXR2","SOX15","MPDU1","CD68","SNORA67","SNORD10","SNORA48","EIF4A1","SENP3","TNFSF13","TNFSF12-TNFSF13","TNFSF12","POLR2A","AMAC1L3","ZBTB4","CHRNB1","FGF11","TMEM102","C17orf74","SPEM1","NLGN2","PLSCR3","TNK1","TMEM95","KCTD11","ACAP1","NEURL4","GPS2","EIF5A","YBX2","SLC2A4","CLDN7","DULLARD","GABARAP","PHF23","DVL2","MIR324","ACADVL","DLG4","ASGR1","RPL7AP64","ASGR2","CLEC10A","SLC16A11","SLC16A13","BCL6B","ALOX12")
P17_all <- c()
for (i in P17_gene_all){
	tmp <- subset(sel,symbol==i)
	P17_all <- rbind(P17_all,tmp)
	print(paste(i, "is done",sep=" "))
}
rownames(P17_all) <- P17_all$symbol
P17_all <- P17_all[,-ncol(P17_all)]
a <- as.data.frame(colSums(P17_all)/nrow(P17_all))
a <- t(a)
rownames(a) <- "P17_del_all"

PTGS1 <- subset(sel,symbol=="PTGS1")
PTGS1 <- PTGS1[,-ncol(PTGS1)]
rownames(PTGS1) <- "PTGS1"
PTGS2 <- subset(sel,symbol=="PTGS2")
PTGS2 <- PTGS2[,-ncol(PTGS2)]
rownames(PTGS2) <- "PTGS2"
c <- rbind(a,PTGS1)
c <- rbind(c,PTGS2)
c <- as.data.frame(t(c))
c$patient_id <- as.factor(rownames(c))
dat <- c

library(GGally)
my_custom_cor <- function(data, mapping, color = I("grey50"), sizeRange = c(1, 5), ...) {

  # get the x and y data to use the other code
  x <- GGally::eval_data_col(data, mapping$x)
  y <- GGally::eval_data_col(data, mapping$y)

  ct <- cor.test(x,y)
  sig <- symnum(
    ct$p.value, corr = FALSE, na = FALSE,
    cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
    symbols = c("***", "**", "*", ".", " ")
  )

  r <- unname(ct$estimate)
  rt <- format(r, digits=2)[1]

  # since we can't print it to get the strsize, just use the max size range
  cex <- max(sizeRange)

  # helper function to calculate a useable size
  percent_of_range <- function(percent, range) {
    percent * diff(range) + min(range, na.rm = TRUE)
  }

  # plot the cor value
  ggally_text(
    label = as.character(rt), 
    mapping = aes(),
    xP = 0.5, yP = 0.5, 
    size = I(percent_of_range(cex * abs(r), sizeRange)),
    color = color,
    ...
  ) + 
    # add the sig stars
    geom_text(
      aes_string(
        x = 0.8,
        y = 0.8
      ),
      label = sig, 
      size = I(cex),
      color = color,
      ...
    ) + 
    # remove all the background stuff and wrap it with a dashed line
    theme_classic() + 
    theme(
      panel.background = element_rect(
        color = color, 
        linetype = "longdash"
      ), 
      axis.line = element_blank(), 
      axis.ticks = element_blank(), 
      axis.text.y = element_blank(), 
      axis.text.x = element_blank()
    )
}
my_custom_cor(dat, aes(P17_del_all,PTGS1,PTGS2))
my_custom_smooth <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_point(color = I("blue")) + 
    geom_smooth(method = "lm", color = I("black"), ...)
}
my_custom_smooth(iris, aes(Sepal.Length, Sepal.Width))

pm <- ggpairs(
  dat[,1:3], 
  upper = list(continuous = wrap(my_custom_cor,sizeRange = c(6,10))),
  lower = list(continuous = my_custom_smooth), 
  axisLabels = "show"
)
pm



dat[,c(1:3)] <- log(dat[,c(1:3)]+1,2)
my_custom_cor <- function(data, mapping, color = I("grey50"), sizeRange = c(1, 5), ...) {

  # get the x and y data to use the other code
  x <- GGally::eval_data_col(data, mapping$x)
  y <- GGally::eval_data_col(data, mapping$y)

  ct <- cor.test(x,y)
  sig <- symnum(
    ct$p.value, corr = FALSE, na = FALSE,
    cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
    symbols = c("***", "**", "*", ".", " ")
  )

  r <- unname(ct$estimate)
  rt <- format(r, digits=2)[1]

  # since we can't print it to get the strsize, just use the max size range
  cex <- max(sizeRange)

  # helper function to calculate a useable size
  percent_of_range <- function(percent, range) {
    percent * diff(range) + min(range, na.rm = TRUE)
  }

  # plot the cor value
  ggally_text(
    label = as.character(rt), 
    mapping = aes(),
    xP = 0.5, yP = 0.5, 
    size = I(percent_of_range(cex * abs(r), sizeRange)),
    color = color,
    ...
  ) + 
    # add the sig stars
    geom_text(
      aes_string(
        x = 0.8,
        y = 0.8
      ),
      label = sig, 
      size = I(cex),
      color = color,
      ...
    ) + 
    # remove all the background stuff and wrap it with a dashed line
    theme_classic() + 
    theme(
      panel.background = element_rect(
        color = color, 
        linetype = "longdash"
      ), 
      axis.line = element_blank(), 
      axis.ticks = element_blank(), 
      axis.text.y = element_blank(), 
      axis.text.x = element_blank()
    )
}
my_custom_cor(dat, aes(P17_del_all,PTGS1,PTGS2))
my_custom_smooth <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_point(color = I("blue")) + 
    geom_smooth(method = "lm", color = I("black"), ...)
}
my_custom_smooth(iris, aes(Sepal.Length, Sepal.Width))

pm <- ggpairs(
  dat[,1:3], 
  upper = list(continuous = wrap(my_custom_cor,sizeRange = c(6,10))),
  lower = list(continuous = my_custom_smooth), 
  axisLabels = "show"
)
pm
ggsave("correlation.svg", plot=pm,width = 10, height = 8,dpi=1080)


##################################################
#Cor of PTGS genes with 17pGenes in TCGA-DLBCL####
#Cor of PTGS genes with 17pGenes in TCGA-DLBCL####
##################################################
sel <- read.csv("./TCGA_DLBCL_ALL_MERGE_FPKM.csv")
rownames(sel) <- sel$X
sel <- sel[,-1]
P17_gene_all <- c("SCO1","MYH3","MYH2","MYH1","MYH4","MYH8","MYH13","RPS27AP1","GAS7","RCVRN","GLP2R","DHRS7C","USP43","WDR16","STX8","NTN1","PIK3R5","PIK3R6","MFSD6L","SPDYE4","CCDC42","RPS26P53","MYH10","NDEL1","RNF222","RPL26","KRBA2","ODF4","ARHGEF15","RANGRF","PFAS","TRNAI24","TRNAS14","TRNAT9","TRNAP8","TRNAD17","TRNAG18","TRNAW4","AURKB","TRNAI23","TRNAT8","TRNAS23","TRNAW3","SNORD118","TMEM107","VAMP2","PER1","TRNAT22","TRNAS6","TRNAG2","TRNAR1","HES7","TRNAL2","TRNAQ1","TRNAK1","ALOXE3","ALOX12B","ALOX15B","GUCY2D","CNTROB","TRAPPC1","KCNAB3","SCARNA21","CHD3","CYB5D1","LSMD1","TMEM88","KDM6B","RPL29P2","DNAH2","EFNB3","WRAP53","TP53","ATP1B2","SAT2","SHBG","FXR2","SOX15","MPDU1","CD68","SNORA67","SNORD10","SNORA48","EIF4A1","SENP3","TNFSF13","TNFSF12-TNFSF13","TNFSF12","POLR2A","AMAC1L3","ZBTB4","CHRNB1","FGF11","TMEM102","C17orf74","SPEM1","NLGN2","PLSCR3","TNK1","TMEM95","KCTD11","ACAP1","NEURL4","GPS2","EIF5A","YBX2","SLC2A4","CLDN7","DULLARD","GABARAP","PHF23","DVL2","MIR324","ACADVL","DLG4","ASGR1","RPL7AP64","ASGR2","CLEC10A","SLC16A11","SLC16A13","BCL6B","ALOX12")
P17_all <- c()
for (i in P17_gene_all){
	tmp <- subset(sel,symbol==i)
	P17_all <- rbind(P17_all,tmp)
	print(paste(i, "is done",sep=" "))
}
rownames(P17_all) <- P17_all$symbol
P17_all <- P17_all[,-ncol(P17_all)]
a <- as.data.frame(colSums(P17_all)/nrow(P17_all))
a <- t(a)
rownames(a) <- "P17_del_all"

PTGS1 <- subset(sel,symbol=="PTGS1")
PTGS1 <- PTGS1[,-ncol(PTGS1)]
rownames(PTGS1) <- "PTGS1"
PTGS2 <- subset(sel,symbol=="PTGS2")
PTGS2 <- PTGS2[,-ncol(PTGS2)]
rownames(PTGS2) <- "PTGS2"
c <- rbind(a,PTGS1)
c <- rbind(c,PTGS2)
c <- as.data.frame(t(c))
c$patient_id <- as.factor(rownames(c))
dat <- c

PTGS1 <- subset(sel,symbol=="PTGS1")
PTGS1 <- PTGS1[,-ncol(PTGS1)]
rownames(PTGS1) <- "PTGS1"
PTGS2 <- subset(sel,symbol=="PTGS2")
PTGS2 <- PTGS2[,-ncol(PTGS2)]
rownames(PTGS2) <- "PTGS2"
c <- rbind(a,PTGS1)
c <- rbind(c,PTGS2)
c <- as.data.frame(t(c))
c$patient_id <- as.factor(rownames(c))
dat <- c

dat2 <- dat
dat2[,c(1:3)] <- log(dat2[,c(1:3)],2)

library(GGally)
my_custom_cor <- function(data, mapping, color = I("grey50"), sizeRange = c(1, 5), ...) {

  # get the x and y data to use the other code
  x <- GGally::eval_data_col(data, mapping$x)
  y <- GGally::eval_data_col(data, mapping$y)

  ct <- cor.test(x,y)
  sig <- symnum(
    ct$p.value, corr = FALSE, na = FALSE,
    cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
    symbols = c("***", "**", "*", ".", " ")
  )

  r <- unname(ct$estimate)
  rt <- format(r, digits=2)[1]

  # since we can't print it to get the strsize, just use the max size range
  cex <- max(sizeRange)

  # helper function to calculate a useable size
  percent_of_range <- function(percent, range) {
    percent * diff(range) + min(range, na.rm = TRUE)
  }

  # plot the cor value
  ggally_text(
    label = as.character(rt), 
    mapping = aes(),
    xP = 0.5, yP = 0.5, 
    size = I(percent_of_range(cex * abs(r), sizeRange)),
    color = color,
    ...
  ) + 
    # add the sig stars
    geom_text(
      aes_string(
        x = 0.8,
        y = 0.8
      ),
      label = sig, 
      size = I(cex),
      color = color,
      ...
    ) + 
    # remove all the background stuff and wrap it with a dashed line
    theme_classic() + 
    theme(
      panel.background = element_rect(
        color = color, 
        linetype = "longdash"
      ), 
      axis.line = element_blank(), 
      axis.ticks = element_blank(), 
      axis.text.y = element_blank(), 
      axis.text.x = element_blank()
    )
}
my_custom_cor(dat, aes(P17_del_all,PTGS1,PTGS2))
my_custom_smooth <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_point(color = I("blue")) + 
    geom_smooth(method = "lm", color = I("black"), ...)
}
my_custom_smooth(iris, aes(Sepal.Length, Sepal.Width))

pm <- ggpairs(
  dat[,1:3], 
  upper = list(continuous = wrap(my_custom_cor,sizeRange = c(6,10))),
  lower = list(continuous = my_custom_smooth), 
  axisLabels = "show"
)
pm


dat2 <- dat
dat2[,c(1:3)] <- log(dat2[,c(1:3)]+1,2)
dat2 <- subset(dat2,PTGS2<= 3)
my_custom_cor <- function(data, mapping, color = I("grey50"), sizeRange = c(1, 5), ...) {

  # get the x and y data to use the other code
  x <- GGally::eval_data_col(data, mapping$x)
  y <- GGally::eval_data_col(data, mapping$y)

  ct <- cor.test(x,y)
  sig <- symnum(
    ct$p.value, corr = FALSE, na = FALSE,
    cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
    symbols = c("***", "**", "*", ".", " ")
  )

  r <- unname(ct$estimate)
  rt <- format(r, digits=2)[1]

  # since we can't print it to get the strsize, just use the max size range
  cex <- max(sizeRange)

  # helper function to calculate a useable size
  percent_of_range <- function(percent, range) {
    percent * diff(range) + min(range, na.rm = TRUE)
  }

  # plot the cor value
  ggally_text(
    label = as.character(rt), 
    mapping = aes(),
    xP = 0.5, yP = 0.5, 
    size = I(percent_of_range(cex * abs(r), sizeRange)),
    color = color,
    ...
  ) + 
    # add the sig stars
    geom_text(
      aes_string(
        x = 0.8,
        y = 0.8
      ),
      label = sig, 
      size = I(cex),
      color = color,
      ...
    ) + 
    # remove all the background stuff and wrap it with a dashed line
    theme_classic() + 
    theme(
      panel.background = element_rect(
        color = color, 
        linetype = "longdash"
      ), 
      axis.line = element_blank(), 
      axis.ticks = element_blank(), 
      axis.text.y = element_blank(), 
      axis.text.x = element_blank()
    )
}
my_custom_cor(dat2, aes(P17_del_all,PTGS1,PTGS2))
my_custom_smooth <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_point(color = I("blue")) + 
    geom_smooth(method = "lm", color = I("black"), ...)
}
my_custom_smooth(iris, aes(Sepal.Length, Sepal.Width))

pm <- ggpairs(
  dat2[,1:3], 
  upper = list(continuous = wrap(my_custom_cor,sizeRange = c(6,10))),
  lower = list(continuous = my_custom_smooth), 
  axisLabels = "show"
)
pm
ggsave("DLBCL_PTGS2.svg", plot=pm,width = 10, height = 8,dpi=1080)




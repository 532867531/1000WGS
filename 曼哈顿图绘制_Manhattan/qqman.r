library(qqman)
library(calibrate)
qqman_<-function (x, chr = "CHR", bp = "BP", p = "P", snp = "SNP", col = "COLOR", chrlabs = NULL, suggestiveline = -log10(1e-05), 
          genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE, 
          annotatePval = NULL, annotateTop = TRUE) 
{
  CHR = BP = P = index = NULL
  if (!(chr %in% names(x))) 
    stop(paste("Column", chr, "not found!"))
  if (!(bp %in% names(x))) 
    stop(paste("Column", bp, "not found!"))
  if (!(p %in% names(x))) 
    stop(paste("Column", p, "not found!"))
  if (!(snp %in% names(x))) 
    warning(paste("No SNP column found. OK unless you're trying to highlight."))
  if (!is.numeric(x[[chr]])) 
    stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
  if (!is.numeric(x[[bp]])) 
    stop(paste(bp, "column should be numeric."))
  if (!is.numeric(x[[p]])) 
    stop(paste(p, "column should be numeric."))
  d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]],col=x[[col]],LABEL=x[['LABEL']])
  if (!is.null(x[[snp]])) 
    d = transform(d, SNP = x[[snp]])
  d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
  d <- d[order(d$CHR, d$BP), ]
  if (logp) {
    d$logp <- -log10(d$P)
  }
  else {
    d$logp <- d$P
  }
  d$pos = NA
  d$index = NA
  ind = 0
  for (i in unique(d$CHR)) {
    ind = ind + 1
    d[d$CHR == i, ]$index = ind
  }
  nchr = length(unique(d$CHR))
  if (nchr == 1) {
    d$pos = d$BP
    ticks = floor(length(d$pos))/2 + 1
    xlabel = paste("Chromosome", unique(d$CHR), "position")
    labs = ticks
  }
  else {
    lastbase = 0
    ticks = NULL
    for (i in unique(d$index)) {
      if (i == 1) {
        d[d$index == i, ]$pos = d[d$index == i, ]$BP
      }
      else {
        lastbase = lastbase + tail(subset(d, index == 
                                            i - 1)$BP, 1)
        d[d$index == i, ]$pos = d[d$index == i, ]$BP + 
          lastbase
      }
      ticks = c(ticks, (min(d[d$index == i, ]$pos) + max(d[d$index == 
                                                             i, ]$pos))/2 + 1)
    }
    xlabel = "Chromosome"
    labs <- unique(d$CHR)
  }
  xmax = ceiling(max(d$pos) * 1.03)
  xmin = floor(max(d$pos) * -0.03)
  def_args <- list(xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", 
                   las = 1, pch = 20, xlim = c(xmin, xmax), ylim = c(0, 
                                                                     ceiling(max(d$logp))), xlab = xlabel, ylab = expression(-log[10](italic(p))))
  dotargs <- list()
  do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% 
                                            names(dotargs)]))
  if (!is.null(chrlabs)) {
    if (is.character(chrlabs)) {
      if (length(chrlabs) == length(labs)) {
        labs <- chrlabs
      }
      else {
        warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
      }
    }
    else {
      warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
    }
  }
  if (nchr == 1) {
    axis(1)
  }
  else {
    axis(1, at = ticks, labels = labs)
  }
  col1 = rep(c("gray","black"), max(d$CHR))
  #添加空白数据点
  # for(onechr in unique(d$index)){
    # temp=d[NA,]
    # temp=na.omit(temp)
    # ones=runif(n = 1000,min = max(d[d$index == onechr,"pos"]),max = max(d[d$index == onechr,"pos"])+3000000)
    # for(rowindex in c(1:length(ones))){
      # temp[rowindex,"pos"]=ones[rowindex]
    # }
    # temp$logp=0
    # temp$index=onechr
    # print(head(temp))
    # d=rbind(d,temp)
  # }
  if (nchr == 1) {
    with(d, points(pos, logp, pch = 20, col = col1[1]))
  }
  else {
    icol = 1
    for (i in unique(d$index)) {
      with(d[d$index == unique(d$index)[i], ], points(pos, 
                                                      logp, col = col1[icol], pch = 20));
      
      icol = icol + 1
    }
  }
  
  if (suggestiveline) 
    abline(h = suggestiveline, col = "blue")
  if (genomewideline) 
    abline(h = genomewideline, col = "red")
  if (!is.null(highlight)) {
    if (any(!(highlight %in% d$SNP))) 
      warning("You're trying to highlight SNPs that don't exist in your results.")
    d.highlight = d[intersect(which(d$SNP %in% highlight),which(d$col == "green3")), ]
	#d.highlight = d[which(d$SNP %in% highlight),]
    #with(d.highlight, points(pos, logp, col = "green3", pch = 4
    message("d.highlight")
    print(head(d.highlight))
    with(d.highlight, 
	{
		points(pos, logp, col = "green3", pch = 4)
		#message(col)
	}
	)
  }
  if (!is.null(annotatePval)) {
    topHits = subset(d, P <= annotatePval)
    message("topHits")
    print(head(topHits))
    par(xpd = TRUE)
    if (annotateTop == FALSE) {
      with(subset(d, P <= annotatePval), textxy(pos, -log10(P), 
                                                offset = 0.625, labs = topHits$LABEL, cex = 0.45), 
      )
    }
    else {
      topHits <- topHits[order(topHits$P), ]
      topSNPs <- NULL
      for (i in unique(topHits$CHR)) {
        chrSNPs <- topHits[topHits$CHR == i, ]
        topSNPs <- rbind(topSNPs, chrSNPs[1, ])
      }
      textxy(topSNPs$pos, -log10(topSNPs$P), offset = 0.625, 
             labs = topSNPs$LABEL, cex = 0.5)
    }
  }
  par(xpd = FALSE)
}
message("qqman_")
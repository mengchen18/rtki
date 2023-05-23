# tumcolor

tumcolor <- list(
  blue = "#0065BD",
  gray = "#999999",
  lightgray = "#dad7cb",
  blue_scale3 = c("#0065BD", "#005293", "#003359"),
  gray_scale3 = c("#CCCCCC", "#808080", "#333333"),
  c4 = c("#DAD7CB", "#E37222", "#A2AD00", "#64A0C8"),
  c8  = c("#3070B3", "#FED702", "#F7811E", "#B55CA5", "#8F81EA", "#EA7237", "#9FBA36", "#6A757E")
)

tcolor <- function(x, alpha = 100) {
  n <- col2rgb(x)
  rgb(n[1], n[2], n[3], alpha = alpha, maxColorValue = 255)
}

tcolorl <- function(x, ...) 
  sapply(x, tcolor, ...)


##############

rescale <- function(x, min = 2, max = 4) {
  x <- x-min(x)
  s <- max - min
  s * (x/max(x)) + min
}


drawSeqs <- function(x, dat, file = NULL, pY = NULL) {
  
  # on exit
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(xpd = oldpar$xpd))
  
  # generate plot area
  rtkname <- as.character(x$summary$proteinName[1])
  nc <- nchar(x$seq)
  par(mar = c(3, 8, 4, 3))
  plot(x = 0, col = NA, xlim = c(0, 2), ylim = c(nc+1, 0), xlab = "", ylab = "", axes = F, main = rtkname)
  
  # regions
  xf <- x$feature
  xf <- xf[xf$type == "domain", ]
  
  # plot1
  csol <- distinctColorPalette(k = min(nrow(xf), 8))
  cc <- col2rgb(csol)
  cc <- rgb(cc[1, ], cc[2, ], cc[3, ], alpha = 150, maxColorValue = 255)
  
  #######################################################
  i <- which(dat$rtk == rtkname)
  d0 <- dat[i, ]
  d0$pYs <- sub(paste0(rtkname, "_"), "", d0$pY.id)
  d0$pos <- as.integer( stringr::str_extract(d0$pYs, "\\d+"))
  d0 <- d0[order(d0$pos), ]
  
  tab <- table(d0$pYs)
  pos <- as.integer( stringr::str_extract(names(tab), "\\d+"))
  ord <- order( pos, decreasing = FALSE)
  tab <- tab[ord]
  pos <- pos[ord]
  
  cx <- rep(0.9, length(tab))
  cy <- seq(1, nc, length.out = length(tab))
  names(cy) <- names(tab)
  
  #######################################################
  ipos <- !names(tab) %in% d0$pYs[d0$pY.type == "IF"]
  pos <- pos[ipos]
  segments(x0 = 0, y0 = pos, x1 = 0.2, y1 = pos, col = tcolor(tumcolor$blue, alpha = 100))
  segments(x0 = 0.2, y0 = pos, x1 = 0.9, y1 = cy[ipos], col = tcolor(tumcolor$blue, alpha = 100))
  
  
  #######################################################
  bp <- unique(d0$bp.gene.name)
  tab2 <- table(d0$bp.gene.name)[bp]
  rx <- rep(2, length(tab2))
  ry <- seq(2, nc, length.out = length(tab2))
  names(ry) <- names(tab2)
  
  #########
  d00 <- d0
  if (!is.null(pY))
    d00 <- d0[d0$pY.id2 %in% pY, , drop = FALSE]
  
  segments(x0 = 1.2, y0 = cy[match(d00$pYs, names(cy))], col = "gray",
           x1 = 2, y1 = ry[match(d00$bp.gene.name, names(ry))])
  
  tcol <- c("WT" = tcolor(tumcolor$c4[4]), 
            "PC" = tcolor(tumcolor$c4[4]), 
            "HSM" = tcolor(tumcolor$c4[2]), 
            "MGY" = tcolor(tumcolor$c4[3]), 
            "IF" = tcolor(tumcolor$c4[1]))
  
  typecol <-  unique(d0[, c("pYs", "pY.type")])
  tcol2 <- tcol[typecol$pY.type]
  names(tcol2) <- typecol$pYs
  
  ############# draw boxes and circles ###################
  # 1 - rect
  rect(xleft = -0.05, xright = .05, ybottom = 1, ytop = nc, col = "gray", border = "gray")
  # 2 - domain rect
  for (i1 in seq_len(nrow(xf))) {
    f1 <- xf[i1, ]
    rect(xleft = -0.1, xright = .1, ybottom = f1$begin, ytop = f1$end, col = cc[i1], border = "gray")
  }
  atm <- (xf$begin + xf$end)/2
  mtext(text = xf$description, side = 2, at = atm, las = 2)
  # 3 - pYs
  points(x = cx,  y = cy, cex = rescale(tab), las = 2, pch = 19, col = tcol2)
  text(cx, cy, labels = tab, cex = 0.8)
  text(cx+0.01, cy, labels = names(tab), cex = 0.9, pos = 4)
  # 4 - binder
  points(x = rx,  y = ry, cex = rescale(tab2), las = 2, pch = 19, col = tcolor(tumcolor$c8[5]))
  text(rx, ry, labels = tab2, cex = 0.8)
  mtext(side = 4, at = ry+0.01, text = names(tab2), cex = 0.9, las = 2)
  
  ############## labels and legends ###################
  mtext(side = 3, at = c(0, 1, 2), 
        text = c("Sequence", "pY site", "Binder"), cex = 0.9, las = 1, font = 2)
  text(0.1, y = c(0, nc), c("AA 0", paste0("AA ", nc)), font = 3, pos = 4)
  
  par(xpd = TRUE)
  tt <- tcol[-2]
  legend(x = 0.4, y = nc*1.03, col = tt, legend = names(tt), pch = 15, bty = "n", pt.cex = 2, ncol = 4, title = "pY type")
  d00
}

########
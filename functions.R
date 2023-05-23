drawSeqPrep <- function(x, dat, pY = NULL) {

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
  
  d00 <- d0
  if (!is.null(pY))
    d00 <- d0[d0$pY.id2 %in% pY, , drop = FALSE]
  list(tab = tab, df1 = d0, df2 = d00, cy = cy, cx = cx, pos = pos, nc = nc, xf = xf, cc = cc)
  }

drawSeq <- function(x) {
  # on exit
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(xpd = oldpar$xpd))
  
  d0 <- x$df1
  d00 <- x$df2
  cy <- x$cy
  cx <- x$cx
  tab <- x$tab
  pos <- x$pos
  nc <- x$nc
  xf <- x$xf
  cc <- x$cc
  
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
  
  if (nrow(d00) > 0) {
    segments(x0 = 1.2, y0 = cy[match(d00$pYs, names(cy))], col = "gray",
             x1 = 2, y1 = ry[match(d00$bp.gene.name, names(ry))])
    }
  
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
}

########



formatDT <- function(tab, sel = c("single", "multiple")[1], pageLength = 15) {
  
  dt <- DT::datatable(
    na2char( tab ),
    selection =  sel,
    rownames = FALSE,
    filter = "top",
    class="table-bordered compact nowrap",
    options = list(
      scrollX = TRUE, pageLength = pageLength, dom = 'tip', columnDefs = list(list(
        targets = unname(which(sapply(tab, inherits, c('factor', "character"))))-1,
        render = DT::JS(
          "function(data, type, row, meta) {",
          "return type === 'display' && data.length > 30 ?",
          "'<span title=\"' + data + '\">' + data.substr(0, 30) + '...</span>' : data;",
          "}")
      )))
  )
  DT::formatStyle(dt, columns = (1:ncol(tab))-1, fontSize = '90%')
}

na2char <- function(x) {
  cc <- sapply(x, inherits, c("character", "factor"))
  for (i in which(cc)) {
    if (!is.character(x[, i]))
      x[, i] <- as.character(x[, i])
    ir <- is.na(x[, i])
    x[ir, i] <- ""
  }
  x
}

formatDTScrollY <- function(tab, sel = c("single", "multiple")[1], height = 500) {
  dt <- DT::datatable( 
    na2char( tab ),
    extensions = 'Scroller',
    selection =  sel,
    rownames = FALSE,
    filter = "top",
    class="table-bordered compact nowrap",
    options = list(
      scrollX = TRUE, scrollY = height, dom = 'ti', scroller = TRUE, 
      pageLength = nrow(tab), columnDefs = list(list(
        targets = unname(which(sapply(tab, inherits, c('factor', "character"))))-1,
        render = DT::JS(
          "function(data, type, row, meta) {",
          "return type === 'display' && data.length > 30 ?",
          "'<span title=\"' + data + '\">' + data.substr(0, 30) + '...</span>' : data;",
          "}")
      ))
    )
  )
  DT::formatStyle(dt, columns = 1:ncol(tab), fontSize = '90%')
}


hm <- function(tab) {
  
  dd <- dcast(data = tab, bp.gene.name ~ pY.id, value.var = "af.filter.batch", fun.aggregate = max)
  dd[-1] <- lapply(dd[-1], function(x) {
    x[is.infinite(x)] <- NA
    x
  })
  dd$domain <- c("", "+")[as.integer(dd$bp.gene.name %in% pl_domain) + 1]
  
  d1 <- dd[dd$domain == "+", ]
  d2 <- dd[dd$domain == "", ]
  
  cc <- rbind(d1, d2)
  ii <- is.na(cc)
  cc[ii] <- 0.6
  
  mm <- apply(cc[, sapply(cc, is.numeric), drop = FALSE], 2, c)
  if (!is.matrix(mm)) {
    cn0 <- names(mm)
    mm <- matrix(mm, nrow = 1)
    colnames(mm) <- cn0
  }
    
  rownames(mm) <- paste(cc$bp.gene.name, cc$domain)
  mm <- t(mm)
  
  sc <- mm
  sc[mm == 0.6] <- "Enrichment factor < 0.6"
  sc[mm > 0.6] <- paste0("Enrichment factor = ", signif(mm[mm > 0.6], 3))
  
  heatmaply(
    mm, scale = "none", 
    color = c("white", colorRampPalette(brewer.pal(n = 7, name ="Blues"))(100)), 
    breaks = c(0.5, seq(0.60001, 3, length.out = 99)),
    Rowv = nrow(mm) > 2,
    Colv = ncol(mm) > 2, 
    custom_hovertext = sc
  )
}

run_plot_heatmap_LD <- function (data.cross = NULL, data.tibble =NULL,  heterozygotes = FALSE){

  nchr.max <- nchr (data.cross)

  list.chr <- 1:nchr.max
  #filt.chr = 1


  lapply(list.chr, function (filt.chr){

    print (filt.chr)
    crossobj <- subset(data.cross, chr=filt.chr)
    chr=filt.chr

    data <- NULL
    for (i in 1:length(chr)) {
      a <- paste("crossobj$geno$'", chr[i], "'$data", sep = "")
      p1 <- eval(parse(text = a))
      data <- cbind(data, p1)
    }

    if (heterozygotes == "FALSE") {
      data[data == 1] <- "A/A"
      data[data == 2] <- "B/B"
    }

    if (heterozygotes == "TRUE") {
      data[data == 1] <- "A/A"
      data[data == 2] <- "B/B"
      data[data == 3] <- "A/B"
    }

    genos <- genotype(data[, 1])
    for (i in 2:dim(data)[2]) {
      g <- genotype(data[, i])
      genos <- data.frame(genos, g)
    }

    #### add marker names calculate LD
    ld <- LD(genos)

    # plot LD heatmap
    plot.hm <- LDheatmap(ld$"R^2", title=str_c("Pairwise LD_Chr.",filt.chr),
                         color = colorRampPalette(c("red4", "red","orangered", "orange","yellow1",  "blue4"))(60))


    LD.soja.cross.matrix <- plot.hm$LDmatrix


    write.table (LD.soja.cross.matrix ,
                 file = str_c("./Data/procdata/LD.soja.cross.matrix.",filt.chr,".txt"),
                 append = FALSE, quote = TRUE, sep = ",",
                 eol = "\n", na = "NA", dec = ".", row.names = TRUE,
                 col.names = TRUE)

    if(filt.chr < 10 ) {

      soja.cross.tibble_snp <- data.tibble  %>%
        dplyr::select (starts_with(str_c("S0", filt.chr)))
    }

    if(filt.chr >= 10 ) {

      soja.cross.tibble_snp <- data.tibble  %>%
        dplyr::select (starts_with(str_c("S", filt.chr)))

    }

    #soja.cross.chr.tibble_snp.1 <- soja.cross.tibble_snp [-c(1,2),]


    soja.cross.tibble_bp <- soja.cross.tibble_snp [2,]

    mkr.chr <- colnames (soja.cross.tibble_snp)

    x.mk <- data.frame (mrks = mkr.chr)

    pos.chr    <- soja.cross.tibble_snp [2,]
    pos.chr.1  <- data.frame (pos = pos.chr[1,])

    pos.chr.2  <- as.numeric(t(pos.chr.1))

    soja.cross.tibble_bp.1 <- x.mk %>%
      dplyr::mutate (pos=pos.chr.2)

    ###
    dt <- soja.cross.tibble_bp.1
    list.pos <- (unique (dt$pos))
    list.mrks <- (unique (dt$mrks))

    Z <- matrix (0, nrow=length(list.pos), ncol=length(list.pos))
    Z <- as.data.frame(Z)
    colnames(Z) <-list.mrks
    rownames(Z) <-list.mrks

    #### aca empieza la distancia en bp
    x.bp <- lapply (list.pos, function (filtro.x1) {
      lapply (list.pos, function (filtro.x2) {
        dt.x1 <- dt %>%
          dplyr::filter (pos == filtro.x1)
        x1    <- dt.x1 [,2]
        id.x1 <- dt.x1 [,1]
        dt.x2 <- dt %>%
          dplyr::filter (pos == filtro.x2)
        x2    <- dt.x2 [,2]
        id.x2 <- dt.x2 [,1]
        dt.z <- abs (x1 - x2)
        Z ["id.x1","id.x2"] <- dt.z
      })
    })

    names(x.bp) <- list.mrks
    XX.x.bp <- as.data.frame (do.call (cbind, x.bp))

    dt.LD.decay <- lapply (list.mrks, function (filtro) {

      XX.x.bp.1 <- XX.x.bp %>%
        dplyr::select (filtro)
      delta.bp <- unlist (XX.x.bp.1 , use.names=FALSE)

      colnames (LD.soja.cross.matrix ) <- list.mrks
      rownames (LD.soja.cross.matrix) <- list.mrks
      LD_soja.cross <-  as_tibble(LD.soja.cross.matrix)
      LD_soja.cross <- LD_soja.cross %>%
        dplyr::mutate (mrk.id =list.mrks)%>%
        dplyr::select (mrk.id, everything())

      id <- colnames (XX.x.bp.1)

      LD_soja.cross.1 <- LD_soja.cross%>%
        dplyr::filter (mrk.id==id)

      id.gather <- colnames (LD_soja.cross.1) [-1]

      LD_soja.cross.2 <- LD_soja.cross.1 %>%
        gather (id.gather, key="mrk.2" , value = "R2") %>%
        dplyr::rename (mrk.1 = mrk.id)

      LD_soja.cross.2$mrk.1 <- as.character(LD_soja.cross.2$mrk.1)


      colnames (XX.x.bp.1) == unique ( LD_soja.cross.2$mrk.1)

      LD_soja.cross.3 <- LD_soja.cross.2 %>%
        dplyr::mutate (delta.bp = delta.bp )

      print (LD_soja.cross.3)



      return ( LD_soja.cross.3)

    })

    df.LD.decay <- as.data.frame (do.call (rbind, dt.LD.decay))

    df.LD.decay <- df.LD.decay  %>%
      dplyr::mutate (chrom = str_c ("chr_", filt.chr))%>%
      dplyr::arrange(delta.bp)%>%
      dplyr::mutate (delta.bp.1 = delta.bp/1e6)

    write_delim (df.LD.decay  , path=str_c("./Data/procdata/plot.LD.decay.chr", filt.chr,".txt"),
                 delim = ",", na = "NA")

    png (filename = str_c("./Figures/plot.r2.decay.chr_", filt.chr,".png"),
         width = 480, height = 480, units = "px", pointsize = 12,
         bg = "white", res = NA)

    plot (x = df.LD.decay$delta.bp.1 , y = df.LD.decay$R2,
          main=str_c("LD.decay.", filt.chr),
          pch = 20,
          type ="n",
          xaxt="none",
          yaxt="none",
          axes = F,
          xlim = c(0, max (df.LD.decay$delta.bp.1)),
          ylim = c(0, max (df.LD.decay$R2, na.rm = TRUE)),
          ylab = expression(LD ~ (r^2)),
          xlab = expression(Distance ~ (Mb)))

    axis(side = 2, las = 1)
    x2 <- max (df.LD.decay$delta.bp.1, na.rm = TRUE)
    axis (side=1,at=seq(0,x2,10),las = 1)


    points (df.LD.decay$delta.bp.1,df.LD.decay$R2,
            pch = 20, cex=1.5, col="gray28")
    box()
    dev.off()

    plot (x = df.LD.decay$delta.bp.1 , y = df.LD.decay$R2,
          main=str_c("LD.decay.", filt.chr),
          pch = 20,
          type ="n",
          xaxt="none",
          yaxt="none",
          axes = F,
          xlim = c(0, max (df.LD.decay$delta.bp.1)),
          ylim = c(0, max (df.LD.decay$R2, na.rm = TRUE)),
          ylab = expression(LD ~ (r^2)),
          xlab = expression(Distance ~ (Mb)))

    axis(side = 2, las = 1)
    x2 <- max (df.LD.decay$delta.bp.1, na.rm = TRUE)
    axis (side=1,at=seq(0,x2,10),las = 1)


    points (df.LD.decay$delta.bp.1,df.LD.decay$R2,
            pch = 20, cex=1.5, col="gray28")
    box()

  })
}

# hay que editar

run_plot_heatmap_LD <- function (id.cross = NULL ,
                                 data.cross = NULL,
                                 heterozygotes = FALSE,
                                 distance.unit = NULL) {

  # Nota: se debe crear la estructura del directorio

  print (str_c("Se encontraron " , nchr (data.cross), " grupos de ligamientos")) # verifico el numero de cromosomas

  if   ( distance.unit != "cM" &  distance.unit != "bp" &  distance.unit != "kb" & distance.unit != "Mb") {

    stop (str_c ("Debe definir una unidad de distancia valida"))


  }

  if (is.null (distance.unit)) {

    stop (str_c ("Falta definir unidad de distancia"))

  }

  if   (distance.unit == "bp" | distance.unit == "Kb" | distance.unit == "Mb"  ) {

    print (str_c ("El mapa es un mapa fisico"))

  }

  if   ( distance.unit == "cM") {
    print (str_c ("El mapa es un mapa genetico"))

  }

  #  numero maximo de grupos de ligamientos

  nchr.max <- nchr (data.cross)

  list.chr <- 1:nchr.max

  ### este es temporal para corre el lapply

  #filt.chr = 1
  ###########

  lapply(list.chr, function (filt.chr){

    print (str_c("Estimando LD LG= " ,filt.chr))

    ## hago un subset de cada cromosoma

    crossobj <- subset(data.cross, chr=filt.chr)

    chr <- filt.chr

    #  Esta seccion es para formatear la matriz para para LDheatmap.
    #  Nota: ver si con los nuevos formatos cambia o no

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
    ### Esta es la funcion que calcula el LD, verificar la documentacion
    ##
    ##############

    #######################3
    # Nota : hay que agregar barra de progreso

    ld <- LD (genos)

    map.cross <- pull.map ( crossobj, as.table = TRUE)

    # plot LD heatmap
    # Nota: las figuras se podrian guardar solas como png.

    if   ( distance.unit == "bp" | distance.unit == "Kb" | distance.unit == "Mb" ) {

      plot.hm <- LDheatmap ( ld$"R^2",
                             genetic.distances= map.cross$pos,
                             distances="physical",
                             title=str_c(id.cross, "_Pairwise LD_LG.",filt.chr),
                             color = colorRampPalette(c("red4", "red","orangered", "orange","yellow1",  "blue4"))(60))

    }

    if   ( distance.unit == "cM") {

      plot.hm <- LDheatmap ( ld$"R^2", genetic.distances= map.cross$pos,
                             distances="genetic",
                             title=str_c(id.cross, "Pairwise LD_LG.", filt.chr),
                             color = colorRampPalette(c("red4", "red","orangered", "orange","yellow1",  "blue4"))(60))
    }

    ### esta es la matriz que sale de la funcion LD

    LD.cross.matrix <- plot.hm$LDmatrix

    write.table (LD.cross.matrix ,
                 file = str_c("./Data/procdata/",id.cross, "_LD.cross.matrix_", filt.chr,".txt"),
                 append = FALSE, quote = TRUE, sep = ",",
                 eol = "\n", na = "NA", dec = ".", row.names = TRUE,
                 col.names = TRUE)

    ### genero el df de los marcadores y sus distancias

    cross.tibble_dist <- map.cross %>%
      dplyr::mutate ( mrks  = rownames(map.cross)) %>%
      dplyr::select ( mrks, pos)


    cross.tibble_dist <- as_tibble (cross.tibble_dist)

    dt <-  cross.tibble_dist
    list.pos <- (unique (dt$pos))
    list.mrks <- (unique (dt$mrks))

    #### aca empieza la distancia

    print (str_c("Estimando diff.dist LG= " ,filt.chr))



    dt.diff.dist <- bind_rows (lapply (list.pos, function (filtro.x1) {

      #filtro.x1 =  277230
      #print (filtro.x1)

      dt.dist.2 <- bind_rows ( lapply (list.pos, function (filtro.x2) {

        #filtro.x2 =  2415604
        #print (filtro.x2)

        dt.x1 <- dt %>%
          dplyr::filter (pos == filtro.x1)

        dt.x2 <- dt %>%
          dplyr::filter (pos == filtro.x2)

        dt.z <- data.frame (dt.x1,  dt.x2)

        dt.z <- dt.z  %>%
          dplyr::mutate (diff.dist = abs (pos - pos.1)) %>%
          dplyr::select ("mrks", "mrks.1", "diff.dist" ) %>%
          dplyr::rename (mrk.1 = mrks) %>%
          dplyr::rename (mrk.2 = mrks.1)

        #print (dt.z)

        #return (dt.z)

      }))

    }))

    start.time <- Sys.time()

    df.LD.decay <- bind_rows ( lapply (list.mrks, function (filt.mrk) {

      #filt.mrk= "JHI-Hv50k-2016-270"

      #print (filt.mrk)

      dt.diff.dist.mrk <-  dt.diff.dist %>%
        dplyr::filter (mrk.1 == filt.mrk)

      colnames (LD.cross.matrix ) <- list.mrks
      rownames (LD.cross.matrix)  <- list.mrks

      LD.cross <-  as_tibble (LD.cross.matrix)
      LD.cross <- LD.cross %>%
        dplyr::mutate (mrk.id = list.mrks)%>%
        dplyr::select (mrk.id, everything())

      LD.mrk.1 <- LD.cross %>%
        dplyr::filter (mrk.id== filt.mrk)

      id.gather <- colnames (LD.mrk.1) [-1]

      LD.cross.2 <- LD.mrk.1 %>%
        gather (all_of (id.gather), key="mrk.2" , value = "R2") %>%
        dplyr::rename (mrk.1 = mrk.id)

      LD.dist.cross.3 <- LD.cross.2 %>%
        dplyr::inner_join ( dt.diff.dist.mrk ,  LD.cross.2, by= c("mrk.1", "mrk.2"))

      #return (LD.dist.cross.3)

    }))




    df.LD.decay <- df.LD.decay  %>%
      dplyr::mutate (LG = str_c ("lg.", filt.chr)) %>%
      dplyr::arrange (diff.dist)

    ### verificar estas medidas

    if   ( distance.unit == "bp" ) {

      df.LD.decay <- df.LD.decay  %>%
        dplyr::mutate (diff.dist = (diff.dist*1e-6))

    }


    if   ( distance.unit == "Kb" ) {

      df.LD.decay <- df.LD.decay  %>%
        dplyr::mutate (diff.dist = diff.dist/1e6)

    }

    write_delim (df.LD.decay  , file =str_c("./Data/procdata/", id.cross, "LD.decay.LG", filt.chr,".txt"),
                 delim = ",", na = "NA")

    if   ( distance.unit == "bp" | distance.unit == "Kb" | distance.unit == "bp" ) {

      png (filename = str_c("./Figures/",id.cross,".LD.decay.LG_", filt.chr,".png"),
           width = 480, height = 480, units = "px", pointsize = 12,
           bg = "white", res = NA)

      plot (x = df.LD.decay$diff.dist , y = df.LD.decay$R2,
            main=str_c(id.cross,".LD.decay.LG_", filt.chr),
            pch = 20,
            type ="n",
            xaxt="none",
            yaxt="none",
            axes = F,
            xlim = c(0, max (df.LD.decay$diff.dist)),
            ylim = c(0, max (df.LD.decay$R2, na.rm = TRUE)),
            ylab = expression(LD ~ (r^2)),
            xlab = expression(Distance ~ (Mb)))

      axis(side = 2, las = 1)
      x2 <- max (df.LD.decay$diff.dist, na.rm = TRUE)
      axis (side=1,at=seq(0,x2,10),las = 1)


      points (df.LD.decay$diff.dist,df.LD.decay$R2,
              pch = 20, cex=1.5, col="gray28")
      box()
      dev.off()

      plot (x = df.LD.decay$diff.dist , y = df.LD.decay$R2,
            main=str_c(id.cross,".LD.decay.LG_", filt.chr),
            pch = 20,
            type ="n",
            xaxt="none",
            yaxt="none",
            axes = F,
            xlim = c(0, max (df.LD.decay$diff.dist)),
            ylim = c(0, max (df.LD.decay$R2, na.rm = TRUE)),
            ylab = expression(LD ~ (r^2)),
            xlab = expression(Distance ~ (Mb)))

      axis(side = 2, las = 1)
      x2 <- max (df.LD.decay$diff.dist, na.rm = TRUE)
      axis (side=1,at=seq(0,x2,10),las = 1)

      points (df.LD.decay$diff.dist, df.LD.decay$R2,
              pch = 20, cex=1.5, col="gray28")
      box()
    }

    #if   ( distance.unit == "cM" ) {

    #png (filename = str_c("./Figures/",id.cross,".LD.decay.LG_", filt.chr,".png"),
    #     width = 480, height = 480, units = "px", pointsize = 12,
    #    bg = "white", res = NA)

    #  plot (x = df.LD.decay$diff.dist , y = df.LD.decay$R2,
    #       main=str_c(id.cross,".LD.decay.LG_", filt.chr),
    #        pch = 20,
    #       type ="n",
    #      xaxt="none",
    #     yaxt="none",
    #    axes = F,
    #   xlim = c(0, max (df.LD.decay$diff.dist)),
    #  ylim = c(0, max (df.LD.decay$R2, na.rm = TRUE)),
    # ylab = expression(LD ~ (r^2)),
    #xlab = expression(Distance ~ (cM)))

    # axis(side = 2, las = 1)
    #x2 <- max (df.LD.decay$diff.dist, na.rm = TRUE)
    #axis (side=1,at=seq(0,x2,1),las = 1)


    #  points (df.LD.decay$diff.dist, df.LD.decay$R2,
    #pch = 20, cex=1.5, col="gray28")
    # box()
    #dev.off()

    #plot (x = df.LD.decay$diff.dist, y = df.LD.decay$R2,
    #     main=str_c(id.cross,".LD.decay.LG_", filt.chr),
    #pch = 20,
    #type ="n",
    #xaxt="none",
    #yaxt="none",
    #axes = F,
    #xlim = c(0, max (df.LD.decay$diff.dist)),
    #ylim = c(0, max (df.LD.decay$R2, na.rm = TRUE)),
    #ylab = expression(LD ~ (r^2)),
    #xlab = expression(Distance ~ (cM)))

    #  axis(side = 2, las = 1)
    #x2 <- max (df.LD.decay$diff.dist, na.rm = TRUE)
    #axis (side=1,at=seq(0,x2,1),las = 1)


    #points (df.LD.decay$diff.dist,df.LD.decay$R2,
    #pch = 20, cex=1.5, col="gray28")
    #box()
    #}#

  })

}

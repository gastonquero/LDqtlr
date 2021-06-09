
run_freq_decay <- function (id.cross =NULL, data, chr=NULL, ini = NULL, l1= 0.25, l2=0.5, l3=0.75, keep.Mb =NULL,
                            prob.HLMK = NULL,
                            prob.HUMK = NULL)
{

  # el primer index p
  #index.chrom <- unique (data$chr)

  # control de la clases de los objetos
  assert_is_data.frame (data)
  #assert_is_numeric (chr)
  #assert_is_numeric(ini)
  assert_is_numeric(l1)
  assert_is_numeric(l2)
  assert_is_numeric(l3)


  # if (any (is_non_positive(ini), na.rm = TRUE)) {
  # Throw an error
  # stop ("ini contains non-positive values, so no puede caminar.")
  #}

  if (any (is_non_positive(c(l1, l2, l3)), na.rm = TRUE)) {
    # Throw an error
    stop ("limit contains non-positive values, so no se puede calcular.")
  }


  if (l1 > 1) {
    # Throw an error
    stop ("l1 contains valor mayor a 1, so no se puede calcular.")
  }


  if (l2 > 1) {
    # Throw an error
    stop ("l2 contains valor mayor a 1, so no se puede calcular.")
  }

  if (l3 > 1) {
    # Throw an error
    stop ("l2 contains valor mayor a 1, so no se puede calcular.")
  }


  index.chrom <- unique (df.geno.crom$chrom)

  #filt.chr =1

  dt.all.chrom <- bind_rows (lapply (index.chrom, function (filt.chr) {

    print (filt.chr)

    x1 <- data %>%
      dplyr::filter (chrom == filt.chr)

    ### reeeplazr NA por 1
    x1.na <- x1 %>%
      dplyr:::filter (R2 != "NA")

    #max_delta_bp <- max(x1.na$diff.dist)
    max_delta_Mb <- max(x1.na$diff.dist)

    max_delta_Mb.r <- round(max_delta_Mb + 1,0)

    #index.Mb <- seq (1, (round(max_delta_Mb + 1,0)), 1)

    index.Mb <- seq (keep.Mb, max_delta_Mb.r, keep.Mb)


    dt.plot.LD <- bind_rows (lapply (index.Mb, function (filt.Mb) {

      ###### verificar esto con sebas
      #filt.Mb = 0.1


      x.ini.Mb <- x1.na %>%
        dplyr::filter (diff.dist <= 0.1)  ### aca me quedo con la primer 0.1 Mb siempre

      QQ <- quantile (x.ini.Mb$R2)

      l1 <- l1
      l2 <- l2
      l3 <- l3
      q1 <- quantile (QQ, l1)
      q2 <- quantile (QQ, l2)
      q3 <- quantile (QQ, l3)



      ###################
      XX <- data.frame( HUMk=q1 [[1]] ,  H1= q2 [[1]], HLMk =q3 [[1]], chrom = filt.chr)
      ###################


      #print (str_c (filt.Mb, "Mb_chr_", filt.chr))


      x2 <- x1.na %>%
        dplyr::filter (diff.dist <= filt.Mb) #este tiene que ser el argumento que cambiar
      # con algun criterio segun la logitud del cromosoma
      #dplyr::mutate (Mb = "0.1")

      num.total <- nrow (x2) # este hay que usar despues.


      # pp <- # el summary de x2$R2

      p1 <- x2 %>%
        dplyr::filter (R2 > q3[[1]]) %>%
        dplyr::mutate (prob="HLMk") #%>% ### este tiene que cambia

      np1 <- nrow (p1)/num.total

      p2 <- x2 %>%
        dplyr::filter (R2 > q2[[1]]) %>% ### este tiene que cambia
        dplyr::filter (R2 <= q3[[1]]) %>%
        dplyr::mutate (prob= "LMk")### este tiene que cambia

      np2 <- nrow (p2)/num.total

      p3 <- x2 %>%
        dplyr::filter (R2 >  q1[[1]]) %>% ### este tiene que cambia
        dplyr::filter (R2 <= q2[[1]]) %>%
        dplyr::mutate (prob= "UMk") ### este tiene que cambia


      np3 <- nrow (p3)/num.total

      p4 <- x2 %>%
        dplyr::filter (R2 <= q1[[1]]) %>%
        dplyr::mutate (prob= "HUMk") ### este tiene que cambia

      np4 <- nrow (p4)/num.total

      #XX <- rbind (p1, p2,p3,p4) ## este para que esta??

      XX.1 <- as.numeric (rbind (np1, np2, np3, np4))

      XX.1.1 <- c("HLMk","LMk","UMk","HUMk")

      XX.2 <- data.frame (ratio = XX.1, class = XX.1.1, Mb=filt.Mb)

      XX.3 <- data.frame (XX.2, num.total = num.total, chrom= filt.chr)

      #XX.2$clase <- factor(df2$Genotype, levels = c("Genotype 2", "Genotype 3", "Genotype 1")).

      # Convert the cyl variable to a factor
      XX.3$class <- as.factor(XX.3$class)
      XX.3$class <- factor(XX.3$class, levels = c("HUMk", "UMk", "LMk", "HLMk"))

      return (XX.3)

    }))

    df.plot.LD <- as_tibble (dt.plot.LD)

    df.plot.LD.W <- df.plot.LD %>%
      dplyr::select (-num.total)

    df.plot.LD.W <- df.plot.LD %>%
      dplyr::select (-num.total)%>%
      pivot_wider(names_from = class, values_from =ratio) %>%
      dplyr::select (chrom, Mb, everything())

    write_delim (df.plot.LD.W, file=str_c("./Data/procdata/df.freq.LD.Wider_", filt.chr,".txt"),
                 delim = ",", na = "NA")

    df.plot.LD <- df.plot.LD %>%
      dplyr::mutate (num.cat = num.total * ratio) %>%
      dplyr::select (ratio, class, Mb, num.cat,num.total, chrom )

    write_csv (df.plot.LD, file= str_c("./Data/procdata/df.freq.LD.Longer_", filt.chr,".csv"),
               na = "NA", append = FALSE)

    df.plot.1 <- df.plot.LD %>%
      dplyr::select (c(Mb, ratio, class))

    bplt <- ggbarplot (df.plot.1 , "Mb",  "ratio",
                       title = str_c("LD.decay by interval chr_", filt.chr),
                       fill = "class",
                       border ="white",
                       #sort.val = "desc",
                       x.text.angle = 90,
                       xlab = str_c ("(bin ", keep.Mb ,"Mb)"),
                       color = "class",
                       palette =c("navyblue","royalblue3","orange", "red4"),
                       label = FALSE, lab.col = "white", lab.pos = "in")

    bplt1 <- bplt +
      rremove ("x.text")

    bplt1 %>%
      ggexport(filename = str_c("./Figures/plot.freq.decay.",filt.chr,".png"))

    print (bplt1)

    #### df para el el grafico de numero de parwise
    clases <- c("HUMk", "HLMk")

    df.num.HX <-  df.plot.LD %>%
      dplyr::filter (class %in% clases) %>%
      dplyr::select (Mb, class, num.cat, chrom )

    df.num.HX$class <- as.character(df.num.HX$class)

    df.num.total <-  df.plot.LD %>%
      #dplyr::filter (class == "HUMk") %>%
      dplyr::select (Mb, class, num.total, chrom )%>%
      dplyr::rename (num.cat = num.total) %>%
      dplyr::distinct_at (vars(Mb,num.cat, chrom)) %>%
      dplyr::mutate (class =  "Ntotal")%>%
      dplyr::select (Mb, class, num.cat, chrom)

    df.num <-  bind_rows (df.num.HX, df.num.total) %>%
      dplyr::arrange (Mb)

    df.num$class <- as.factor (df.num$class)

    scatter.num <- ggscatter (df.num, x = "Mb", y = "num.cat",
                              color="class",
                              palette = c ("red", "navyblue", "gray48"),
                              title = str_c("num.parwise by interval chr_", filt.chr),
                              xlab = str_c ("(bin ", keep.Mb ,"Mb)"),
                              ylab = "num.parwise")

    scatter.num %>%
      ggexport(filename = str_c("./Figures/plot.scatter.num.chr_",filt.chr,".png"))

    df.num.HLMk <- df.num %>%
      dplyr::filter (class == "HLMk")

    ##### esta figura no se si vale la pena
    scatter.num.HLMk <- ggscatter (df.num.HLMk, x = "Mb", y = "num.cat",
                                   color="red",
                                   title = str_c("num.parwise by interval_HLMk_chr_", filt.chr),
                                   xlab = str_c ("(bin ", keep.Mb ,"Mb)"),
                                   ylab = "num.parwise")

    scatter.num.HLMk %>%
      ggexport(filename = str_c("./Figures/plot.scatter.HLMk_chr_",filt.chr,".png"))

    print (scatter.num)
    print (scatter.num.HLMk)

    ## datos para el grafico de proporciones

    df.plot.LD.HUMk <- df.plot.LD %>%
      dplyr::filter (class == "HUMk")

    df.plot.LD.HLMk <- df.plot.LD %>%
      dplyr::filter (class == "HLMk")

    ### revisar el ratio si esta bien o hay otra mejor forma

    df.plot.prob.HLMk_HUMk <-  df.plot.LD.HLMk %>%
      dplyr::mutate (ratio.HLHU = num.cat/df.plot.LD.HUMk$num.cat) %>%
      dplyr::select (ratio.HLHU, Mb,  chrom)

    df.plot.prob.HLMk_HUMk$ratio.HLHU [which(!is.finite( df.plot.prob.HLMk_HUMk$ratio.HLHU))] <- NA

    df.plot.LD.HLMk_HUMk <- bind_rows (df.plot.LD.HLMk, df.plot.LD.HUMk) %>%
      dplyr::arrange (Mb)

    df.plot.LD.HLMk_HUMk <- df.plot.LD.HLMk_HUMk %>%
      dplyr::select ( ratio, class, Mb,chrom )


    max.pHLMk <- round (max (df.plot.LD.W$HLMk ), 2)

    print (str_c("La maxima probalidad para HLMk es del ", max.pHLMk*100, "%" ))

    if (prob.HLMK/100 > max.pHLMk) {

      print ( str_c("La probalidad para HLMk usada fue del ", max.pHLMk*100, "%" ))

      umbral.H <- df.plot.LD.W %>%
        dplyr::filter (HUMk < prob.HUMK/100 & HLMk  >= max.pHLMk ) %>%
        dplyr::filter (Mb == max(Mb))

    }

    if ( prob.HLMK/100 < max.pHLMk ) {

      print (str_c("La probalidad para HLMk usada fue del ", prob.HLMK, "%" ))

      umbral.H <- df.plot.LD.W %>%
        dplyr::filter (HUMk < prob.HUMK/100 & HLMk  >= prob.HLMK/100) %>%
        dplyr::filter (Mb == max(Mb))

    }

    scatter.prop.total <- ggscatter (df.plot.LD.HLMk_HUMk,  x = "Mb", y = "ratio",
                                     ylim=c(0,1),
                                     title = str_c("proption by interval_chr_", filt.chr),
                                     color = "class",
                                     palette = c( "navyblue",  "red"),
                                     xlab = str_c ("(bin ", keep.Mb ,"Mb)")) +
      geom_hline(yintercept = umbral.H$HUMk, color = "navyblue") +
      geom_hline(yintercept = umbral.H$HLMk, color ="red") +
      geom_vline(xintercept = umbral.H$Mb, color ="black")

    scatter.prop.total %>%
      ggexport(filename = str_c("./Figures/scatter.prop.total.chr_",filt.chr,".png"))

    print (scatter.prop.total)


    umbral.HH <- df.plot.prob.HLMk_HUMk %>%
      dplyr::filter (Mb == umbral.H$Mb )

    scatter.propHLMk_HUMk <- ggscatter (df.plot.prob.HLMk_HUMk, x = "Mb", y = "ratio.HLHU",
                                        title = str_c("proportion HLMk/HUMk_chr_", filt.chr),
                                        #ylim=c(0,1),
                                        color = "black",
                                        xlab =  str_c ("(bin ", keep.Mb ,"Mb)")) +
      geom_vline(xintercept = umbral.H$Mb, color ="black") +
      geom_hline(yintercept = umbral.HH$ratio.HLHU, color ="black")


    scatter.propHLMk_HUMk %>%
      ggexport(filename = str_c("./Figures/scatter.propHLMk_HUMk_chr",filt.chr,".png"))

    print (scatter.propHLMk_HUMk)


    umbral.H <- umbral.H %>%
      dplyr::mutate (HLMk_HUMk = umbral.HH$ratio.HLHU)

    return (umbral.H)

  }))

  dt.all.chrom <- dt.all.chrom  %>%
    dplyr::mutate ( Mb = round   ( Mb, 2))   %>%
    dplyr::mutate ( HLMk = round ( HLMk, 2)) %>%
    dplyr::mutate ( LMk = round ( LMk, 2)) %>%
    dplyr::mutate ( UMk = round ( UMk, 2)) %>%
    dplyr::mutate ( HUMk = round ( HUMk, 2)) %>%
    dplyr::mutate ( HLMk_HUMk = round ( HLMk_HUMk, 2))


  write_csv (dt.all.chrom , file= str_c("./Data/procdata/umbral.HLMk_HUMk_", id.cross, ".csv"),
             na = "NA", append = FALSE)

}



## cargar los paquetes
library("Matrix")
library("RColorBrewer")
library(lme4)
library(nlme)
library(car)
library("ggplot2")
library("lattice")
library("latticeExtra")
library(multcompView)
library(dplyr)
library(xtable)
library(tidyverse)
library (emmeans)
library("qtl")
library(stringr)
library(data.table)
library(svMisc)
library(ggpubr)
library("ggsci")
library("FactoMineR")
library("factoextra")
library("corrplot")
library(gameofthrones)
library(LDheatmap)
library(genetics)
library("assertive")

## Cargo las matrices de datos ############

#start_time <- Sys.time()

# Time difference of 1.000327 min
## esto lo tenes setear
index.chr <- 1:4

# se cargan las matrices generadas antes
# se genear una lista con 20 tibbles
GENO.crom <- lapply (index.chr, function (filtro) {

  G <- read_delim (file=paste("./Data/EEMAC.cross.1LD.decay.LG",filtro,".txt",
                              sep=""), delim = ",",
                   na = "NA", quote = "\"",col_names = TRUE)


  G <- G %>%
       dplyr::mutate (chrom = filtro)

  #G <- G %>%
  #    dplyr::mutate (chr = str_c("chr_",filtro))
  print (G)
  return (G)
})

df.geno.crom <- do.call (rbind, GENO.crom)

###############################
#### funciones para analisis de LD
max (df.geno.crom$diff.dist)


#### pruebo solo el LG.1

  df.geno.crom.2 <- read_delim ("./Data/EEMAC.cross.1LD.decay.LG2.txt", delim = ",",
                   na = "NA", quote = "\"",col_names = TRUE)


  summary ( df.geno.crom.2 )


  df.geno.crom.1 <- df.geno.crom.1 %>%
                    dplyr::mutate (chrom = 1)

  max (df.geno.crom.1$diff.dist)

#########
# Argumentos anteriores

id.cross = "EEMAC.cross.1"
data= df.geno.crom
ini = 1
l1= 0.25
l2=0.5
l3=0.75
seq1 = c(0.5, 1,10, 50)
distance.unit = "cM"

run_table_quantiles <- function (id.cross=NULL,data= NULL, ini = 1, l1= 0.25, l2=0.5, l3=0.75, seq1= NULL, distance.unit = NULL) {

  if   ( distance.unit != "cM" &  distance.unit != "Mb" & distance.unit != "Kb" & distance.unit != "bp" ) {

    stop (str_c ("Debe definir una unidad de distancia valida"))


  }

  if (is.null (distance.unit)) {

    stop (str_c ("Falta definir unidad de distancia"))

  }

  if   ( distance.unit == "Mb" | distance.unit == "Kb" | distance.unit == "bp") {

    print (str_c ("El mapa es un mapa fisico"))

  }

  if   ( distance.unit == "cM") {
    print (str_c ("El mapa es un mapa genetico"))

  }


  list.chrom <- unique (data$chrom)
  #filt.crom =2  ###

  df.qq <- bind_rows (lapply (list.chrom, function (filt.crom){

    print(filt.crom)

    dat.1 <- data %>%
      dplyr::filter (chrom == filt.crom) %>%
      dplyr:::filter (R2 != "NA")

    list.dist <- seq1

    # filt.dist = 1 ######!!!!!!!

    df.hist.plot <- bind_rows (lapply (list.dist, function (filt.dist){

      print (filt.dist)

      ############### hay que revisar estos numeros ##############
      # if   ( distance.unit == "Mb" ) {

      x.dist.Mb <- dat.1 %>%
        dplyr::filter (diff.dist <= filt.dist * 1e5)

      bin <- (filt.dist* 1e5)/1e6

      x.dist.Mb.1 <- x.dist.Mb %>%
        dplyr::mutate (bin = str_c (bin, "Mb"))

      #print (x.dist.Mb.1)

      #}

      #########################################################################

      ############### hay que revisar estos numeros ##############
      #if   ( distance.unit == "cM" ) {

      # x.dist.cM <- dat.1 %>%
      #  dplyr::filter ( diff.dist <= filt.dist)

      #  bin <- filt.dist

      # x.dist.cM <-x.dist.cM %>%
      #            dplyr::mutate (bin = str_c (bin, "cM"))
      #}

      ############### hay que revisar estos numeros ##############

    }))

    # ggv <-  ggviolin(df.violin.plot, x = "bin", y = "R2", title = filt.crom,
    #          color = "black", fill = "azure")

    #print (ggv)


    ggh <- gghistogram (df.hist.plot, y = "..density..", x = "R2", facet.by = "bin",title =filt.crom,
                        #fill = "lightgray",
                        fill = "bin", palette =   "RdBu",
                        add = "median", rug = TRUE)
    print (ggh)

    ggh  %>%
      ggexport(filename = str_c("./Figures/ggh_",id.cross, "_", filt.crom,".png"))


    qqunif <- ggplot(df.hist.plot, aes(R2, colour =  bin , fill = bin), size=1.2) + stat_ecdf() +
      theme_bw()+
      stat_function(fun=punif,args=list(0,1))+
      scale_color_manual(values=c("green", "black", "blue", "red"))+
      labs(title=str_c("ECDF and theoretical CDF","_", filt.crom)) +
      labs(y = "Theoretical Quantiles", x = "Sample Quantiles")

    print (qqunif)

    qqunif %>%
      ggexport(filename = str_c("./Figures/qqunif_",id.cross, "_", filt.crom,".png"))

    # if   ( distance.unit == "cM" ) {

    # Nota: aca los estoy estoy tomando cada 0.5 cM hay que ver si esto reemplaza a 0.1 Mb

    #list.seqCM <- seq (from = 0.5, to = 20, by = 0.5)
    #filt.distCM <- 0.5
    #df.QQ.seq.cM <- bind_rows (lapply (list.seqCM, function (filt.distCM){

    # x.ini.cM <- dat.1 %>%
    #  dplyr::filter (diff.dist <= ini * filt.distCM)

    #QQ <- quantile (x.ini.cM$R2)

    #l1 <- l1
    #l2 <- l2
    #l3 <- l3
    #  q1 <- quantile (QQ, l1)
    # q2 <- quantile (QQ, l2)
    #  q3 <- quantile (QQ, l3)

    # XX <- data.frame( HLE=round (q1 [[1]],2) ,  H1= round (q2 [[1]],2), HLD = round (q3 [[1]],2), chrom = filt.crom)

    #dt.QQ.seq.cM <- XX %>%
    # dplyr::mutate (inter.cM = filt.distCM ) %>%
    #dplyr::select (chrom, inter.cM, HLE,   H1,  HLD )

    #}))


    #df.QQ.seq.cM.long <- df.QQ.seq.cM %>%
    # pivot_longer( ! c(chrom,inter.cM), names_to = "cat.LD", values_to = "unqq.R2")


    #qq.scatt <- ggscatter (df.QQ.seq.cM.long, x = "inter.cM", y = "unqq.R2",
    #                      title = str_c("Caida de qqR2 en funcion del bin", id.cross, "_", filt.crom),
    #                     color = "cat.LD", shape = "cat.LD",
    #                    palette = c("navyblue", "gray48", "darkorange"),
    #                   add = "loess", rug= TRUE)

    #print (qq.scatt)

    #qq.scatt  %>%
    # ggexport(filename = str_c("./Figures/qq.scatt_",id.cross, "_", filt.crom,".png"))

    #  }

  }))
  return (df.qq)
}




quantiles.EEMAC <- run_table_quantiles (id.cross = "EEMAC.cross.1",
                                        data= df.geno.crom,
                                        ini = 1,
                                        l1= 0.25,
                                        l2=0.5,
                                        l3=0.75,
                                        seq1 = c(0.5, 1,10, 50),
                                        distance.unit = "cM")


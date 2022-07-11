## Evan Newell 2015 Edited by Tan Yong Kee

#' tSNE and OneSENSE algorithm for FCS data
#'
#' @param LoaderPATH Path where FCS file is located
#' @param ceil Maximum number of cells to sample from each fcs sample/file
#' @param FNnames .csv file generated when markers from each category
#'     are selected
#' @param OutputSuffix suffix to name output folder
#' @param DotSNE boolean, if TRUE do tSNE, if FALSE skip tSNE
#' @param DoOneSENSE boolean, if TRUE do OneSENSE, if FALSE skip OneSENSE
#' @param Bins number of bins to put the cell data into, DEFAULT = 250
#'
#' @return FCS files, tSNE histograms, OneSENSE plot
#'
#' @importFrom Rtsne Rtsne
#' @importFrom graphics hist
#' @importFrom flowCore read.flowSet exprs keyword write.FCS logicleTransform
#'     inverseLogicleTransform identifier exprs<- identifier<-
#' @importFrom methods cbind2
#'
#' @examples
#' #dir <- system.file('extdata/fcs', package='oneSENSE')
#' #FCStSNE(LoaderPATH=dir, FNnames=fnnames) #remove hash symbol to run
FCStSNE <- function(LoaderPATH = "fcs",
                    ceil = 5000,
                    FNnames = "names.csv",
                    local = NULL,
                    OutputSuffix = "Out",
                    DotSNE = TRUE,
                    DoOneSENSE = TRUE,
                    multiopt = F,
                    Bins = 250) {
  fs <- read.flowSet(path = LoaderPATH, pattern = ".fcs$")
  FcsFileNames <- rownames(keyword(fs, "FILENAME"))
  NumBC <- length(fs)  #3
  FFdata <- NULL  #FlowFrame data

  for (FFs in 1:NumBC) {
    # iterate through each FCS file
    FFt <- exprs(fs[[FFs]])
    ## Downsample ##
    if (nrow(FFt) <= ceil)
      FFa <- FFt
    else
      FFa <- FFt[sample(nrow(FFt), ceil,
                        replace = FALSE),]
    colnames(FFa) <- fs[[FFs]]@parameters$desc
    empties <-
      which(is.na(colnames(FFa)) | colnames(FFa) == " ")
    colnames(FFa)[empties] <-
      fs[[FFs]]@parameters$name[empties]
    fs[[FFs]]@parameters$desc <- colnames(FFa)
    FFa <- cbind(FFa, rep(FFs, dim(FFa)[1]))
    colnames(FFa)[dim(FFa)[2]] <- "InFile"
    # Concatenate
    FFdata <- rbind(FFdata, FFa)
  }
  message("FCS Files Read")
  if (is.null(local)) {
    keeptable <-
      read.csv(
        paste(
          "/media/data/html/source/CIPHEoneSENSE/OUTPUT",
          "names.csv",
          sep = .Platform$file.sep
        )
      )
  }
  else {
    keeptable <- read.csv(paste(local,
                                "names.csv",
                                sep = .Platform$file.sep))
  }
  keeprowbool <- sapply(keeptable[, 2], function(x)
    any(x == "Y"))
  keeprows <- subset(keeptable, keeprowbool)
  data <- FFdata[, which(colnames(FFdata) %in% keeprows[, 1])]
  data <- as.data.frame(data)
  lgcl <- logicleTransform(w = 0.25,
                           t = 16409,
                           m = 4.5,
                           a = 0)
  ilgcl <- inverseLogicleTransform(trans = lgcl)
  data1 <- apply(data, 2, lgcl)
  FFdata1 <- apply(FFdata, 2, lgcl)
  score <- NULL
  tSNEmat <- NULL
  if (DotSNE) {
    message("Doing tSNE")
    if (!multiopt) {
      tSNEdata3 <- Rtsne(data1, dims = 2, check_duplicates = FALSE)
      tSNEmat <- tSNEdata3$Y
      colnames(tSNEmat) <- c("tSNE1", "tSNE2")
    }
    else{
      ###Added functionality for multiopt analysis
      #file_to_read<-read.FCS("/media/data/html/source/CIPHEoneSENSE/raw_mass/c01_21-Aug-04_EXT186Naquet-tube1_MG_02_1_Beads_2-.fcs")@exprs
      roots.temp <-
        "/media/data/html/source/CIPHEoneSENSE/raw_mass/"
      rep <- 'test'
      # analyzed_markers <- c("Dy163Di", "Gd156Di", "Ho165Di")
      write.csv(
        data1,
        paste0(roots.temp, rep, ".csv"),
        quote = FALSE,
        row.names = FALSE
      )
      a <- paste0(roots.temp, rep, ".csv")
      b <- paste0(roots.temp, rep, "_tsne.csv")
      temp <- paste0(roots.temp, rep, ".log")
      #iter <- 1000
      iter <- 100
      run <-
        "/media/data/cyto/MultiOptTSNE/Multicore-opt-SNE/MulticoreTSNE/run/run_optsne.py"
      print('running python')
      cmd <- paste0(
        "python2 ",
        run,
        " --optsne --data ",
        a
        ,
        " --outfile ",
        b,
        " --n_components 2 --n_threads 35 --perp 30 --early_exaggeration 12 --n_iter ",
        iter,
        " ",
        "> ",
        temp,
        " 2>&1"
      )
      system(cmd)
      #roots.temp <- "/media/data/cyto/MultiOptTSNE/TEMP/"
      c <- paste0(roots.temp, rep, "_tsne.csv")
      mat <- read.csv(c, header = FALSE)
      # tester 2-3 essais
      colnames(mat) <- c('TSNE1', 'TSNE2')
      tSNEmat <- mat
    }

    plot(
      tSNEmat[, 1],
      tSNEmat[, 2],
      pch = ".",
      xlab = "tSNE1",
      ylab = "tSNE2",
      cex = 0.1
    )
  } else
    tSNEmat <- NULL


  if (DoOneSENSE) {
    message("Doing oneSENSE")
    Xx1DtSNEmat <- NULL
    if (dim(keeptable)[2] == 4) {
      for (factor in 2:(dim(keeptable)[2])) {
        # for loop from 2 to 4
        OneDtSNEname <- colnames(keeptable)[factor]
        keeprowbool <- sapply(keeptable[, factor],
                              function(x)
                                any(x == "Y"))
        keeprows <- subset(keeptable, keeprowbool)
        dataX <-
          FFdata1[, which(colnames(FFdata1) %in% keeprows[, 1])]
        dataX <- as.matrix(dataX)
        print('datax')
        print(head(dataX))
        if (!multiopt) {
          tSNEdata3 <- Rtsne(dataX,
                             dims = 1,
                             check_duplicates = FALSE)
          tSNEmat1 <- tSNEdata3$Y
          colnames(tSNEmat1) <- OneDtSNEname
        }
        else{
          ###Added functionality for multiopt analysis
          #file_to_read<-read.FCS("/media/data/html/source/CIPHEoneSENSE/raw_mass/c01_21-Aug-04_EXT186Naquet-tube1_MG_02_1_Beads_2-.fcs")@exprs
          roots.temp <-
            "/media/data/html/source/CIPHEoneSENSE/raw_mass/"
          rep <- paste0(factor, '_datax')
          # analyzed_markers <- c("Dy163Di", "Gd156Di", "Ho165Di")
          write.csv(
            dataX,
            paste0(roots.temp, rep, ".csv"),
            quote = FALSE,
            row.names = FALSE
          )
          print('csv written')
          a <- paste0(roots.temp, rep, ".csv")
          b <- paste0(roots.temp, rep, "_tsne.csv")
          temp <- paste0(roots.temp, rep, ".log")
          #iter <- 1000
          iter <- 100
          run <-
            "/media/data/cyto/MultiOptTSNE/Multicore-opt-SNE/MulticoreTSNE/run/run_optsne.py"
          print('running python')
          cmd <- paste0(
            "python2 ",
            run,
            " --optsne --data ",
            a
            ,
            " --outfile ",
            b,
            " --n_components 1 --n_threads 35 --perp 30 --early_exaggeration 12 --n_iter ",
            iter,
            " ",
            "> ",
            temp,
            " 2>&1"
          )
          system(cmd)
          #roots.temp <- "/media/data/cyto/MultiOptTSNE/TEMP/"
          c <- paste0(roots.temp, rep, "_tsne.csv")
          mat <- read.csv(c, header = FALSE)
          # tester 2-3 essais
          colnames(mat) <- c('TSNE1', 'TSNE2')
          mat <- as.matrix(mat)
          tSNEmat1 <- mat
        }
        #     tSNEmat1 <- tSNEdata3$Y
        #    colnames(tSNEmat1) <- OneDtSNEname
        hist(tSNEmat1, 100)
        Xx1DtSNEmat <- cbind(Xx1DtSNEmat, tSNEmat1)
      }
      scatterplot3d(
        x = Xx1DtSNEmat[, 1],
        y = Xx1DtSNEmat[, 2],
        z = Xx1DtSNEmat[, 3],
        pch = ".",
        xlab = paste("tSNE1 ", colnames(Xx1DtSNEmat)[1],
                     sep = ""),
        ylab = paste("tSNE2 ", colnames(Xx1DtSNEmat)[2],
                     sep = ""),
        zlab = paste("tSNE3 ", colnames(Xx1DtSNEmat)[3],
                     sep = "")
      )
    } else {
      # loop from 2 to 4
      for (factor in 2:(dim(keeptable)[2])) {
        print('factor')
        print(factor)
        OneDtSNEname <- colnames(keeptable)[factor]
        keeprowbool <- sapply(keeptable[, factor],
                              function(x)
                                any(x == "Y"))
        keeprows <- subset(keeptable, keeprowbool)
        dataX <-
          FFdata1[, which(colnames(FFdata1) %in% keeprows[, 1])]
        print('datax')
        print(head(dataX))
        dataX <- as.matrix(dataX)
        if(!multiopt){
        tSNEdata3 <-
          Rtsne(dataX, dims = 1, check_duplicates = FALSE)
        tSNEmat1 <- tSNEdata3$Y
        colnames(tSNEmat1) <- OneDtSNEname}
        else{
          ###Added functionality for multiopt analysis
          #file_to_read<-read.FCS("/media/data/html/source/CIPHEoneSENSE/raw_mass/c01_21-Aug-04_EXT186Naquet-tube1_MG_02_1_Beads_2-.fcs")@exprs
          roots.temp <-
            "/media/data/html/source/CIPHEoneSENSE/raw_mass/"
          rep <- paste0(factor, '_datax')
          # analyzed_markers <- c("Dy163Di", "Gd156Di", "Ho165Di")
          write.csv(
            dataX,
            paste0(roots.temp, rep, ".csv"),
            quote = FALSE,
            row.names = FALSE
          )
          print('csv written')
          a <- paste0(roots.temp, rep, ".csv")
          b <- paste0(roots.temp, rep, "_tsne.csv")
          temp <- paste0(roots.temp, rep, ".log")
          #iter <- 1000
          iter <- 100
          run <-
            "/media/data/cyto/MultiOptTSNE/Multicore-opt-SNE/MulticoreTSNE/run/run_optsne.py"
          print('running python')
          cmd <- paste0(
            "python2 ",
            run,
            " --optsne --data ",
            a
            ,
            " --outfile ",
            b,
            " --n_components 1 --n_threads 35 --perp 30 --early_exaggeration 12 --n_iter ",
            iter,
            " ",
            "> ",
            temp,
            " 2>&1"
          )
          system(cmd)
          print('system cmd done')
          #roots.temp <- "/media/data/cyto/MultiOptTSNE/TEMP/"
          c <- paste0(roots.temp, rep, "_tsne.csv")
          mat <- read.csv(c, header = FALSE)
          # tester 2-3 essais
          colnames(mat) <- OneDtSNEname
          print('mat')
          mat <- as.matrix(mat)
          print('head mat')
          tSNEmat1 <- mat
        }
        print('head tsnemat1')
        print(tSNEmat1)
        hist(tSNEmat1, 100)
        Xx1DtSNEmat <- cbind(Xx1DtSNEmat, tSNEmat1)
        print('Xx1DtSNEmat')
        print(head(Xx1DtSNEmat))
      }
      plot(
        Xx1DtSNEmat[, 1],
        Xx1DtSNEmat[, 2],
        pch = ".",
        xlab = paste("tSNE1 ", colnames(Xx1DtSNEmat)[1], sep = ""),
        ylab = paste("tSNE2 ",
                     colnames(Xx1DtSNEmat)[2], sep = ""),
        cex = 1
      )
    }
  } else
    Xx1DtSNEmat <- NULL
  print('b4 nxx1')
  NXx1 <- apply(Xx1DtSNEmat, 2,
                function(x)
                  ((x - min(x)) / (max(x) - min(x))) * 10000)
  tSNEmat <- as.matrix(tSNEmat)
  colnames(tSNEmat) <- c('TSNE1', 'TSNE2')
  score2 <- cbind(score, tSNEmat)
  Nscore <- apply(score2, 2,
                  function(x)
                    ((x - min(x)) / (max(x) - min(x))) * 3.7)
  ilgcl <- inverseLogicleTransform(trans = lgcl)
  NIscore <- apply(Nscore, 2, ilgcl)
  colnames(NIscore) <- colnames(score2)
  NIscore <- cbind(NIscore, NXx1)
  message("Writing FCS Files")
  # output new FCS files
  for (FFs in 1:NumBC) {
    newFF <- fs[[1]]
    newBaseData <-
      FFdata[FFdata[, dim(FFdata)[2]] == FFs,-dim(FFdata)[2]]
    colnames(newBaseData) <- colnames(exprs(newFF))
    exprs(newFF) <- newBaseData
    subsetNIscore <- NIscore[FFdata[, dim(FFdata)[2]] == FFs,]
    colnames(subsetNIscore)[colnames(subsetNIscore) == 'input'] <-
      'xOneSense'
    colnames(subsetNIscore)[colnames(subsetNIscore) == 'input2'] <-
      'yOneSense'
    colnames(subsetNIscore)[colnames(subsetNIscore) == 'input3'] <-
      'zOneSense'
    # newFF <- cbind2(newFF, subsetNIscore)

    newFF <- flowCore::fr_append_cols(newFF, subsetNIscore)
    #  newFF@parameters$desc <- colnames(cbind(newBaseData, subsetNIscore))
    newFF@parameters$desc <-
      c(fs[[FFs]]@parameters$desc, colnames(subsetNIscore))
    suppressWarnings(dir.create(paste0(LoaderPATH, "/", "Output")))
    BaseFN <-
      sapply(strsplit(FcsFileNames[FFs], split = "\\."), "[", 1)
    # FNresult <- paste0(LoaderPATH, "_", OutputSuffix,
    #                 "/", BaseFN, "_", OutputSuffix, ".fcs")


    if (is.null(local)) {
      FNresult <- paste0(LoaderPATH,
                         "/",
                         "Output",
                         "/",
                         BaseFN,
                         "_",
                         OutputSuffix,
                         ".fcs")
    }
    else {
      FNresult <- paste0(local,
                         "/", BaseFN, "_", OutputSuffix, ".fcs")
    }


    newFF@description$"$FIL" <-
      paste0(BaseFN, "_", OutputSuffix, ".fcs")
    newFF@description$FILENAME <-
      paste0(BaseFN, "_", OutputSuffix, ".fcs")
    identifier(newFF) <- paste0(BaseFN, "_", OutputSuffix)
    write.FCS(newFF, FNresult)
  }
  message("New FCS Files Written")
}

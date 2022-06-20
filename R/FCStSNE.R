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
            FFa <- FFt else FFa <- FFt[sample(nrow(FFt), ceil,
                                                replace = FALSE), ]
            colnames(FFa) <- fs[[FFs]]@parameters$desc
            empties <- which(is.na(colnames(FFa)) | colnames(FFa) == " ")
            colnames(FFa)[empties] <- fs[[FFs]]@parameters$name[empties]
            fs[[FFs]]@parameters$desc <- colnames(FFa)
            FFa <- cbind(FFa, rep(FFs, dim(FFa)[1]))
            colnames(FFa)[dim(FFa)[2]] <- "InFile"
            # Concatenate
            FFdata <- rbind(FFdata, FFa)
    }
    message("FCS Files Read")
    if (is.null(local)){
    keeptable <- read.csv(paste(LoaderPATH,
                                "Output/names.csv",
                                sep = .Platform$file.sep))}
    else {
      keeptable <- read.csv(paste(local,
                                  "names.csv",
                                  sep = .Platform$file.sep))
    }
    keeprowbool <- sapply(keeptable[, 2], function(x) any(x == "Y"))
    keeprows <- subset(keeptable, keeprowbool)
    data <- FFdata[, which(colnames(FFdata) %in% keeprows[, 1])]
    data <- as.data.frame(data)
    lgcl <- logicleTransform(w = 0.25, t = 16409, m = 4.5, a = 0)
    ilgcl <- inverseLogicleTransform(trans = lgcl)
    data1 <- apply(data, 2, lgcl)
    FFdata1 <- apply(FFdata, 2, lgcl)
    score <- NULL
    tSNEmat <- NULL
    if (DotSNE) {
        message("Doing tSNE")
        tSNEdata3 <- Rtsne(data1, dims = 2,check_duplicates = FALSE)
        tSNEmat <- tSNEdata3$Y
        colnames(tSNEmat) <- c("tSNE1", "tSNE2")
        plot(tSNEmat[, 1], tSNEmat[, 2], pch = ".",
            xlab = "tSNE1", ylab = "tSNE2", cex = 0.1)
    } else tSNEmat <- NULL


    if (DoOneSENSE) {
        message("Doing oneSENSE")
        Xx1DtSNEmat <- NULL
        if (dim(keeptable)[2] == 4) {
        for (factor in 2:(dim(keeptable)[2])) {
            # for loop from 2 to 4
            OneDtSNEname <- colnames(keeptable)[factor]
            keeprowbool <- sapply(keeptable[, factor],
                                function(x) any(x == "Y"))
            keeprows <- subset(keeptable, keeprowbool)
            dataX <- FFdata1[, which(colnames(FFdata1) %in% keeprows[, 1])]
            dataX <- as.matrix(dataX)
            tSNEdata3 <- Rtsne(dataX, dims = 1, check_duplicates = FALSE)
            tSNEmat1 <- tSNEdata3$Y
            colnames(tSNEmat1) <- OneDtSNEname
            hist(tSNEmat1, 100)
            Xx1DtSNEmat <- cbind(Xx1DtSNEmat, tSNEmat1)
            }
        scatterplot3d(x = Xx1DtSNEmat[, 1], y = Xx1DtSNEmat[, 2],
                        z = Xx1DtSNEmat[, 3], pch = ".",
                        xlab = paste("tSNE1 ", colnames(Xx1DtSNEmat)[1],
                                    sep = ""),
                        ylab = paste("tSNE2 ", colnames(Xx1DtSNEmat)[2],
                                    sep = ""),
                        zlab = paste("tSNE3 ", colnames(Xx1DtSNEmat)[3],
                                    sep = ""))
    } else {
        # loop from 2 to 4
        for (factor in 2:(dim(keeptable)[2])) {
        OneDtSNEname <- colnames(keeptable)[factor]
        keeprowbool <- sapply(keeptable[, factor],
                                function(x) any(x == "Y"))
        keeprows <- subset(keeptable, keeprowbool)
        dataX <- FFdata1[, which(colnames(FFdata1) %in% keeprows[, 1])]
        dataX <- as.matrix(dataX)
        tSNEdata3 <- Rtsne(dataX, dims = 1, check_duplicates = FALSE)
        tSNEmat1 <- tSNEdata3$Y
        colnames(tSNEmat1) <- OneDtSNEname
        hist(tSNEmat1, 100)
        Xx1DtSNEmat <- cbind(Xx1DtSNEmat, tSNEmat1)
    }
    plot(Xx1DtSNEmat[, 1], Xx1DtSNEmat[, 2], pch = ".",
        xlab = paste("tSNE1 ", colnames(Xx1DtSNEmat)[1], sep = ""),
                    ylab = paste("tSNE2 ",
                                colnames(Xx1DtSNEmat)[2], sep = ""), cex = 1)
        }
    } else Xx1DtSNEmat <- NULL
    NXx1 <- apply(Xx1DtSNEmat, 2,
                    function(x) ((x - min(x))/(max(x) - min(x))) * 10000)

    score2 <- cbind(score, tSNEmat)
    Nscore <- apply(score2, 2,
                    function(x) ((x - min(x))/(max(x) - min(x))) * 3.7)
    ilgcl <- inverseLogicleTransform(trans = lgcl)
    NIscore <- apply(Nscore, 2, ilgcl)
    colnames(NIscore) <- colnames(score2)
    NIscore <- cbind(NIscore, NXx1)
    message("Writing FCS Files")
    # output new FCS files
    for (FFs in 1:NumBC) {
        newFF <- fs[[1]]
        newBaseData <- FFdata[FFdata[, dim(FFdata)[2]] == FFs, -dim(FFdata)[2]]
        colnames(newBaseData) <- colnames(exprs(newFF))
        exprs(newFF) <- newBaseData
        subsetNIscore <- NIscore[FFdata[, dim(FFdata)[2]] == FFs, ]
        colnames(subsetNIscore)[colnames(subsetNIscore) == 'input'] <- 'xOneSense'
        colnames(subsetNIscore)[colnames(subsetNIscore) == 'input2'] <- 'yOneSense'
        newFF <- flowCore::fr_append_cols(newFF, subsetNIscore)
        newFF@parameters$desc <- c(fs[[FFs]]@parameters$desc, colnames(subsetNIscore))
        suppressWarnings(dir.create(paste0(LoaderPATH, "/", "Output")))
        BaseFN <- sapply(strsplit(FcsFileNames[FFs], split = "\\."), "[", 1)

        if (is.null(local)){
          FNresult <- paste0(LoaderPATH, "/", "Output",
                             "/", BaseFN, "_", OutputSuffix, ".fcs")}
        else {
          FNresult <- paste0(local,
                             "/", BaseFN, "_", OutputSuffix, ".fcs")
        }


        newFF@description$"$FIL" <- paste0(BaseFN, "_", OutputSuffix, ".fcs")
        newFF@description$FILENAME <- paste0(BaseFN, "_", OutputSuffix, ".fcs")
        identifier(newFF) <- paste0(BaseFN, "_", OutputSuffix)
        write.FCS(newFF, FNresult)
    }
    message("New FCS Files Written")
}

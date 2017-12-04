#  =============================================================================
#  TCGA-Assembler version 2
#
#  Copyright (C) <2017>  <Yitan Zhu>
#  This file is part of TCGA-Assembler.
#
#  TCGA-Assembler is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  TCGA-Assembler is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with TCGA-Assembler.  If not, see <http://www.gnu.org/licenses/>.
#  =============================================================================


#  =============================================================================
#  TCGA-Assembler Version 2 Module B
#  =============================================================================


#  =============================================================================
#  variable prefix example
#  =============================================================================
#  s : string
#  n : number
#  v : vector
#  d : data.frame
#  m : matrix
#  l : list
#  lv : list of vector
#  ld : list of data.frame


#  =============================================================================
#  library
#  =============================================================================

library(RCurl)
library(httr)
library(stringr)
#library(HGNChelper)
#library(rjson)


#  =============================================================================
#  Main Functions of Module B
#  =============================================================================

CombineMultiPlatformData <- function(inputDataList,
                                     combineStyle = "Intersect") {
    options(warn = -1)
    vl <- inputDataList
    for (i in 1:length(vl)) {
        # Keep only one sample for a tissue type of a patient.
        # Usually, there is only one sample of a tissue type of a patient existing in data.
        vl[[i]]$Data <- ToPatientData(vl[[i]]$Data)
        if (i == 1) {
            sampleId <- colnames(vl[[i]]$Data)
        }
        # For each genomic feature, keep only one row of data.
        #lResult <- CombineRedundantFeature(Data = vl[[i]]$Data, Des = vl[[i]]$Des)
        lResult <- list(Data = vl[[i]]$Data,
                        Des = vl[[i]]$Des)
        vl[[i]]$Des <- lResult$Des
        vl[[i]]$Data <- lResult$Data
        if (combineStyle == "Intersect") {
            sampleId <- sort(intersect(sampleId,
                                       colnames(vl[[i]]$Data)))
        }
        if (combineStyle == "Union") {
            sampleId <- sort(union(sampleId,
                                   colnames(vl[[i]]$Data)))
        }
        if (dim(vl[[i]]$Des)[2] == 1) {
            vl[[i]]$Des <- cbind(vl[[i]]$Des,
                                 Description = cbind(rep("",
                                                         dim(vl[[i]]$Des)[1])))
        }
        else {
            if ((dim(vl[[i]]$Des)[2] == 3) && (vl[[i]]$dataType == "copyNumber")) {
                vl[[i]]$Des <- cbind(vl[[i]]$Des[,
                                     1,
                                     drop = FALSE],
                Description = paste(vl[[i]]$Des[, 2],
                                    vl[[i]]$Des[, 3],
                                    sep = ""))
            }
        }
        # combining methylation data at CpG site level.
        if (vl[[i]]$dataType == "methylation") {
            if (sum(toupper(c("REF",
                              "geneSymbolId",
                              "ChromosomeID",
                              "CoordinateID")) %in% toupper(colnames(vl[[i]]$Des))) == 4) {
                vl[[i]]$Des <- cbind(vl[[i]]$Des[, "geneSymbolId", drop = FALSE],
                                     Description = paste(vl[[i]]$Des[, "REF"],
                                                         vl[[i]]$Des[, "ChromosomeID"],
                                                         vl[[i]]$Des[, "CoordinateID"],
                                                         sep = "|"))
                vl[[i]]$Des[, "Description"] <- gsub(" ", "", vl[[i]]$Des[, "Description"])
            }
        }
        vl[[i]]$Des <- switch(vl[[i]]$dataType,
                              methylation     = cbind(vl[[i]]$Des[, "geneSymbolId", drop = FALSE],
                                                      Platform = rep("methylation"    ,
                                                                     dim(vl[[i]]$Des)[1]),
                                                      Description = vl[[i]]$Des[, 2]),
                              copyNumber      = cbind(vl[[i]]$Des[, "geneSymbolId", drop = FALSE],
                                                      Platform = rep("copyNumber"     ,
                                                                     dim(vl[[i]]$Des)[1]),
                                                      Description = vl[[i]]$Des[, 2]),
                              somaticMutation = cbind(vl[[i]]$Des[, "geneSymbolId", drop = FALSE],
                                                      Platform = rep("somaticMutation",
                                                                     dim(vl[[i]]$Des)[1]),
                                                      Description = vl[[i]]$Des[, 2]),
                              geneExp         = cbind(vl[[i]]$Des[, "geneSymbolId", drop = FALSE],
                                                      Platform = rep("geneExp"        ,
                                                                     dim(vl[[i]]$Des)[1]),
                                                      Description = vl[[i]]$Des[, 2]),
                              miRNAExp        = cbind(vl[[i]]$Des[, "geneSymbolId", drop = FALSE],
                                                      Platform = rep("miRNAExp"       ,
                                                                     dim(vl[[i]]$Des)[1]),
                                                      Description = vl[[i]]$Des[, 2]),
                              protein_RPPA    = cbind(vl[[i]]$Des[, "geneSymbolId", drop = FALSE],
                                                      Platform = rep("protein_RPPA"   ,
                                                                     dim(vl[[i]]$Des)[1]),
                                                      Description = vl[[i]]$Des[, 2]),
                              protein_iTRAQ   = cbind(vl[[i]]$Des[, "geneSymbolId", drop = FALSE],
                                                      Platform = rep("protein_iTRAQ"  ,
                                                                     dim(vl[[i]]$Des)[1]),
                                                      Description = vl[[i]]$Des[, 2]))
        # check NA gene. which(vl[[i]]$Des[, "geneSymbolId"] != "NA")
        ID <- intersect(which(!is.na(vl[[i]]$Des[, "geneSymbolId"])),
                        which(vl[[i]]$Des[, "geneSymbolId"] != "?"))
        vl[[i]]$Des <- vl[[i]]$Des[ID, , drop = FALSE]
        vl[[i]]$Data <- vl[[i]]$Data[ID, , drop = FALSE]
    }
    for (i in 1:length(vl)) {
        # Get the data of samples that should be kept
        if (combineStyle == "Intersect") {
            vl[[i]]$Data <- vl[[i]]$Data[, sampleId, drop = FALSE]
        }
        if (combineStyle == "Union") {
            tempData <- matrix(NA,
                               dim(vl[[i]]$Data)[1],
                               length(sampleId))
            colnames(tempData) <- sampleId
            tempData[, colnames(vl[[i]]$Data)] <- vl[[i]]$Data
            vl[[i]]$Data <- tempData
        }
        rownames(vl[[i]]$Data) <- NULL
        rownames(vl[[i]]$Des) <- NULL
    }
    # Combine the datasets into matrix format
    Data <- vl[[1]]$Data
    Des <- vl[[1]]$Des
    for (i in 2:length(vl)) {
        Data <- rbind(Data, vl[[i]]$Data)
        Des  <- rbind(Des , vl[[i]]$Des)
    }
    OrderID <- order(as.character(Des[, "geneSymbolId"]),
                     as.character(Des[, "Platform"]),
                     as.character(Des[, "Description"]),
                     na.last = TRUE,
                     decreasing = FALSE)
    Data <- Data[OrderID, , drop = FALSE]
    Des <- Des[OrderID, , drop = FALSE]
    rownames(Data) <- NULL
    rownames(Des) <- NULL
    # filter the all zero line of mutation
    vbDes <- Des[, "Platform"] == "somaticMutation"
    vbData <- apply(Data,
                    1,
                    function(x){
                        all(unique(x) == 0)
                    })
    vb <- !(vbDes & vbData)
    Data <- Data[vb, , drop = FALSE]
    Des  <-  Des[vb, , drop = FALSE]
    vbDataRowNA <- apply(Data,
                         1,
                         function(x){all(is.na(x))})
    Data <- Data[!vbDataRowNA, , drop = FALSE]
    Des  <-  Des[!vbDataRowNA, , drop = FALSE]
    mode(Data) <- "numeric"
    lResult <- list(Des = Des,
                    Data = Data)
    options(warn = 0)
    lResult
}

ExtractTissueSpecificSamples <- function(inputData,
                                         tissueType,
                                         singleSampleFlag,
                                         sampleTypeFile = "./SupportingFiles/TCGASampleType.txt") {
    options(warn = -1)
    TCGABarcode <- colnames(inputData)
    TCGABarcode <- sapply(strsplit(TCGABarcode,
                                   split = "-"),
                          function(x){
                              Str = paste(x[1],
                                          x[2],
                                          x[3],
                                          substr(x[4], 1, 2),
                                          sep = "-")
                              return(Str)
                          })
    if (singleSampleFlag == TRUE) {
        DuplicatedLabel <- duplicated(TCGABarcode)
        ID <- which(DuplicatedLabel == FALSE)
        inputData <- inputData[, ID, drop = FALSE]
        TCGABarcode <- TCGABarcode[ID]
    }
    TCGADataLabel <- sapply(strsplit(TCGABarcode,
                                     split = "-"),
                            function(x){
                                return(x[4])
                            })
    TCGADataLabelValue <- as.numeric(substr(TCGADataLabel, 1, 2))
    sampleType <- read.table(file = sampleTypeFile,
                             header = TRUE,
                             sep = "\t",
                             na.string = "",
                             quote = "",
                             stringsAsFactors = FALSE)
    Code <- sampleType[which(sampleType[, "Options"] %in% tissueType), "Code"]
    if (length(Code) == 0) {
        writeLines("Error: no valid tissue or cell type specified")
        inputData <- NULL
    }
    else {
        ID <- which(TCGADataLabelValue %in% Code)
        if (length(ID) == 0) {
            writeLines(paste("There is no sample of ",
                             paste(tissueType,
                                   collapse = ", "),
                             sep = ""))
            inputData <- NULL
        }
        else {
            inputData <- inputData[, ID, drop = FALSE]
            writeLines(paste(length(ID),
                             " samples of ",
                             paste(tissueType,
                                   collapse = ", "),
                             " are extracted.",
                             sep = ""))
        }
    }
    options(warn = 0)
    inputData
}

CalculateSingleValueMethylationData <- function(input,
                                                regionOption,
                                                DHSOption,
                                                outputFileName,
                                                outputFileFolder = ".",
                                                chipAnnotationFile = "./SupportingFiles/dMethyInfo.rda") {
    options(warn = -1)
    if (outputFileFolder != ".") {dir.create(outputFileFolder, recursive = TRUE)}
    if (sum(! regionOption %in% c("TSS1500",
                                  "TSS200",
                                  "5'UTR",
                                  "1stExon",
                                  "Body",
                                  "3'UTR",
                                  "All")) > 0) {
        stop("Error: invalid regionOption.")
    }
    if (sum(! DHSOption %in% c("DHS",
                               "notDHS",
                               "Both")) > 0) {
        stop("Error: invalid DHSOption")
    }
    regionOptionO <- regionOption
    #   if ("TSS1500" %in% regionOption)
    #   {
    #     regionOption <- unique(c(regionOption, "TSS200"))
    #   }
    #
    if ((!("All" %in% regionOption)) || (DHSOption != "Both")) {
        load(chipAnnotationFile)
        ID <- 1:dim(dMethyInfo)[1]
        if (!("All" %in% regionOption)) {
            ID <- sort(unique(which(dMethyInfo[,
                                    "UCSC_RefGene_Group"] %in% regionOption)))
        }
        if (DHSOption == "DHS") {
            ID <- sort(intersect(ID,
                                 which(dMethyInfo[, "DHS"] == TRUE)))
        }
        else {
            if (DHSOption == "notDHS") {
                ID <- sort(intersect(ID,
                                     setdiff(1:dim(dMethyInfo)[1],
                                             which(dMethyInfo[, "DHS"] == TRUE))))
            }
        }
        dMethyInfo <- dMethyInfo[ID, , drop = FALSE]
        ID <- which(input$Des[, "REF"] %in% dMethyInfo[, "IlmnID"])
        #ID <- which(input$Des[, "REF"] %in% rownames(dMethyInfo))
        input$Des <- input$Des[ID, , drop = FALSE]
        input$Data <- input$Data[ID, , drop = FALSE]
    }
    # Check whether conflict with NA gene
    ID <- sort(intersect(which(!is.na(input$Des[, "geneSymbolId"])),
                         which(input$Des[, "geneSymbolId"] != "")),
               decreasing = FALSE)
    input$Des <- input$Des[ID, , drop = FALSE]
    input$Data <- input$Data[ID, , drop = FALSE]
    UniGene <- unique(input$Des[, "geneSymbolId"])
    Data <- matrix(NA, length(UniGene), dim(input$Data)[2])
    for (i in 1:length(UniGene)) {
        if ((i %% 2000 == 0) || (i == length(UniGene))) {
            writeLines(paste("Calculation ",
                             round(i/length(UniGene)*100,
                                   digits = 2),
                             "% done.",
                             sep = ""))
        }
        Data[i, ] <- colMeans(input$Data[which(input$Des[, "geneSymbolId"] == UniGene[i]), , drop = FALSE],
                              na.rm = TRUE)
    }
    colnames(Data) <- colnames(input$Data)
    Des <- cbind(geneSymbolId = UniGene,
                 SingleValueType = rep(paste(paste(regionOptionO,
                                                   collapse = "-"),
                                             DHSOption,
                                             sep = "|"),
                                       length(UniGene)))
    # draw and save a box plot of data
    #:png(filename = paste(outputFileFolder,
    #:                     ifelse(endsWith(outputFileFolder, "/"), "", "/"),
    #:                     outputFileName,
    #:                     "__",
    #:                     paste(regionOptionO,
    #:                           collapse = "-"),
    #:                     "__",
    #:                     DHSOption,
    #:                     "__boxplot.png",
    #:                     sep = ""),
    #:    width = 30*dim(Data)[2]+300,
    #:    height = 1500,
    #:    units = "px")
    #:par(las = 2,
    #:    mai = c(4.5, 1, 0.5, 1),
    #:    cex = 1.5)
    #:boxplot(Data)
    #:dev.off()
    # save data files.
    rownames(Data) <- NULL
    rownames(Des) <- NULL
    ID <- which(is.nan(Data))
    Data[ID] <- NA
    save(Des, Data,
         file = paste(outputFileFolder,
                      ifelse(endsWith(outputFileFolder, "/"), "", "/"),
                      outputFileName,
                      "__",
                      paste(regionOptionO,
                            collapse = "-"),
                      "__",
                      DHSOption,
                      ".rda",
                      sep = ""))
    #:write.table(cbind(Des, Data),
    #:            file = paste(outputFileFolder,
    #:                         ifelse(endsWith(outputFileFolder, "/"), "", "/"),
    #:                         outputFileName,
    #:                         "__",
    #:                         paste(regionOptionO,
    #:                               collapse = "-"),
    #:                         "__",
    #:                         DHSOption,
    #:                         ".txt",
    #:                         sep = ""),
    #:            quote = FALSE,
    #:            sep = "\t",
    #:            na = "",
    #:            col.names = TRUE,
    #:            row.names = FALSE)
    options(warn = 0)
    return(list(Data = Data,
                Des = Des))
}

MergeMethylationData <- function(input1,
                                 input2,
                                 outputFileName,
                                 outputFileFolder = ".") {
    options(warn = -1)
    if (outputFileFolder != ".") {dir.create(outputFileFolder, recursive = TRUE)}
    input1_O <- input1
    input2_O <- input2
    rownames(input1$Des) <- paste(input1$Des[, "REF"],
                                  input1$Des[, "geneSymbolId"],
                                  sep = "||")
    rownames(input1$Data) <- paste(input1$Des[, "REF"],
                                   input1$Des[, "geneSymbolId"],
                                   sep = "||")
    rownames(input2$Des) <- paste(input2$Des[, "REF"],
                                  input2$Des[, "geneSymbolId"],
                                  sep = "||")
    rownames(input2$Data) <- paste(input2$Des[, "REF"],
                                   input2$Des[, "geneSymbolId"],
                                   sep = "||")
    CommonID <- intersect(rownames(input1$Des),
                          rownames(input2$Des))
    input1$Des <- input1$Des[CommonID, , drop = FALSE]
    input1$Data <- input1$Data[CommonID, , drop = FALSE]
    input2$Des <- input2$Des[CommonID, , drop = FALSE]
    input2$Data <- input2$Data[CommonID, , drop = FALSE]
    if (dim(input1_O$Des)[1] >= dim(input2_O$Des)[1]) {
        Des <- input1$Des
    }
    else {
        Des <- input2$Des
    }
    Data <- cbind(input1$Data,
                  input2$Data)
    # draw and save a box plot of combined data before normalization
    #:png(filename = paste(outputFileFolder,
    #:                     ifelse(endsWith(outputFileFolder, "/"), "", "/"),
    #:                     outputFileName,
    #:                     "__BeforeNormalizationBoxplot.png",
    #:                     sep = ""),
    #:    width = 30*dim(Data)[2]+300,
    #:    height = 1500,
    #:    units = "px")
    #:par(las = 2,
    #:    mai = c(4.5, 1, 0.5, 1),
    #:    cex = 1.5)
    #:boxplot(Data)
    #:dev.off()
    Data <- Normalize.quantiles(Data)
    # draw and save a box plot of combined data after normalization
    #:png(filename = paste(outputFileFolder,
    #:                     ifelse(endsWith(outputFileFolder, "/"), "", "/"),
    #:                     outputFileName,
    #:                     "__AfterNormalizationBoxplot.png",
    #:                     sep = ""),
    #:    width = 30*dim(Data)[2]+300,
    #:    height = 1500,
    #:    units = "px")
    #:par(las = 2,
    #:    mai = c(4.5, 1, 0.5, 1),
    #:    cex = 1.5)
    #:boxplot(Data)
    #:dev.off()
    # save data files.
    rownames(Data) <- NULL
    rownames(Des) <- NULL
    save(Des, Data,
         file = paste(outputFileFolder,
                      ifelse(endsWith(outputFileFolder, "/"), "", "/"),
                      outputFileName,
                      ".rda",
                      sep = ""))
    #:write.table(cbind(Des, Data),
    #:            file = paste(outputFileFolder,
    #:                         ifelse(endsWith(outputFileFolder, "/"), "", "/"),
    #:                         outputFileName,
    #:                         ".txt",
    #:                         sep = ""),
    #:            quote = FALSE,
    #:            sep = "\t",
    #:            na = "",
    #:            col.names = TRUE,
    #:            row.names = FALSE)
    options(warn = 0)
    return(list(Data = Data,
                Des = Des))
}

ProcessCNAData <- function(inputFilePath,
                           outputFileName,
                           outputFileFolder = ".",
                           refGenomeFile = "SupportingFiles/dGeneHg19.rda") {
    options(warn = -1)
    if (outputFileFolder != ".") {dir.create(outputFileFolder, recursive = TRUE)}
    # Load reference genomes, which give the genomic locations of genes
    #dGeneHg19 <- read.table(file = refGenomeFile, header = TRUE, sep = "\t", na.string = "", quote = "", stringsAsFactors = FALSE)
    load(refGenomeFile)
    NumGene <- dim(dGeneHg19)[1]
    dGeneHg19$chr    <- toupper(as.character(dGeneHg19$chr))
    dGeneHg19$strand <- as.character(dGeneHg19$strand)
    # Load copy number data
    #InData <- read.table(file = inputFilePath, header = TRUE, sep = "\t", na.string = "", quote = "", stringsAsFactors = FALSE, check.names = FALSE)
    load(gsub("\\.txt",
              ".rda",
              inputFilePath))
    InData <- dPiped
    InData[, 1] <- as.character(InData[, 1])
    InData[, "Sample"] <- toupper(InData[, "Sample"])
    TCGAID <- InData[, "Sample"]
    UniqueTCGAID <- unique(TCGAID)
    InData[, "Chromosome"] <- as.character(InData[, "Chromosome"])
    InData[, "Chromosome"] <- paste("CHR", InData[, "Chromosome"],
                                    sep = "")
    IDID <- which(InData[, "Chromosome"] == "CHR23")
    InData[IDID, "Chromosome"] <- "CHRX"
    IDID <- which(InData[, "Chromosome"] == "CHR24")
    InData[IDID, "Chromosome"] <- "CHRY"
    InData[, "Start"] <- as.numeric(InData[, "Start"])
    InData[, "End"] <- as.numeric(InData[, "End"])
    InData[, "Segment_Mean"] <- as.numeric(InData[, "Segment_Mean"])
    # Calculate gene-level copy number
    Data <- matrix(NA, NumGene, length(UniqueTCGAID))
    for (i in 1:length(UniqueTCGAID)) {
        IDi <- which(InData[, "Sample"] == UniqueTCGAID[i])
        InDatai <- InData[IDi, c("Chromosome", "Start", "End", "Segment_Mean"), drop = FALSE]
        for (j in 1:NumGene) {
            IDj <- which(InDatai[, "Chromosome"] == dGeneHg19[j, "chr"])
            if (length(IDj) > 0) {
                StartP <- InDatai[IDj, "Start"]
                EndP <- InDatai[IDj, "End"]
                Ratio <- InDatai[IDj, "Segment_Mean"]
                Overlap <- ifelse((EndP >= dGeneHg19[j, "start"]) & (dGeneHg19[j, "end"] >= StartP),
                                  pmin(EndP, dGeneHg19[j, "end"]) - pmax(StartP, dGeneHg19[j, "start"]) + 1,
                                  0)
                IDNonNA <- which(!is.na(Ratio))
                Ratio <- Ratio[IDNonNA]
                Overlap <- Overlap[IDNonNA]
                if (sum(Overlap, na.rm = TRUE) > 0) {
                    Data[j, i] = sum(Overlap * Ratio, na.rm = TRUE) / sum(Overlap, na.rm = TRUE)
                }
            }
        }
        if ((i %% 20 == 0) & (i < length(UniqueTCGAID))) {
            writeLines(paste("Calculating gene copy number, ",
                             round(i/length(UniqueTCGAID)*100,
                                   digits = 2),
                             "% done.",
                             sep = ""))
        }
    }
    writeLines("Calculating gene copy number, 100% done.")
    colnames(Data) <- UniqueTCGAID
    Des <- as.matrix(dGeneHg19[, 1:3, drop = FALSE])
    #lResult <- CombineRedundantFeature(Data = Data, Des = Des)
    lResult <- list(Data = Data, Des = Des)
    Data <- lResult$Data
    Des <- lResult$Des
    lRow <- split(seq(dim(Data)[1]),
                  Des[, 1])
    lData <- lapply(lRow,
                    function(x) {
                        if (length(x) == 1) {
                            Data[x, ]
                        } else {
                            colMeans(Data[x, ])
                    }})
    Data <- do.call(rbind, lData)
    load("./SupportingFiles/lGeneInfo.rda")
    load("./SupportingFiles/dGeneHg19.rda")
    sSymbolId2ChrMap <- as.character(dGeneHg19$strand)
    names(sSymbolId2ChrMap) <- dGeneHg19$geneSymbolId
    Des <- as.matrix(cbind(names(lData),
                           paste("chr",
                                 lGeneInfo$sId2ChrMap[sapply(strsplit(names(lData),
                                                                      split = "\\|"),
                                                             function(x){x[2]})],
                                 sep = ""),
                           sSymbolId2ChrMap[names(lData)]))
    colnames(Des) <- c("geneSymbolId", "chr", "strand")
    # Draw and save a box plot
    #:png(filename = paste(outputFileFolder,
    #:                     ifelse(endsWith(outputFileFolder, "/"), "", "/"),
    #:                     outputFileName,
    #:                     "__boxplot.png",
    #:                     sep = ""),
    #:    width = 30*dim(Data)[2]+300,
    #:    height = 1500,
    #:    units = "px")
    #:par(las = 2,
    #:    mai = c(4.5, 1, 0.5, 1),
    #:    cex = 1.5)
    #:boxplot(Data)
    #:dev.off()
    # save data files.
    rownames(Data) <- NULL
    rownames(Des) <- NULL
    save(Des, Data,
         file = paste(outputFileFolder,
                      ifelse(endsWith(outputFileFolder, "/"), "", "/"),
                      outputFileName,
                      ".rda",
                      sep = ""))
    #:write.table(cbind(Des, Data),
    #:            file = paste(outputFileFolder,
    #:                         ifelse(endsWith(outputFileFolder, "/"), "", "/"),
    #:                         outputFileName,
    #:                         ".txt",
    #:                         sep = ""),
    #:            quote = FALSE,
    #:            sep = "\t",
    #:            na = "",
    #:            col.names = TRUE,
    #:            row.names = FALSE)
    options(warn = 0)
    return(list(Data = Data,
                Des = Des))
}

ProcessRNASeqData <- function(inputFilePath,
                              outputFileName,
                              outputFileFolder = ".",
                              dataType,
                              verType = "RNASeqV2") {
    options(warn = -1)
    if (outputFileFolder != ".") {dir.create(outputFileFolder, recursive = TRUE)}
    # Read in data.
    #InData <- read.table(file = inputFilePath, header = TRUE, sep = "\t", na.string = "", stringsAsFactors = FALSE, quote = "", check.names = FALSE)
    load(gsub("\\.txt",
              ".rda",
              inputFilePath))
    InData <- dPiped
    InData[, 1] <- as.character(InData[, 1])
    #InData <- InData[2:dim(InData)[1], , drop = FALSE]
    if ((dataType == "geneExp") & (verType == "RNASeqV1")) {
        REF <- InData[, 1]
        RPKM <- as.matrix(InData[, seq(4, dim(InData)[2], 3), drop = FALSE])
        mode(RPKM) <- "numeric"
        geneSymbolId <- sapply(strsplit(REF,
                                        split = "\\|"),
                               function(x){x[1]})
        EntrezID <- sapply(strsplit(REF,
                                    split = "\\|"),
                           function(x){x[2]})
        Des <- cbind(geneSymbolId = geneSymbolId,
                     EntrezID = EntrezID)
        Data <- RPKM
    }
    if (dataType == "exonExp") {
        REF <- InData[, 1]
        # RPKM <- as.matrix(InData[, seq(4, dim(InData)[2], 3), drop = FALSE])
        RPKM <- as.matrix(InData[, 2:dim(InData)[2], drop = FALSE])
        mode(RPKM) <- "numeric"
        Des <- cbind(ExonID = REF)
        Data <- RPKM
    }
    if ((dataType == "geneExp") & (verType == "RNASeqV2")) {
        REF <- InData[, 1]
        NormalizedCount <- as.matrix(InData[, 2:dim(InData)[2], drop = FALSE])
        mode(NormalizedCount) <- "numeric"
        #geneSymbolId <- sapply(strsplit(REF, split = "\\|"), function(x)x[1])
        geneSymbolId <- REF
        EntrezID <- sapply(strsplit(REF,
                                    split = "\\|"),
                           function(x){x[2]})
        Des <- cbind(geneSymbolId = geneSymbolId,
                     EntrezID = EntrezID)
        Data <- NormalizedCount
    }
    if (verType == "Microarray") {
        dataType <- "geneExp"
        REF <- InData[, 1]
        Data <- as.matrix(InData[, 2:dim(InData)[2], drop = FALSE])
        mode(Data) <- "numeric"
        EntrezID <- sapply(strsplit(REF,
                                    split = "\\|"),
                           function(x)x[2])
        Des <- cbind(geneSymbolId = geneSymbolId,
                     EntrezID = EntrezID)
    }
    #lResult <- CombineRedundantFeature(Data = Data, Des = Des)
    lResult <- list(Data = Data,
                    Des = Des)
    Data <- lResult$Data
    Des <- lResult$Des
    if (dataType == "geneExp") {
        # Draw and save a box plot
        #:png(filename = paste(outputFileFolder,
        #:                     ifelse(endsWith(outputFileFolder, "/"), "", "/"),
        #:                     outputFileName,
        #:                     "__boxplot.png",
        #:                     sep = ""),
        #:    width = 30*dim(Data)[2]+300,
        #:    height = 1500,
        #:    units = "px")
        #:par(las = 2,
        #:    mai = c(4.5, 1, 0.5, 1),
        #:    cex = 1.5)
        #:nonNaID <- which(!is.na(Data))
        #:if ((sum(Data[nonNaID]<0) == 0) && (max(Data[nonNaID])>50)) {
        #:  boxplot(log2(Data),
        #:          main = "Boxplot drawn based on log2 tranformed data.")
        #:} else {
        #:  boxplot(Data)
        #:}
        #:dev.off()
    }
    # save data files
    rownames(Data) <- NULL
    rownames(Des) <- NULL
    save(Des, Data,
         file = paste(outputFileFolder,
                      ifelse(endsWith(outputFileFolder, "/"), "", "/"),
                      outputFileName,
                      ".rda",
                      sep = ""))
    #:write.table(cbind(Des,
    #:                  Data),
    #:            file = paste(outputFileFolder,
    #:                         ifelse(endsWith(outputFileFolder, "/"), "", "/"),
    #:                         outputFileName,
    #:                         ".txt",
    #:                         sep = ""),
    #:            quote = FALSE,
    #:            sep = "\t",
    #:            na = "",
    #:            col.names = TRUE,
    #:            row.names = FALSE)
    options(warn = 0)
    return(list(Data = Data,
                Des = Des))
}

ProcessmiRNASeqData <- function(inputFilePath,
                                outputFileName,
                                outputFileFolder = ".",
                                fileSource = "TCGA-Assembler") {
    options(warn = -1)
    if (outputFileFolder != ".") {dir.create(outputFileFolder, recursive = TRUE)}
    # Read in data.
    #InData <- read.table(file = inputFilePath, header = TRUE, sep = "\t", na.string = "", quote = "", stringsAsFactors = FALSE, check.names = FALSE)
    load(gsub("\\.txt",
              ".rda",
              inputFilePath))
    #dPiped <- dGt23k
    InData <- dPiped
    InData[, 1] <- as.character(InData[, 1])
    #InData <- InData[2:dim(InData)[1], , drop = FALSE]
    # divide read cout data and RPM data
    REF <- InData[, 1]
    if (fileSource == "Firehose") {
        Count <- as.matrix(InData[, seq(2, dim(InData)[2], 3), drop = FALSE])
        RPM <- as.matrix(InData[, seq(3, dim(InData)[2], 3), drop = FALSE])
    }
    if (fileSource == "TCGA-Assembler") {
        Count <- as.matrix(InData[, seq(2, dim(InData)[2], 2), drop = FALSE])
        RPM <- as.matrix(InData[, seq(3, dim(InData)[2], 2), drop = FALSE])
        colnames(Count) <- sapply(strsplit(colnames(Count),
                                           split = ":"),
                                  function(x){x[1]})
        colnames(RPM  ) <- sapply(strsplit(colnames(RPM  ),
                                           split = ":"),
                                  function(x){x[1]})
        #colnames(Count) <- sapply(strsplit(InData[, seq(2, dim(InData)[2], 2), drop = FALSE], split = ":"), function(x){x[1]})
        #colnames(RPM)   <- sapply(strsplit(InData[, seq(3, dim(InData)[2], 2), drop = FALSE], split = ":"), function(x){x[1]})
    }
    #  rownames(Count) <- REF
    mode(Count) <- "numeric"
    #  rownames(RPM) <- REF
    mode(RPM) <- "numeric"
    # Draw and save box plot
    #:   png(filename = paste(outputFileFolder,
    #:                        ifelse(endsWith(outputFileFolder, "/"), "", "/"),
    #: outputFileName, "__ReadCount.boxplot.png", sep = ""), width = 30*dim(Count)[2]+300, height = 1500, units = "px")
    #:   par(las = 2, mai = c(4.5, 1, 0.5, 1), cex = 1.5)
    #:   boxplot(log2(Count))
    #:   dev.off()
    #:   png(filename = paste(outputFileFolder,
    #:                     ifelse(endsWith(outputFileFolder, "/"), "", "/"),
    #: outputFileName, "__RPM.boxplot.png", sep = ""), width = 30*dim(RPM)[2]+300, height = 1500, units = "px")
    #:   par(las = 2, mai = c(4.5, 1, 0.5, 1), cex = 1.5)
    #:   boxplot(log2(RPM))
    #:   dev.off()
    # save data files
    Des <- cbind(geneSymbolId = REF,
                 EntrezID = sapply(strsplit(REF,
                                            split = "\\|"),
                                   function(x){x[2]}))
    Data <- Count
    #lResult <- CombineRedundantFeature(Data = Data, Des = Des)
    lResult <- list(Data = Data,
                    Des = Des)
    Data <- lResult$Data
    Des <- lResult$Des
    ID <- which(is.na(Des[, "geneSymbolId"]))
    Des[ID, "geneSymbolId"] <- ""
    rownames(Data) <- NULL
    rownames(Des) <- NULL
    save(Des, Data,
         file = paste(outputFileFolder,
                      ifelse(endsWith(outputFileFolder, "/"), "", "/"),
                      outputFileName,
                      "__ReadCount.rda",
                      sep = ""))
    #:write.table(cbind(Des,
    #:                  Data),
    #:            file = paste(outputFileFolder,
    #:                         ifelse(endsWith(outputFileFolder, "/"), "", "/"),
    #:                         outputFileName,
    #:                         "__ReadCount.txt",
    #:                         sep = ""),
    #:            quote = FALSE,
    #:            sep = "\t",
    #:            na = "",
    #:            col.names = TRUE,
    #:            row.names = FALSE)
    Des <- cbind(geneSymbolId = REF,
                 EntrezID = sapply(strsplit(REF,
                                            split = "\\|"),
                                   function(x){x[2]}))
    Data <- RPM
    #lResult <- CombineRedundantFeature(Data = Data, Des = Des)
    lResult <- list(Data = Data,
                    Des = Des)
    Data <- lResult$Data
    Des <- lResult$Des
    ID <- which(is.na(Des[, "geneSymbolId"]))
    Des[ID, "geneSymbolId"] <- ""
    rownames(Data) <- NULL
    rownames(Des) <- NULL
    save(Des, Data,
         file = paste(outputFileFolder,
                      ifelse(endsWith(outputFileFolder, "/"), "", "/"),
                      outputFileName,
                      "__RPM.rda",
                      sep = ""))
    #:write.table(cbind(Des,
    #:                  Data),
    #:            file = paste(outputFileFolder,
    #:                         ifelse(endsWith(outputFileFolder, "/"), "", "/"),
    #:                         outputFileName,
    #:                         "__RPM.txt",
    #:                         sep = ""),
    #:            quote = FALSE,
    #:            sep = "\t",
    #:            na = "",
    #:            col.names = TRUE,
    #:            row.names = FALSE)
    options(warn = 0)
    return(list(Data = Data,
                Des = Des))
}

ProcessMethylation450Data <- function(inputFilePath,
                                      outputFileName,
                                      outputFileFolder = ".",
                                      fileSource = "TCGA-Assembler") {
    options(warn = -1)
    if (outputFileFolder != ".") {dir.create(outputFileFolder, recursive = TRUE)}
    writeLines("Loading data.")
    if (fileSource == "Firehose") {
        # read in data
        GetMethylation450Data(inputFilePath,
                              outputFileName,
                              outputFileFolder)
        load(paste(outputFileFolder,
                   ifelse(endsWith(outputFileFolder, "/"), "", "/"),
                   outputFileName,
                   ".TempFiles/",
                   outputFileName,
                   ".TempData.rda",
                   sep = ""))
        unlink(x = paste(outputFileFolder,
                         ifelse(endsWith(outputFileFolder, "/"), "", "/"),
                         outputFileName,
                         ".TempFiles",
                         sep = ""),
               recursive = TRUE,
               force = TRUE)
        Methy450Des <- as.matrix(Methy450Des)
    }
    if (fileSource == "TCGA-Assembler") {
        # read in data
        #InData <- read.table(file = inputFilePath, header = TRUE, sep = "\t", na.string = "", quote = "", stringsAsFactors = FALSE, check.names = FALSE)
        load(gsub("\\.txt",
                  ".rda",
                  inputFilePath))
        InData <- dPiped
        InData[, 1] <- as.character(InData[, 1])
        Methy450Des <- as.matrix(InData[, 1:4, drop = FALSE])
        colnames(Methy450Des) <- c("REF",
                                   "geneSymbolId",
                                   "ChromosomeID",
                                   "CoordinateID")
        ID <- which(Methy450Des[, "geneSymbolId"] == "DISCONTUNUED")
        Methy450Des[ID, "geneSymbolId"] <- ""
        Methy450Data <- as.matrix(InData[, 5:dim(InData)[2], drop = FALSE])
        mode(Methy450Data) <- "numeric"
        rm("InData")
    }
    # some probes have more than one gene symbols, collpase the data
    # AddDes <- matrix("", dim(Methy450Des)[1], 4)
    AddDes <- matrix("",
                     500000,
                     4)
    colnames(AddDes) <- colnames(Methy450Des)
    # AddData <- matrix(NA, dim(Methy450Data)[1], dim(Methy450Data)[2])
    AddData <- matrix(NA,
                      500000,
                      dim(Methy450Data)[2])
    colnames(AddData) <- colnames(Methy450Data)
    AddIndex <- 0
    for (ProbeID in 1:dim(Methy450Des)[1]) {
        GrepID <- grep(";",
                       Methy450Des[ProbeID, "geneSymbolId"])
        if ((ProbeID %% 50000 == 0) || (ProbeID == dim(Methy450Des)[1])) {
            writeLines(paste("Process CpG sites corresponding to multiple genes, ",
                             ceiling(ProbeID/dim(Methy450Des)[1]*100),
                             "% done.",
                             sep = ""))
        }
        if (length(GrepID) > 0) {
            TempStr <- strsplit(Methy450Des[ProbeID, "geneSymbolId"],
                                split = ";")
            Methy450Des[ProbeID, "geneSymbolId"] <- TempStr[[1]][1]
            for (TempStri in 2:length(TempStr[[1]])) {
                TempAddDes <- Methy450Des[ProbeID, , drop = FALSE]
                TempAddDes[1, "geneSymbolId"] <- TempStr[[1]][TempStri]
                AddIndex <- AddIndex + 1
                AddDes[AddIndex, ] <- TempAddDes
                AddData[AddIndex, ] <- Methy450Data[ProbeID, ]
            }
        }
    }
    # writeLines(paste("Number of added rows is ", AddIndex, sep = ""))
    Methy450Des <- rbind(Methy450Des,
                         AddDes[1:AddIndex, , drop = FALSE])
    Methy450Data <- rbind(Methy450Data, AddData[1:AddIndex, , drop = FALSE])
    Methy450Des[, "ChromosomeID"] <- gsub(" ",
                                          "",
                                          Methy450Des[, "ChromosomeID"])
    Methy450Des[, "CoordinateID"] <- gsub(" ",
                                          "",
                                          Methy450Des[, "CoordinateID"])
    OrderID <- order(as.numeric(Methy450Des[, "ChromosomeID"]),
                     as.numeric(Methy450Des[, "CoordinateID"]),
                     Methy450Des[, "geneSymbolId"],
                     na.last = TRUE,
                     decreasing = FALSE)
    Des <- Methy450Des[OrderID, , drop = FALSE]
    Data <- Methy450Data[OrderID, , drop = FALSE]
    #lResult <- CombineRedundantFeature(Data = Data, Des = Des)
    lResult <- list(Data = Data,
                    Des = Des)
    Data <- lResult$Data
    Des <- lResult$Des
    # Draw and save box plot
    #:png(filename = paste(outputFileFolder,
    #:                     ifelse(endsWith(outputFileFolder, "/"), "", "/"),
    #:                     outputFileName,
    #:                     "__boxplot.png",
    #:                     sep = ""),
    #:    width = 30*dim(Data)[2]+300,
    #:    height = 1500,
    #:    units = "px")
    #:par(las = 2,
    #:    mai = c(4.5, 1, 0.5, 1),
    #:    cex = 1.5)
    #:boxplot(Data)
    #:dev.off()
    # save output data files
    rownames(Data) <- NULL
    rownames(Des) <- NULL
    save(Des, Data,
         file = paste(outputFileFolder,
                      ifelse(endsWith(outputFileFolder, "/"), "", "/"),
                      outputFileName,
                      ".rda",
                      sep = ""))
    #:write.table(cbind(Des,
    #:                  Data),
    #:            file = paste(outputFileFolder,
    #:                     ifelse(endsWith(outputFileFolder, "/"), "", "/"),
    #:                         outputFileName,
    #:                         ".txt",
    #:                         sep = ""),
    #:            quote = FALSE,
    #:            sep = "\t",
    #:            na = "",
    #:            col.names = TRUE,
    #:            row.names = FALSE)
    options(warn = 0)
    return(list(Data = Data,
                Des = Des))
}

ProcessRPPADataWithGeneAnnotation <- function(inputFilePath,
                                              outputFileName,
                                              outputFileFolder = ".") {
    options(warn = -1)
    if (outputFileFolder != ".") {dir.create(outputFileFolder, recursive = TRUE)}
    # Read in data
    #InData <- read.table(file = inputFilePath, header = TRUE, sep = "\t", na.string = "", quote = "", stringsAsFactors = FALSE, check.names = FALSE)
    load(gsub("\\.txt",
              ".rda",
              inputFilePath))
    InData <- dPiped
    InData[, 1] <- as.character(InData[, 1])
    #OutData <- matrix(NA, 0, dim(InData)[2]-1)
    ## Split the gene symbol and protein antibody into two separate columns
    ## Duplicate rows for proteins coded by multiple genes.
    #Gene <- c()
    #Antibody <- c()
    #for (i in 1:dim(InData)[1]) {
    #  GeneAntibodyStr <- InData[i, 1]
    #  GeneStr <- strsplit(GeneAntibodyStr, split = "\\|")[[1]][c(1,2)]
    #  AntibodyStr <- strsplit(GeneAntibodyStr, split = "\\|")[[1]][3]
    #  GeneStr <- unlist(strsplit(GeneStr, split = " "))
    #  Gene <- c(Gene, GeneStr)
    #  Antibody <- c(Antibody, rep(AntibodyStr, length = length(GeneStr)))
    #  OutData <- rbind(OutData, InData[rep(i, length = length(GeneStr)), 2:dim(InData)[2], drop = FALSE])
    #}
    #colnames(OutData) <- colnames(InData)[2:dim(InData)[2]]
    #Des <- cbind(geneSymbolId = Gene, ProteinAntibody = Antibody)
    #Data <- as.matrix(OutData)
    #mode(Data) <- "numeric"
    lGeneAb <- strsplit(InData[, 1],
                        split = "\\|")
    Des <- cbind(geneSymbolId = sapply(lGeneAb,
                                       function(x){
                                           return(paste(x[seq(2)],
                                                        collapse = "|"))
                                       }),
                 Antibody = sapply(lGeneAb,
                                   function(x){
                                       return(x[3])
                                       #return(paste(x[seq(3)],
                                       #             collapse = "|"))
                                   }))
    Data <- as.matrix(InData[, -1])
    #lResult <- CombineRedundantFeature(Data = Data, Des = Des)
    lResult <- list(Data = Data,
                    Des = Des)
    Data <- lResult$Data
    Des <- lResult$Des
    mode(Data) <- "numeric"
    # save output data files
    rownames(Data) <- NULL
    rownames(Des) <- NULL
    save(Des, Data,
         file = paste(outputFileFolder,
                      ifelse(endsWith(outputFileFolder, "/"), "", "/"),
                      outputFileName,
                      ".rda",
                      sep = ""))
    #:write.table(cbind(Des, Data),
    #:            file = paste(outputFileFolder,
    #:                         ifelse(endsWith(outputFileFolder, "/"), "", "/"),
    #:                         outputFileName,
    #:                         ".txt",
    #:                         sep = ""),
    #:            quote = FALSE,
    #:            sep = "\t",
    #:            na = "",
    #:            col.names = TRUE,
    #:            row.names = FALSE)
    # draw and save boxplot
    #:png(filename = paste(outputFileFolder,
    #:                     ifelse(endsWith(outputFileFolder, "/"), "", "/"),
    #:                     outputFileName,
    #:                     "__boxplot.png",
    #:                     sep = ""),
    #:    width = 30*dim(Data)[2]+300,
    #:    height = 1500,
    #:    units = "px")
    #:par(las = 2,
    #:    mai = c(4.5, 1, 0.5, 1),
    #:    cex = 1.5)
    #:boxplot(Data)
    #:dev.off()
    # collapse the data of same antibody
    sSplit <- "="
    vGeneAb <- apply(Des,
                     1,
                     function(x){
                         paste(x, collapse = sSplit)
                     })
    rownames(Data) <- vGeneAb
    rownames(Des) <- vGeneAb
    vGeneAb1 <- sapply(vGeneAb,
                       function(x){substr(x,
                                          1,
                                          nchar(x) - 4)})
    vGeneAb2 <- names(table(vGeneAb1)[table(vGeneAb1) > 1])
    for (sGeneAb in vGeneAb2) {
        vb <- vGeneAb1 == sGeneAb
        mGeneAb <- Data[vb, ]
        m1 <- ifelse(is.na(mGeneAb),
                     0,
                     1)
        #stopifnot(all(apply(m1, 2, sum) == 1))
        if (all(apply(m1,
                      2,
                      sum) != 1)) {
            print(c(sGeneAb, "count > 1"))
            stop(sGeneAb, "count > 1")
        }
        vRow <- row(m1)[m1==1]
        vCollapsed <- vector(length = length(vRow))
        for (i in seq(length(vRow))) {
            vCollapsed[i] <- mGeneAb[vRow[i], i]
        }
        for (s in vGeneAb[vGeneAb1 == sGeneAb]) {
            Data[s, ] <- vCollapsed
            Des[s, 2] <- paste(strsplit(sGeneAb,
                                        split = sSplit)[[1]][2],
                               "-?-?",
                               sep = "")
        }
    }
    Data <- Data[!duplicated(Des), ] # 1st
    Des  <-  Des[!duplicated(Des), ] # 2nd
    rownames(Data) <- NULL
    rownames(Des) <- NULL
    save(Des, Data,
         file = paste(outputFileFolder,
                      ifelse(endsWith(outputFileFolder, "/"), "", "/"),
                      outputFileName,
                      "__AbCollapsed.rda",
                      sep = ""))
    #:write.table(cbind(Des, Data),
    #:            file = paste(outputFileFolder,
    #:                         ifelse(endsWith(outputFileFolder, "/"), "", "/"),
    #:                         outputFileName,
    #:                         "__AbCollapsed.txt",
    #:                         sep = ""),
    #:            quote = FALSE,
    #:            sep = "\t",
    #:            na = "",
    #:            col.names = TRUE,
    #:            row.names = FALSE)
    # draw and save boxplot
    #:png(filename = paste(outputFileFolder,
    #:                     ifelse(endsWith(outputFileFolder, "/"), "", "/"),
    #:                     outputFileName,
    #:                     "__AbCollapsed__boxplot.png",
    #:                     sep = ""),
    #:    width = 30*dim(Data)[2]+300,
    #:    height = 1500,
    #:    units = "px")
    #:par(las = 2,
    #:    mai = c(4.5, 1, 0.5, 1),
    #:    cex = 1.5)
    #:boxplot(Data)
    #:dev.off()
    options(warn = 0)
    return(list(Data = Data,
                Des = Des))
}

ProcessMethylation27Data <- function(inputFilePath,
                                     outputFileName,
                                     outputFileFolder = ".",
                                     fileSource = "TCGA-Assembler") {
    options(warn = -1)
    if (outputFileFolder != ".") {dir.create(outputFileFolder, recursive = TRUE)}
    if (fileSource == "TCGA-Assembler") {
        # read in data
        #InData <- read.table(file = inputFilePath, header = TRUE, sep = "\t", na.string = "", quote = "", stringsAsFactors = FALSE, check.names = FALSE)
        load(gsub("\\.txt",
                  ".rda",
                  inputFilePath))
        InData <- dPiped
        InData[, 1] <- as.character(InData[, 1])
        Methy27Des <- as.matrix(InData[, 1:4, drop = FALSE])
        colnames(Methy27Des) <- c("REF",
                                  "geneSymbolId",
                                  "ChromosomeID",
                                  "CoordinateID")
        ID <- which(Methy27Des[, "geneSymbolId"] == "DISCONTUNUED")
        Methy27Des[ID, "geneSymbolId"] <- ""
        Methy27Data <- as.matrix(InData[, 5:dim(InData)[2], drop = FALSE])
        mode(Methy27Data) <- "numeric"
    }
    if (fileSource == "Firehose") {
        # read in data
        InData <- read.table(file = inputFilePath,
                             header = TRUE,
                             sep = "\t",
                             na.string = "",
                             quote = "",
                             stringsAsFactors = FALSE,
                             check.names = FALSE)
        InData <- InData[2:dim(InData)[1], , drop = FALSE]
        # remove redudant columns of CpG site descriptions
        geneSymbolId <- as.matrix(InData[, seq(3, dim(InData)[2], 4), drop = FALSE])
        ID <- which(geneSymbolId == "DISCONTUNUED")
        geneSymbolId[ID] <- ""
        geneSymbolId <- t(unique(t(geneSymbolId)))
        ChromID <- as.matrix(InData[, seq(4, dim(InData)[2], 4), drop = FALSE])
        ChromID <- t(unique(t(ChromID)))
        CoordID <- as.matrix(InData[, seq(5, dim(InData)[2], 4), drop = FALSE])
        CoordID <- t(unique(t(CoordID)))
        Methy27Data <- as.matrix(InData[, seq(2, dim(InData)[2], 4), drop = FALSE])
        mode(Methy27Data) <- "numeric"
        REF <- InData[, 1]
        # Check whether the probe information is the same between samples
        if (dim(geneSymbolId)[2] != 1) {
            stop("Some probe's gene symbols are not consistent between samples")
        }
        if (dim(ChromID)[2] != 1) {
            stop("Some probe's chromosome IDs are not consistent between samples")
        }
        if (dim(CoordID)[2] != 1) {
            stop("Some probe's genmoe coordinates are not consistent between samples")
        }
        ChromID <- toupper(ChromID)
        Methy27Des <- cbind(REF,
                            geneSymbolId,
                            ChromID,
                            CoordID)
        colnames(Methy27Des) <- c("REF",
                                  "geneSymbolId",
                                  "ChromosomeID",
                                  "CoordinateID")
    }
    # some probes have more than one gene symbols, collpase the data
    AddGeneDes <- matrix("",
                         dim(Methy27Des)[1]*3,
                         4)
    AddData <- matrix(NA,
                      dim(Methy27Data)[1]*3,
                      dim(Methy27Data)[2])
    AddIndex <- 0
    for (ProbeID in 1:dim(Methy27Des)[1]) {
        GrepID <- grep(";",
                       Methy27Des[ProbeID, "geneSymbolId"])
        if (length(GrepID) > 0) {
            TempStr <- strsplit(Methy27Des[ProbeID, "geneSymbolId"],
                                split = ";")
            Methy27Des[ProbeID, "geneSymbolId"] <- TempStr[[1]][1]
            for (TempStri in 2:length(TempStr[[1]])) {
                TempAddGeneDes <- c(Methy27Des[ProbeID, "REF"],
                                    TempStr[[1]][TempStri],
                                    Methy27Des[ProbeID, "ChromosomeID"],
                                    Methy27Des[ProbeID, "CoordinateID"])
                AddIndex <- AddIndex + 1
                AddGeneDes[AddIndex, ] <- TempAddGeneDes
                AddData[AddIndex, ] <- Methy27Data[ProbeID, ]
            }
        }
    }
    colnames(AddGeneDes) <- colnames(Methy27Des)
    colnames(AddData) <- colnames(Methy27Data)
    Methy27Des <- rbind(Methy27Des,
                        AddGeneDes[1:AddIndex, , drop = FALSE])
    Methy27Data <- rbind(Methy27Data,
                         AddData[1:AddIndex, , drop = FALSE])
    Methy27Des[, "ChromosomeID"] <- gsub(" ",
                                         "",
                                         Methy27Des[, "ChromosomeID"])
    Methy27Des[, "CoordinateID"] <- gsub(" ",
                                         "",
                                         Methy27Des[, "CoordinateID"])
    OrderID <- order(as.numeric(Methy27Des[, "ChromosomeID"]),
                     as.numeric(Methy27Des[, "CoordinateID"]),
                     Methy27Des[, "geneSymbolId"],
                     na.last = TRUE,
                     decreasing = FALSE)
    Methy27Des <- Methy27Des[OrderID, , drop = FALSE]
    Methy27Data <- Methy27Data[OrderID, , drop = FALSE]
    Data <- Methy27Data
    Des <- Methy27Des
    #lResult <- CombineRedundantFeature(Data = Data, Des = Des)
    lResult <- list(Data = Data,
                    Des = Des)
    Data <- lResult$Data
    Des <- lResult$Des
    # Draw box plot
    #:png(filename = paste(outputFileFolder,
    #:                     ifelse(endsWith(outputFileFolder, "/"), "", "/"),
    #:                     outputFileName,
    #:                     "__boxplot.png",
    #:                     sep = ""),
    #:    width = 30*dim(Data)[2]+300,
    #:    height = 1500,
    #:    units = "px")
    #:par(las = 2,
    #:    mai = c(4.5, 1, 0.5, 1),
    #:    cex = 1.5)
    #:boxplot(Data)
    #:dev.off()
    # save output data files
    rownames(Data) <- NULL
    rownames(Des) <- NULL
    save(Des, Data,
         file = paste(outputFileFolder,
                      ifelse(endsWith(outputFileFolder, "/"), "", "/"),
                      outputFileName,
                      ".rda",
                      sep = ""))
    #:write.table(cbind(Des,
    #:                  Data),
    #:            file = paste(outputFileFolder,
    #:                         ifelse(endsWith(outputFileFolder, "/"), "", "/"),
    #:                         outputFileName,
    #:                         ".txt",
    #:                         sep = ""),
    #:            quote = FALSE,
    #:            sep = "\t",
    #:            na = "",
    #:            col.names = TRUE,
    #:            row.names = FALSE)
    options(warn = 0)
    return(list(Data = Data,
                Des = Des))
}

ProcessSomaticMutationData <- function(inputFilePath,
                                       outputFileName,
                                       outputFileFolder = ".") {
    options(warn = -1)
    if (outputFileFolder != ".") {dir.create(outputFileFolder, recursive = TRUE)}
    # Read in data
    #InData <- read.table(file = inputFilePath, header = TRUE, sep = "\t", na.string = "", quote = "", stringsAsFactors = FALSE, check.names = FALSE)
    load(gsub("\\.txt",
              ".rda",
              inputFilePath))
    InData <- dPiped
    InData[, 1] <- as.character(InData[, 1])
    InputData <- InData
    InputData[, "Tumor_Sample_Barcode"] <- paste(InputData[, "Tumor_Sample_Barcode"],
                                                 InputData[, "Matched_Norm_Sample_Barcode"],
                                                 sep = "/")
    colnames(InputData)[16] <- "Sample_Pair_Barcode"
    InputData <- InputData[InputData[, "Mutation_Status"] == "Somatic",
                           c("Hugo_Symbol",
                             "Entrez_Gene_Id",
                             "NCBI_Build",
                             "Chromosome",
                             "Start_Position",
                             "End_Position",
                             "Strand",
                             "Variant_Classification",
                             "Variant_Type",
                             "Reference_Allele",
                             "Tumor_Seq_Allele1",
                             "Tumor_Seq_Allele2",
                             "dbSNP_RS",
                             "dbSNP_Val_Status",
                             "Sample_Pair_Barcode",
                             "Match_Norm_Seq_Allele1",
                             "Match_Norm_Seq_Allele2",
                             "Tumor_Validation_Allele1",
                             "Tumor_Validation_Allele2",
                             "Match_Norm_Validation_Allele1",
                             "Match_Norm_Validation_Allele2",
                             "Verification_Status",
                             "Validation_Status")]
    samplePair <- unique(InputData[, "Sample_Pair_Barcode"])
    TumorAllele1 <- "Tumor_Seq_Allele1"
    TumorAllele2 <- "Tumor_Seq_Allele2"
    MatchNormAllele1 <- "Match_Norm_Seq_Allele1"
    MatchNormAllele2 <- "Match_Norm_Seq_Allele2"
    Input <- InputData[, c("Hugo_Symbol",
                           "NCBI_Build",
                           "Chromosome",
                           "Start_Position",
                           "End_Position",
                           "Strand",
                           "Variant_Classification",
                           "Variant_Type",
                           "Reference_Allele",
                           TumorAllele1,
                           TumorAllele2,
                           MatchNormAllele1,
                           MatchNormAllele2)]
    IDLarge <- which(Input[, TumorAllele1] > Input[, TumorAllele2])
    Y <- Input[IDLarge, TumorAllele1]
    Input[IDLarge, TumorAllele1] <- Input[IDLarge, TumorAllele2]
    Input[IDLarge, TumorAllele2] <- Y
    rm(IDLarge)
    IDLarge <- which(Input[, MatchNormAllele1] > Input[, MatchNormAllele2])
    Y <- Input[IDLarge, MatchNormAllele1]
    Input[IDLarge, MatchNormAllele1] <- Input[IDLarge, MatchNormAllele2]
    Input[IDLarge, MatchNormAllele2] <- Y
    rm(IDLarge)
    A <- apply(Input,
               1,
               paste,
               collapse = "|||")
    DupID <- duplicated(A)
    UniqueA <- A[sort(which(DupID == FALSE))]
    Feature <- InputData[sort(which(DupID == FALSE)),
                         -which(colnames(InputData) == "Sample_Pair_Barcode")]
    FeatureMatrix <- matrix(0,
                            length(UniqueA),
                            length(samplePair))
    colnames(FeatureMatrix) <- samplePair
    for (i in 1:length(A)) {
        RowIDi <- which(UniqueA == A[i])
        ColIDi <- which(samplePair == InputData[i, "Sample_Pair_Barcode"])
        FeatureMatrix[RowIDi, ColIDi] <- FeatureMatrix[RowIDi, ColIDi] + 1
    }
    colnames(Feature)[1] <- "geneSymbolId"
    Des <- Feature
    Des <- cbind(Des[, 1, drop = FALSE],
                 Description = apply(Des,
                                     1,
                                     function(x){
                                         paste("chr",
                                               x[4],
                                               ":",
                                               as.numeric(x[5]),
                                               "-",
                                               as.numeric(x[6]),
                                               ",",
                                               x[7],
                                               ",",
                                               x[8],
                                               ",",
                                               x[9],
                                               ",",
                                               x[10],
                                               x[10],
                                               ">",
                                               x[11],
                                               x[12],
                                               sep = "")}),
                 Des[, seq(2, ncol(Des))])
    Data <- FeatureMatrix
    #lResult <- CombineRedundantFeature(Data = Data, Des = Des)
    lResult <- list(Data = Data,
                    Des = Des)
    Data <- lResult$Data
    Des <- lResult$Des
    # save output data files
    rownames(Data) <- NULL
    rownames(Des) <- NULL
    save(Des, Data,
         file = paste(outputFileFolder,
                      ifelse(endsWith(outputFileFolder, "/"), "", "/"),
                      outputFileName,
                      "_mutationLevel.rda",
                      sep = ""))
    #:write.table(cbind(Des,
    #:                  Data),
    #:            file = paste(outputFileFolder,
    #:                         ifelse(endsWith(outputFileFolder, "/"), "", "/"),
    #:                         outputFileName,
    #:                         "_mutationLevel.txt",
    #:                         sep = ""),
    #:            quote = FALSE,
    #:            sep = "\t",
    #:            na = "",
    #:            col.names = TRUE,
    #:            row.names = FALSE)
    # draw and save boxplot
    ## png(filename = paste(outputFileFolder,
    #                     ifelse(endsWith(outputFileFolder, "/"), "", "/"),
    # outputFileName, "__boxplot.png", sep = ""), width = 30*dim(Data)[2]+300, height = 1500, units = "px")
    ## par(las = 2, mai = c(4.5, 1, 0.5, 1), cex = 1.5)
    ## boxplot(Data)
    ## dev.off()
    # output data files geneLevel (each row is a gene)
    m_Data_geneLevel <- sapply(apply(Data,
                                     2,
                                     function(x){
                                         split(x, Des[, "geneSymbolId"])
                                     }),
                               function(x){
                                   sapply(x, sum)
                               })
    v_variant <- sapply(split(Des$Variant_Classification,
                              Des[, "geneSymbolId"]),
                        function(x){
                            paste(sort(unique(x)),
                                  collapse = '|')
                        })
    Data <- m_Data_geneLevel
    Des <- cbind('geneSymbolId' = rownames(m_Data_geneLevel),
                 'Variant_Classification' = v_variant)
    save(Des, Data,
         file = paste(outputFileFolder,
                      ifelse(endsWith(outputFileFolder, "/"), "", "/"),
                      outputFileName,
                      "_geneLevel.rda",
                      sep = ""))
    #:write.table(cbind(Des,
    #:                  Data),
    #:            file = paste(outputFileFolder,
    #:                         ifelse(endsWith(outputFileFolder, "/"), "", "/"),
    #:                         outputFileName,
    #:                         "_geneLevel.txt",
    #:                         sep = ""),
    #:            quote = FALSE,
    #:            sep = "\t",
    #:            na = "",
    #:            col.names = T,
    #:            row.names = F)
    # draw and save boxplot
    ## png(filename = paste(outputFileFolder,
    #                     ifelse(endsWith(outputFileFolder, "/"), "", "/"),
    # outputFileName, "__boxplot_geneLevel.png", sep = ""), width = 30*dim(Data)[2]+300, height = 1500, units = "px")
    ## par(las = 2, mai = c(4.5, 1, 0.5, 1), cex = 1.5)
    ## boxplot(Data)
    ## dev.off()
    options(warn = 0)
    return(list(Data = Data,
                Des = Des))
}

ProcessCPTACData <- function(inputFilePath,
                             outputFileName,
                             outputFileFolder = ".") {
    options(warn = -1)
    if (outputFileFolder != ".") {dir.create(outputFileFolder, recursive = TRUE)}
    out <- inputFilePath
    d <- read.csv(out,
                  skip = 0,
                  sep = "\t",
                  as.is = T,
                  na.strings = "",
                  check.names = FALSE)
    colnames(d)[1] <- "geneSymbolId"
    vColNote <- c("geneSymbolId",
                  "Description",
                  "Organism",
                  "Chromosome",
                  "Locus")
    mInfo <- as.matrix(d[, vColNote])
    mData <- as.matrix(d[, setdiff(colnames(d), vColNote)])
    #colnames(mData) <- gsub("\\.", "-", colnames(mData))
    colnames(mData) <- gsub(" Spectral Counts",
                            " Log Ratio",
                            colnames(mData))
    for (n in grep("^OVARIAN-CONTROL:", colnames(mData))) {
        substr(colnames(mData)[n], 16, 16) <- ':'
    }
    for (n in grep("-[0-9][0-9]* (Log|Unshared|Spectral)", colnames(mData))) {
        substr(colnames(mData)[n], 17, 17) <- ':'
    }
    vnUnsh <- grep(" Unshared Log Ratio", colnames(mData))
    #vnShUn <- grep(" Log Ratio", colnames(mData))
    vnShUn <- grep("[^d] Log Ratio", colnames(mData))
    mDataShUn <- mData[, vnShUn]  # barcode:Log-Ratio
    mDataUnsh <- mData[, vnUnsh]  # barcode:Unshared-Log-Ratio
    colnames(mDataShUn) <- unlist(strsplit(colnames(mDataShUn),
                                           split = " Log Ratio"))
    colnames(mDataUnsh) <- unlist(strsplit(colnames(mDataUnsh),
                                           split = " Unshared Log Ratio"))
    mDataShUn <- mDataShUn[, colnames(mDataShUn)]
    mDataUnsh <- mDataUnsh[, colnames(mDataUnsh)]
    stopifnot(all(colnames(mDataUnsh) == colnames(mDataShUn)))
    Des <- mInfo
    mode(mDataShUn) <- "numeric"
    #lResult <- CombineRedundantFeature(Data = mDataShUn, Des = Des)
    lResult <- list(Data = mDataShUn, Des = Des)
    Data <- lResult$Data
    Des <- lResult$Des
    # save output data files
    rownames(Data) <- NULL
    rownames(Des) <- NULL
    save(Des, Data,
         file = paste(outputFileFolder,
                      ifelse(endsWith(outputFileFolder, "/"), "", "/"),
                      outputFileName,
                      "_allPeptides.rda",
                      sep = ""))
    #:write.table(cbind(Des,
    #:                  Data),
    #:            file = paste(outputFileFolder,
    #:                         ifelse(endsWith(outputFileFolder, "/"), "", "/"),
    #:                         outputFileName,
    #:                         "_allPeptides.txt",
    #:                         sep = ""),
    #:            quote = FALSE,
    #:            sep = "\t",
    #:            na = "",
    #:            col.names = TRUE,
    #:            row.names = FALSE)
    # draw and save boxplot
    #:png(filename = paste(outputFileFolder,
    #:                     ifelse(endsWith(outputFileFolder, "/"), "", "/"),
    #:                     outputFileName,
    #:                     "_allPeptides__boxplot.png",
    #:                     sep = ""),
    #:    width = 30*dim(Data)[2]+300,
    #:    height = 1500,
    #:    units = "px")
    #:par(las = 2,
    #:    mai = c(4.5, 1, 0.5, 1),
    #:    cex = 1.5)
    #:boxplot(Data)
    #:dev.off()
    lShUn <- list(Data = Data,
                  Des = Des)
    Des <- mInfo
    mode(mDataUnsh) <- "numeric"
    #lResult <- CombineRedundantFeature(Data = mDataUnsh, Des = Des)
    lResult <- list(Data = mDataUnsh,
                    Des = Des)
    Data <- lResult$Data
    Des <- lResult$Des
    # save output data files
    rownames(Data) <- NULL
    rownames(Des) <- NULL
    save(Des, Data,
         file = paste(outputFileFolder,
                      ifelse(endsWith(outputFileFolder, "/"), "", "/"),
                      outputFileName,
                      "_unsharedPeptides.rda",
                      sep = ""))
    #:write.table(cbind(Des,
    #:                  Data),
    #:            file = paste(outputFileFolder,
    #:                     ifelse(endsWith(outputFileFolder, "/"), "", "/"),
    #:                         outputFileName,
    #:                         "_unsharedPeptides.txt",
    #:                         sep = ""),
    #:            quote = FALSE,
    #:            sep = "\t",
    #:            na = "",
    #:            col.names = TRUE,
    #:            row.names = FALSE)
    # draw and save boxplot
    #:png(filename = paste(outputFileFolder,
    #:                     ifelse(endsWith(outputFileFolder, "/"), "", "/"),
    #:                     outputFileName,
    #:                     "_unsharedPeptides__boxplot.png",
    #:                     sep = ""),
    #:    width = 30*dim(Data)[2]+300,
    #:    height = 1500,
    #:    units = "px")
    #:par(las = 2,
    #:    mai = c(4.5, 1, 0.5, 1),
    #:    cex = 1.5)
    #:boxplot(Data)
    #:dev.off()
    lUnsh <- list(Data = Data,
                  Des = Des)
    options(warn = 0)
    return(list(allPeptides = lShUn,
                unsharedPeptides = lUnsh))
}


#  =============================================================================
#  Auxiliary Functions of Module B
#  =============================================================================
ToPatientData <- function(TCGAData) {
    TCGABarcode <- colnames(TCGAData)
    TCGABarcode <- sapply(strsplit(TCGABarcode,
                                   split = "-"),
                          function(x){
                              return(paste(x[1],
                                           x[2],
                                           x[3],
                                           substr(x[4],
                                                  1,
                                                  2),
                                           sep = "-"))
                          })
    DuplicatedLabel <- duplicated(TCGABarcode)
    ID <- which(DuplicatedLabel == FALSE)
    TCGAData <- TCGAData[, ID, drop = FALSE]
    TCGABarcode <- TCGABarcode[ID]
    colnames(TCGAData) <- TCGABarcode
    TCGAData
}

CombineRedundantFeature <- function(Data, Des) {
    RowName <- Des[, 1]
    if (dim(Des)[2] > 1) {
        for (i in 2:dim(Des)[2]) {
            RowName <- paste(RowName,
                             Des[, i],
                             sep = "||")
        }
    }
    UniqueRowName <- unique(RowName)
    if (length(RowName) == length(UniqueRowName)) {
        lResult <- list(Data = Data,
                        Des = Des)
    }
    else {
        dupRowName <- unique(RowName[which(duplicated(RowName))])
        IDKeep <- which(!duplicated(RowName))
        for (i in 1:length(dupRowName)) {
            IDi <- which(RowName == dupRowName[i])
            IDKeepi <- intersect(IDKeep, IDi)
            Data[IDKeepi, ] <- apply(Data[IDi, , drop = FALSE],
                                     2,
                                     mean,
                                     na.rm = TRUE)
        }
        IDKeep <- sort(IDKeep)
        lResult <- list(Data = Data[IDKeep, , drop = FALSE],
                        Des = Des[IDKeep, , drop = FALSE])
    }
    lResult
}

Normalize.quantiles <- function(Data) {
    VectorData <- as.vector(Data)
    VectorData <- VectorData[!is.na(VectorData)]
    VectorData <- sort(VectorData,
                       decreasing = FALSE)
    NumVectorData <- length(VectorData)
    NormData <- matrix(NA,
                       dim(Data)[1],
                       dim(Data)[2])
    for (i in 1:dim(Data)[2]) {
        if ((i %% 10 == 0) || (i == dim(Data)[2])) {
            writeLines(paste("Normalizing data, ",
                             round(i/dim(Data)[2]*100,
                                   digits = 2),
                             "% done.",
                             sep = ""))
        }
        IDNonNA <- which(!is.na(Data[, i]))
        RankNonNA <- rank(Data[IDNonNA, i],
                          ties.method = "average")
        NumNonNA <- length(IDNonNA)
        Sep <- round(seq(1,
                         NumVectorData,
                         (NumVectorData-1)/NumNonNA))
        NormData[IDNonNA, i] <- sapply(RankNonNA,
                                       function(x) {
                                           NumTie = sum(RankNonNA == x)
                                           mean(VectorData[(Sep[ceiling(x-NumTie/2)]+1):Sep[floor(x+NumTie/2)+1]])
                                       })
    }
    rownames(NormData) <- rownames(Data)
    colnames(NormData) <- colnames(Data)
    NormData
}

GetMethylation450Data <- function(inputFilePath,
                                  outputFileName,
                                  outputFileFolder = ".") {
    OriginalOutputFileFolder <- outputFileFolder
    outputFileFolder <- paste(outputFileFolder,
                              ifelse(endsWith(outputFileFolder, "/"), "", "/"),
                              outputFileName,
                              ".TempFiles/",
                              sep = "")
    if (outputFileFolder != ".") {dir.create(outputFileFolder, recursive = TRUE)}
    FileHandle <- file(inputFilePath, "rt")
    on.exit(close(FileHandle))
    ReadIn <- readLines(FileHandle, 1)
    ColumnName1 <- strsplit(ReadIn, split = "\t")[[1]]
    ReadIn <- readLines(FileHandle, 1)
    ColumnName2 <- strsplit(ReadIn, split = "\t")[[1]]
    LineNum <- 10000
    Iter <- 0
    Flag <- TRUE
    while (Flag) {
        ReadIn <- readLines(FileHandle, LineNum)
        if (length(ReadIn) > 0) {
            Iter <- Iter + 1
            if (Iter %% 5 == 0) {
                writeLines(paste("Load ",
                                 floor(Iter/49*100),
                                 "% of the humanmethylation450 data set.",
                                 sep = ""))
            }
            if (length(ReadIn) == 1) {
                TempMatrix <- rbind(strsplit(ReadIn,
                                             split = "\t")[[1]])
                colnames(TempMatrix) <- ColumnName1
                TempMethyData <- rbind(TempMatrix[, seq(2, dim(TempMatrix)[2], 4), drop = FALSE])
                geneSymbolId <- unique(TempMatrix[, seq(3, dim(TempMatrix)[2], 4)])
                ChromID <- unique(TempMatrix[, seq(4, dim(TempMatrix)[2], 4)])
                CoordID <- unique(TempMatrix[, seq(5, dim(TempMatrix)[2], 4)])
                REF <- unique(TempMatrix[, 1])
                ChromID <- toupper(ChromID)
                # Check whether the probe information is the same between samples
                if (length(geneSymbolId) != 1) {
                    stop("Some probe's gene symbols are not consistent between samples")
                }
                if (length(ChromID) != 1) {
                    stop("Some probe's chromosome IDs are not consistent between samples")
                }
                if (length(CoordID) != 1) {
                    stop("Some probe's genmoe coordinates are not consistent between samples")
                }
                TempDes <- cbind(REF,
                                 geneSymbolId,
                                 ChromID,
                                 CoordID)
            }
            else {
                ReadInSep <- strsplit(ReadIn,
                                      split = "\t")
                TempMatrix <- matrix(NA,
                                     length(ReadIn),
                                     length(ColumnName1))
                for (i in 1:length(ReadInSep)) {
                    TempMatrix[i, ] <- ReadInSep[[i]]
                }
                colnames(TempMatrix) <- ColumnName1
                geneSymbolId <- TempMatrix[, seq(3, dim(TempMatrix)[2], 4), drop = FALSE]
                ChromID <- TempMatrix[, seq(4, dim(TempMatrix)[2], 4), drop = FALSE]
                CoordID <- TempMatrix[, seq(5, dim(TempMatrix)[2], 4), drop = FALSE]
                TempMethyData <- as.matrix(TempMatrix[, seq(2, dim(TempMatrix)[2], 4), drop = FALSE])
                geneSymbolId <- t(unique(t(geneSymbolId)))
                ChromID <- t(unique(t(ChromID)))
                CoordID <- t(unique(t(CoordID)))
                REF <- TempMatrix[, 1]
                ChromID <- toupper(ChromID)
                # Check whether the probe information is the same between samples
                if (dim(geneSymbolId)[2] != 1) {
                    stop("Some probe's gene symbols are not consistent between samples")
                }
                if (dim(ChromID)[2] != 1) {
                    stop("Some probe's chromosome IDs are not consistent between samples")
                }
                if (dim(CoordID)[2] != 1) {
                    stop("Some probe's genmoe coordinates are not consistent between samples")
                }
                TempDes <- cbind(REF,
                                 geneSymbolId,
                                 ChromID,
                                 CoordID)
            }
            colnames(TempDes) <- c("REF",
                                   "geneSymbolId",
                                   "ChromosomeID",
                                   "CoordinateID")
            mode(TempMethyData) <- "numeric"
            TempMethyData <- round(x = TempMethyData,
                                   digits = 4)
            save(TempMethyData,
                 TempDes,
                 file = paste(outputFileFolder,
                              ifelse(endsWith(outputFileFolder, "/"), "", "/"),
                              "Methy450Data",
                              Iter,
                              ".rda",
                              sep = ""))
        }
        if (length(ReadIn) < LineNum) {
            Flag <- FALSE
        }
    }
    writeLines(paste("Load 100% of the humanmethylation450 data set.",
                     sep = ""))
    rm(ReadIn,
       ReadInSep,
       TempMatrix,
       geneSymbolId,
       ChromID,
       CoordID,
       TempMethyData,
       REF,
       TempDes)
    load(paste(outputFileFolder,
               ifelse(endsWith(outputFileFolder, "/"), "", "/"),
               "Methy450Data",
               1,
               ".rda",
               sep = ""))
    Methy450Des <- TempDes
    Methy450Data <- TempMethyData
    for (i in 2:Iter) {
        load(paste(outputFileFolder,
                   ifelse(endsWith(outputFileFolder, "/"), "", "/"),
                   "Methy450Data",
                   i,
                   ".rda",
                   sep = ""))
        Methy450Des  <- rbind(Methy450Des , TempDes)
        Methy450Data <- rbind(Methy450Data, TempMethyData)
    }
    save(Methy450Data,
         Methy450Des,
         file = paste(outputFileFolder,
                      ifelse(endsWith(outputFileFolder, "/"), "", "/"),
                      outputFileName,
                      ".TempData.rda",
                      sep = ""))
}


#  =============================================================================
#  Check whether this is the most updated version of TCGA-Assembler
#  =============================================================================

#' CheckVersion check the newset version of TCGA-Assembler
#'
#' @return a message if not the newest version
CheckVersion <- function() {
    library(httr)
    library(stringr)
    s1Content <- try(content(GET("http://www.compgenome.org/TCGA-Assembler/"),
                             as = "text"),
                     silent = TRUE)
    if (class(s1Content) == "try-error") {
        rm(s1Content)
    } else {
        lmiStartId <- str_locate_all(string = s1Content,
                                     pattern = "CheckVersionNumber1")
        if (dim(lmiStartId[[1]])[1] == 0) {
            rm(lmiStartId, s1Content)
        } else {
            lmiStartId <- lmiStartId[[1]][1, "end"]+1
            s1Content <- substr(s1Content, lmiStartId, nchar(s1Content))
            lmiStartId <- str_locate_all(string = s1Content,
                                         pattern = "\">")[[1]][1, "end"]+1
            lmiEndId <- str_locate_all(string = s1Content,
                                       pattern = "</span>")[[1]][1, "start"]-1
            sVersion <- substr(s1Content, lmiStartId, lmiEndId)
            if (sVersion != "2.1.0") {
                writeLines("\n")
                writeLines("***************************************************************")
                writeLines("A new version of TCGA-Assembler is available!")
                writeLines(paste("Please download version ", sVersion,
                                 " at www.compgenome.org/TCGA-Assembler/",
                                 sep = ""))
                writeLines("***************************************************************")
                writeLines("\n")
            }
            rm(lmiStartId, s1Content, lmiEndId, sVersion)
        }
    }
}
CheckVersion()

#  end

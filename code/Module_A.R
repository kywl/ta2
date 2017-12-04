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
#  TCGA-Assembler Version 2 Module A
#  =============================================================================


#  =============================================================================
#  variable prefix example
#  =============================================================================
#  s : string, character
#  s1: character vector with only 1 element
#  b : bool
#  n : numeric value
#  v : vector
#  d : data.frame
#  m : matrix
#  l : list
#  ld : list of data.frame
#  lm : list of matrix
#  ln : list of numeric value
#  ls : list of character
#  lv : list of vector


#  =============================================================================
#  internal functions, NOT used directly by user
#  =============================================================================

#' Get fields names of GDC entities
#'
#' @param sArchive String of archive type: "legacy" or "".
#' @param sEndpoint String of endpoint: "files".
#' @return Character vector of all available field names.
#' @examples
#' sFieldList <- FieldsList("legacy")
#' sFieldList <- FieldsList("")
FieldsList <- function(sArchive = "legacy",
                       sEndpoint = "files") {
    library(rjson)
    sUrl <- paste("https://gdc-api.nci.nih.gov/",
                  sArchive,
                  ifelse(sArchive == "", "", "/"),
                  sEndpoint,
                  "/_mapping",
                  sep = "")
    sOption <- "--silent --show-error"
    sArgument <- paste(sOption, sUrl)
    sJson <- paste(system2("curl", sArgument, stdout = TRUE), collapse = "")
    return(fromJSON(sJson)$fields)
}

#' Define the default fields used in metadata file
#'
#' @return Character vector of selected fields names for metadata.
FieldsMeta <- function() {
    return(c(# "access",
             "archive.file_name",
             # "cases.case_id",
             "cases.project.project_id",
             "cases.samples.portions.analytes.aliquots.submitter_id",  # 1st
             # "cases.samples.portions.analytes.submitter_id",
             "cases.samples.portions.submitter_id",  # 2nd
             "cases.samples.sample_type",
             "cases.samples.sample_type_id",  #  sprintf("%02d", c(seq(14), 20, 40, 50, 60, 61))
             "data_category",
             "data_type",
             "experimental_strategy",
             "file_id",
             "file_name",
             "file_size",
             # "md5sum",
             "platform",
             "updated_datetime"))
}

#' Get count of GDC entities in /archive/endpoint
#'
#' @param sArchive String of archive type: "legacy" or "".
#' @param sEndpoint String of endpoint: "files".
#' @return Count of entities in specified archive & endpoint.
#' @examples
#' nEntityCount <- EntityCount("legacy")
#' nEntityCount <- EntityCount("")
EntityCount <- function(sArchive = "legacy",
                        sEndpoint = "files") {
    library(rjson)
    sUrl <- paste("https://gdc-api.nci.nih.gov/",
                  sArchive,
                  ifelse(sArchive == "", "", "/"),
                  sEndpoint,
                  sep = "")
    sOption <- "--silent --show-error"
    sArgument <- paste(sOption, sUrl)
    sJson <- paste(system2("curl", sArgument, stdout = TRUE), collapse = "")
    # stopifnot(length(fromJSON(sJson)$warnings) == 0)
    return(fromJSON(sJson)$data$pagination$total)
}

#' Get a string of current time (YYYYMMDDhhmmss)
#'
#' @return String of current time (YYYYMMDDhhmmss).
#' @examples
#' TimeNow()
TimeNow <- function() {
    return(gsub("[- :]", "", as.character(Sys.time())))
}

#' Get index of the entity with newest archive file version
#'
#' @param sArchiveNames Character vector of "archive.file_name". All of
#' these entities have same "file_name".
#' @return Index of the newest one.
#' @examples
#' v <- c("mdanderson.org_BRCA.MDA_RPPA_Core.Level_3.114.1.0.tar.gz",
#'        "mdanderson.org_BRCA.MDA_RPPA_Core.Level_3.2.1.0.tar.gz",
#'        "mdanderson.org_BRCA.MDA_RPPA_Core.Level_3.10.1.0.tar.gz")
#' n <- ArchiveNewest(v)
ArchiveNewest <- function(sArchiveNames) {
    mVersion <- sapply(strsplit(sArchiveNames, split = "\\."),
                       function(x){rev(x)[5 : 3]})
    mode(mVersion) <- "numeric"  # OR class(mVersion) <- "numeric"
    nOrder <- do.call(order,
                      lapply(seq(nrow(mVersion)),
                             function(x){mVersion[x, ]}))
    return(nOrder[length(nOrder)])
}

#' Get vector of bool indicate the one with newest archive file version or not
#' in one group
#'
#' @param sArchiveNames Character vector of "archive.file_name". All of
#' these entities have same "file_name".
#' @return Bool vector to indicate the newest one in one group
ArchiveNewestInGroup <- function(sArchiveNames) {
    bNewest <- seq(length(sArchiveNames)) == ArchiveNewest(sArchiveNames)
    return(ifelse(bNewest, TRUE, FALSE))
}

#' Get vector of bool indicate the one with newest archive file version or not
#' in every group
#'
#' @param sArchiveNames Character vector of "archive.file_name", could be
#' split into groups, each group with a unique "file_name".
#' @param sSplitFactor Character vector of "file_name" (sorted).
#' @return Bool vector to indicate the newest one in every group.
ArchiveNewestInGroups <- function(sArchiveNames,
                                  sSplitFactor) {
    # sSplitFactor should be sorted before split
    stopifnot(all(sSplitFactor == sSplitFactor[order(sSplitFactor)]))
    return(unlist(lapply(split(sArchiveNames, sSplitFactor),
                         ArchiveNewestInGroup)))
}

#' Choose columns from a file
#'
#' @param sFilename String of filename.
#' @param sFileId String of sFileId with same length of fileName.
#' @param sColNames Character vector of colname, specify the columns choosed.
#' @param sSortByCol Name one column which acts as the rownames.
#' @param sSkipLines Number of lines skipped in \code{read.csv}.
#' @param sNa String indicating the \code{NA} in \code{read.csv}.
#' @return A \code{data.frame} with specified columns from one file.
ColumnsFromFile <- function(sFilename,
                            sFileId,
                            sColNames,
                            sSortByCol = "",
                            sSkipLines = 0,
                            sNa) {
    d <- read.csv(sFilename,
                  sep = "\t",
                  row.names = NULL,
                  as.is = TRUE,
                  skip = sSkipLines,
                  na.strings = sNa,
                  check.names = FALSE)
    if (sSortByCol != "") {
        # probe duplicated, check same value for each probe group
        if (length(d[, sSortByCol]) != length(unique(d[, sSortByCol]))) {
            stopifnot(all(tapply(Reduce(paste, d[, sColNames]), d[, sSortByCol],
                                 function(x){length(unique(x)) == 1})))
            d <- d[!duplicated(d[, sSortByCol]), ]
        }
        stopifnot(length(d[, sSortByCol]) == length(unique(d[, sSortByCol])))
        rownames(d) <- d[, sSortByCol]
        d <- d[order(d[, sSortByCol]), sColNames]
    } else {
        d <- cbind("sFileId" = rep(sFileId, dim(d)[1]), d[, sColNames])
    }
    return(d)
}

#' Choose columns from each file of filenames
#'
#' @param sFilename String of filename, named with file_id.
#' @param sColNames Character vector of colname, specify the columns choosed.
#' @param sSortByCol Name one column which acts as the rownames.
#' @param sSkipLines Number of lines skipped in \code{read.csv}.
#' @param sNa String indicating the \code{NA} in \code{read.csv}.
#' @return List of \code{data.frame} with specified columns from each file.
ColumnsFromFiles <- function(sFileNameById,
                             sColNames,
                             sSortByCol = "",
                             sSkipLines = 0,
                             sNa = "NA") {
    l <- lapply(names(sFileNameById),
                function(s){
                    ColumnsFromFile(sFileNameById[s],
                                    s,
                                    sColNames,
                                    sSortByCol,
                                    sSkipLines,
                                    sNa)
                })
    names(l) <- names(sFileNameById)
    return(l)
}

#' Strip characters at left or right end of each string in a vector
#'
#' @param sUnstripped Character vector of unstripped strings.
#' @param nStrip Number of characters need to be stripped.
#' @param sStripEnd Sting indicating the terminal: "right" or "left".
#' @return Character vector of stripped strings.
StripEnd <- function(sUnstripped,
                     nStrip,
                     sStripEnd = "right") {
    if (nStrip == 0) {
        sStripped <- sUnstripped
    } else  {
        if (sStripEnd %in% c("r", "right")) {
            sStripped <- sapply(strsplit(sUnstripped, split = ""),
                                function(x){
                                    paste(x[-((length(x) - nStrip + 1) : length(x))],
                                          collapse = "")
                                })
        } else if (sStripEnd %in% c("l", "left")) {
            sStripped <- sapply(strsplit(sUnstripped, split = ""),
                                function(x){
                                    paste(x[-(1 : nStrip)], collapse = "")
                                })
        } else {
            print("sStripEnd should be one of 'right', 'r', 'l' or 'left'")
        }
    }
    return(sStripped)
}

#' Make a vector of values named with probes
#'
#' @param dProbeValue A \code{data.frame}, usually read from filename.
#' @param sColProbe String of column name of probe, used as name.
#' @param sColValue String of column name of value.
#' @param nStrip Number of characters need to be stripped from "probe".
#' @param sStripEnd Sting indicating the terminal: "right" or "left".
#' @return Named vector.
ProbeValue <- function(dProbeValue,
                       sColProbe,
                       sColValue,
                       nStrip = 0,
                       sStripEnd = "right") {
    sValue <- dProbeValue[, sColValue]
    names(sValue) <- StripEnd(dProbeValue[, sColProbe], nStrip, sStripEnd)
    return(sValue)
}

#' Cut the "*\\.1\\.*" columns and rows with \code{NA} only
#'
#' @param dMetaData A \code{data.frame} read from metadata file.
#' @param sColPattern String of regular expression pattern to filter colname.
#' @return A \code{data.frame} after cutting.
MetaCut <- function(dMetaData,
                    sColPattern = "\\.1\\.") {
    nColIdx <- grep(sColPattern, colnames(dMetaData))
    nRowIdx <- sapply(seq(dim(dMetaData)[1]),
                      function(x){
                          ifelse(all(is.na(dMetaData[x, nColIdx])), TRUE, FALSE)
                      })
    dMetaCut <- dMetaData[nRowIdx, seq(dim(dMetaData)[2])[-nColIdx]]
    return(dMetaCut)
}

#' Define the filters with specified assay platform
#'
#' @param s1Assay String of assay platform.
#' @return List of filter.
#' @example
#' lFilter <- Filter("methylation_450")
Filter <- function(s1Assay) {
    lFilter <- list()
    sAssay <- c(# Copy number segmentation
                "cna_cnv.hg18",
                "cna_cnv.hg19",
                "cna_nocnv.hg18",
                "cna_nocnv.hg19",
                # Exon junction quantification
                "exonJunction_RNAseq",
                "exonJunction_TotalRNAseq",
                # Exon quantification
                "exon_RNAseq",
                "exon_TotalRNAseq",
                # Gene expression quantification
                "gene_Array",
                "gene.normalized_RNAseq",
                "gene_RNAseq",
                "gene.normalized_TotalRNAseq",
                "gene_TotalRNAseq",
                # Isoform expression quantification
                "isoform.normalized_RNAseq",
                "isoform_RNAseq",
                "isoform.normalized_TotalRNAseq",
                "isoform_TotalRNAseq",
                # Methylation beta value
                "methylation_27",
                "methylation_450",
                # miRNA gene quantification
                "mir_GA.hg18",
                "mir_GA.hg19",
                "mir_GA.hg19.mirbase20",
                # miRNA gene quantification
                "mir_HiSeq.hg18",
                "mir_HiSeq.hg19",
                "mir_HiSeq.hg19.mirbase20",
                # miRNA isoform quantification
                "mirIsoform_GA.hg18",
                "mirIsoform_GA.hg19",
                "mirIsoform_GA.hg19.mirbase20",
                # miRNA isoform quantification
                "mirIsoform_HiSeq.hg18",
                "mirIsoform_HiSeq.hg19",
                "mirIsoform_HiSeq.hg19.mirbase20",
                # Protein expression quantification
                "protein_RPPA",
                # Simple somatic mutation
                "somaticMutation_DNAseq",
                # CPTAC
                "glycoproteome_iTRAQ",
                "phosphoproteome_iTRAQ",
                "proteome_iTRAQ")
    if (!s1Assay %in% sAssay) {  # assayPlatform = s1Assay
        print(paste(c("assayPlatform should be one of:", sAssay), collapse = " "))
    } else if (s1Assay == "cna_cnv.hg18") {
        lFilter$data_category         <- c("Copy number variation")
        lFilter$data_type             <- c("Copy number segmentation")
        lFilter$experimental_strategy <- "Genotyping array"
        lFilter$platform              <- "Affymetrix SNP Array 6.0"
        lFilter$file_name             <- "\\.hg18\\.seg\\.txt"
        lFilter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (s1Assay == "cna_cnv.hg19") {
        lFilter$data_category         <- c("Copy number variation")
        lFilter$data_type             <- c("Copy number segmentation")
        lFilter$experimental_strategy <- "Genotyping array"
        lFilter$platform              <- "Affymetrix SNP Array 6.0"
        lFilter$file_name             <- "\\.hg19\\.seg\\.txt"
        lFilter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (s1Assay == "cna_nocnv.hg18") {
        lFilter$data_category         <- c("Copy number variation")
        lFilter$data_type             <- c("Copy number segmentation")
        lFilter$experimental_strategy <- "Genotyping array"
        lFilter$platform              <- "Affymetrix SNP Array 6.0"
        lFilter$file_name             <- "\\.nocnv_hg18\\.seg\\.txt"
        lFilter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (s1Assay == "cna_nocnv.hg19") {
        lFilter$data_category         <- c("Copy number variation")
        lFilter$data_type             <- c("Copy number segmentation")
        lFilter$experimental_strategy <- "Genotyping array"
        lFilter$platform              <- "Affymetrix SNP Array 6.0"
        lFilter$file_name             <- "\\.nocnv_hg19\\.seg\\.txt"
        lFilter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (s1Assay == "exon_RNAseq") {
        lFilter$data_category         <- c("Gene expression")
        lFilter$data_type             <- "Exon quantification"
        lFilter$experimental_strategy <- "RNA-Seq"
        lFilter$platform              <- "Illumina HiSeq"
        lFilter$file_name             <- "\\.bt\\.exon_quantification\\.txt"
        lFilter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (s1Assay == "exon_TotalRNAseq") {
        lFilter$data_category         <- c("Gene expression")
        lFilter$data_type             <- "Exon quantification"
        lFilter$experimental_strategy <- "Total RNA-Seq"
        lFilter$platform              <- "Illumina HiSeq"
        lFilter$file_name             <- "\\.bt\\.exon_quantification\\.txt"
        lFilter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (s1Assay == "exonJunction_RNAseq") {
        lFilter$data_category         <- c("Gene expression")
        lFilter$data_type             <- "Exon junction quantification"
        lFilter$experimental_strategy <- "RNA-Seq"
        lFilter$platform              <- "Illumina HiSeq"
        lFilter$file_name             <- "\\.junction_quantification\\.txt"
        lFilter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (s1Assay == "exonJunction_TotalRNAseq") {
        lFilter$data_category         <- c("Gene expression")
        lFilter$data_type             <- "Exon junction quantification"
        lFilter$experimental_strategy <- "Total RNA-Seq"
        lFilter$platform              <- "Illumina HiSeq"
        lFilter$file_name             <- "\\.junction_quantification\\.txt"
        lFilter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (s1Assay == "gene_Array") {
        lFilter$data_category         <- c("Gene expression")
        lFilter$data_type             <- "Gene expression quantification"
        lFilter$experimental_strategy <- "Gene expression array"
        lFilter$platform              <- "AgilentG4502A_07_3"
        lFilter$file_name             <- "\\.txt_lmean\\.out\\.logratio\\.gene\\.tcga_level3\\.data\\.txt"
        lFilter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (s1Assay == "gene.normalized_RNAseq") {
        lFilter$data_category         <- c("Gene expression")
        lFilter$data_type             <- "Gene expression quantification"
        lFilter$experimental_strategy <- "RNA-Seq"
        lFilter$platform              <- "Illumina HiSeq"
        lFilter$file_name             <- "\\.rsem\\.genes\\.normalized_results"
        lFilter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (s1Assay == "gene_RNAseq") {
        lFilter$data_category         <- c("Gene expression")
        lFilter$data_type             <- "Gene expression quantification"
        lFilter$experimental_strategy <- "RNA-Seq"
        lFilter$platform              <- "Illumina HiSeq"
        lFilter$file_name             <- "\\.rsem\\.genes\\.results"
        lFilter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (s1Assay == "gene.normalized_TotalRNAseq") {
        lFilter$data_category         <- c("Gene expression")
        lFilter$data_type             <- "Gene expression quantification"
        lFilter$experimental_strategy <- "Total RNA-Seq"
        lFilter$platform              <- "Illumina HiSeq"
        lFilter$file_name             <- "\\.rsem\\.genes\\.normalized_results"
        lFilter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (s1Assay == "gene_TotalRNAseq") {
        lFilter$data_category         <- c("Gene expression")
        lFilter$data_type             <- "Gene expression quantification"
        lFilter$experimental_strategy <- "Total RNA-Seq"
        lFilter$platform              <- "Illumina HiSeq"
        lFilter$file_name             <- "\\.rsem\\.genes\\.results"
        lFilter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (s1Assay == "isoform.normalized_RNAseq") {
        lFilter$data_category         <- c("Gene expression")
        lFilter$data_type             <- "Isoform expression quantification"
        lFilter$experimental_strategy <- "RNA-Seq"
        lFilter$platform              <- "Illumina HiSeq"
        lFilter$file_name             <- "\\.rsem\\.isoforms\\.normalized_results"
        lFilter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (s1Assay == "isoform_RNAseq") {
        lFilter$data_category         <- c("Gene expression")
        lFilter$data_type             <- "Isoform expression quantification"
        lFilter$experimental_strategy <- "RNA-Seq"
        lFilter$platform              <- "Illumina HiSeq"
        lFilter$file_name             <- "\\.rsem\\.isoforms\\.results"
        lFilter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (s1Assay == "isoform.normalized_TotalRNAseq") {
        lFilter$data_category         <- c("Gene expression")
        lFilter$data_type             <- "Isoform expression quantification"
        lFilter$experimental_strategy <- "Total RNA-Seq"
        lFilter$platform              <- "Illumina HiSeq"
        lFilter$file_name             <- "\\.rsem\\.isoforms\\.normalized_results"
        lFilter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (s1Assay == "isoform_TotalRNAseq") {
        lFilter$data_category         <- c("Gene expression")
        lFilter$data_type             <- "Isoform expression quantification"
        lFilter$experimental_strategy <- "Total RNA-Seq"
        lFilter$platform              <- "Illumina HiSeq"
        lFilter$file_name             <- "\\.rsem\\.isoforms\\.results"
        lFilter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (s1Assay == "mir_GA.hg18") {
        lFilter$data_category         <- c("Gene expression")
        lFilter$data_type             <- "miRNA gene quantification"
        lFilter$experimental_strategy <- "miRNA-Seq"
        lFilter$platform              <- "Illumina GA"
        lFilter$file_name             <- "^[^\\.]*\\.mirna\\.quantification\\.txt"
        lFilter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (s1Assay == "mir_GA.hg19") {
        lFilter$data_category         <- c("Gene expression")
        lFilter$data_type             <- "miRNA gene quantification"
        lFilter$experimental_strategy <- "miRNA-Seq"
        lFilter$platform              <- "Illumina GA"
        lFilter$file_name             <- "^[^\\.]*\\.hg19\\.mirna\\.quantification\\.txt"
        lFilter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (s1Assay == "mir_GA.hg19.mirbase20") {
        lFilter$data_category         <- c("Gene expression")
        lFilter$data_type             <- "miRNA gene quantification"
        lFilter$experimental_strategy <- "miRNA-Seq"
        lFilter$platform              <- "Illumina GA"
        lFilter$file_name             <- "^[^\\.]*\\.hg19\\.mirbase20\\.mirna\\.quantification\\.txt"
        lFilter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (s1Assay == "mir_HiSeq.hg18") {
        lFilter$data_category         <- c("Gene expression")
        lFilter$data_type             <- "miRNA gene quantification"
        lFilter$experimental_strategy <- "miRNA-Seq"
        lFilter$platform              <- "Illumina HiSeq"
        lFilter$file_name             <- "^[^\\.]*\\.mirna\\.quantification\\.txt"
        lFilter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (s1Assay == "mir_HiSeq.hg19") {
        lFilter$data_category         <- c("Gene expression")
        lFilter$data_type             <- "miRNA gene quantification"
        lFilter$experimental_strategy <- "miRNA-Seq"
        lFilter$platform              <- "Illumina HiSeq"
        lFilter$file_name             <- "^[^\\.]*\\.hg19\\.mirna\\.quantification\\.txt"
        lFilter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (s1Assay == "mir_HiSeq.hg19.mirbase20") {
        lFilter$data_category         <- c("Gene expression")
        lFilter$data_type             <- "miRNA gene quantification"
        lFilter$experimental_strategy <- "miRNA-Seq"
        lFilter$platform              <- "Illumina HiSeq"
        lFilter$file_name             <- "^[^\\.]*\\.hg19\\.mirbase20\\.mirna\\.quantification\\.txt"
        lFilter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (s1Assay == "mirIsoform_GA.hg18") {  # rows different
        lFilter$data_category         <- c("Gene expression")
        lFilter$data_type             <- "miRNA isoform quantification"
        lFilter$experimental_strategy <- "miRNA-Seq"
        lFilter$platform              <- "Illumina GA"
        lFilter$file_name             <- "^[^\\.]*\\.isoform\\.quantification\\.txt"
        lFilter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (s1Assay == "mirIsoform_GA.hg19") {  # rows different
        lFilter$data_category         <- c("Gene expression")
        lFilter$data_type             <- "miRNA isoform quantification"
        lFilter$experimental_strategy <- "miRNA-Seq"
        lFilter$platform              <- "Illumina GA"
        lFilter$file_name             <- "^[^\\.]*\\.hg19\\.isoform\\.quantification\\.txt"
        lFilter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (s1Assay == "mirIsoform_GA.hg19.mirbase20") {  # rows different
        lFilter$data_category         <- c("Gene expression")
        lFilter$data_type             <- "miRNA isoform quantification"
        lFilter$experimental_strategy <- "miRNA-Seq"
        lFilter$platform              <- "Illumina GA"
        lFilter$file_name             <- "^[^\\.]*\\.hg19\\.mirbase20\\.isoform\\.quantification\\.txt"
        lFilter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (s1Assay == "mirIsoform_HiSeq.hg18") {  # rows different
        lFilter$data_category         <- c("Gene expression")
        lFilter$data_type             <- "miRNA isoform quantification"
        lFilter$experimental_strategy <- "miRNA-Seq"
        lFilter$platform              <- "Illumina HiSeq"
        lFilter$file_name             <- "^[^\\.]*\\.isoform\\.quantification\\.txt"
        lFilter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (s1Assay == "mirIsoform_HiSeq.hg19") {  # rows different
        lFilter$data_category         <- c("Gene expression")
        lFilter$data_type             <- "miRNA isoform quantification"
        lFilter$experimental_strategy <- "miRNA-Seq"
        lFilter$platform              <- "Illumina HiSeq"
        lFilter$file_name             <- "^[^\\.]*\\.hg19\\.isoform\\.quantification\\.txt"
        lFilter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (s1Assay == "mirIsoform_HiSeq.hg19.mirbase20") {  # rows different
        lFilter$data_category         <- c("Gene expression")
        lFilter$data_type             <- "miRNA isoform quantification"
        lFilter$experimental_strategy <- "miRNA-Seq"
        lFilter$platform              <- "Illumina HiSeq"
        lFilter$file_name             <- "^[^\\.]*\\.hg19\\.mirbase20\\.isoform\\.quantification\\.txt"
        lFilter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (s1Assay == "methylation_27") {
        lFilter$data_category         <- "DNA methylation"
        lFilter$data_type             <- "Methylation beta value"
        lFilter$experimental_strategy <- "Methylation array"
        lFilter$platform              <- "Illumina Human Methylation 27"
        lFilter$file_name             <- "jhu-usc\\.edu_.*\\.HumanMethylation27\\."
        lFilter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (s1Assay == "methylation_450") {
        lFilter$data_category         <- "DNA methylation"
        lFilter$data_type             <- "Methylation beta value"
        lFilter$experimental_strategy <- "Methylation array"
        lFilter$platform              <- "Illumina Human Methylation 450"
        lFilter$file_name             <- "jhu-usc\\.edu_.*\\.HumanMethylation450\\."
        lFilter$submitter_id          <- "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"
    } else if (s1Assay == "protein_RPPA") {  # differnt "archive.file_name"
        lFilter$data_category         <- "Protein expression"
        lFilter$data_type             <- "Protein expression quantification"
        lFilter$experimental_strategy <- "Protein expression array"
        lFilter$platform              <- "MDA_RPPA_Core"
        lFilter$file_name             <-
            "mdanderson\\.org_.*\\.MDA_RPPA_Core\\.protein_expression\\.Level_3\\."
        lFilter$submitter_id          <- "cases.0.samples.0.portions.0.submitter_id"
    } else if (s1Assay == "somaticMutation_DNAseq") {
        lFilter$data_category         <- "Simple nucleotide variation"
        lFilter$data_type             <- "Simple somatic mutation"
        lFilter$experimental_strategy <- "DNA-Seq"
        lFilter$platform              <- c("Illumina GA", "Illumina HiSeq", "Mixed platforms")
        lFilter$file_name             <- "\\.somatic\\.maf"
        lFilter$submitter_id          <- c("cases.0.samples.0.portions.0.submitter_id",  # both exist
                                           "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id")
    } else if (s1Assay == "glycoproteome_iTRAQ") {  # from CTPAC not GDC
        lFilter$data_category         <- NA
        lFilter$data_type             <- NA
        lFilter$experimental_strategy <- NA
        lFilter$platform              <- "glycoproteome_iTRAQ"
        lFilter$file_name             <- "_Glycoproteome\\.glycosite\\.itraq\\.tsv"
        lFilter$submitter_id          <- NA
    } else if (s1Assay == "phosphoproteome_iTRAQ") {  # from CTPAC not GDC
        lFilter$data_category         <- NA
        lFilter$data_type             <- NA
        lFilter$experimental_strategy <- NA
        lFilter$platform              <- "phosphoproteome_iTRAQ"
        lFilter$file_name             <- "_Phosphoproteome\\.phosphosite\\.itraq\\.tsv"
        lFilter$submitter_id          <- NA
    } else if (s1Assay == "proteome_iTRAQ") {  # from CTPAC not GDC
        lFilter$data_category         <- NA
        lFilter$data_type             <- NA
        lFilter$experimental_strategy <- NA
        lFilter$platform              <- "proteome_iTRAQ"
        lFilter$file_name             <- "_Proteome.(itraq|spectral_counts)\\.tsv"
        lFilter$submitter_id          <- NA
    }
    return(lFilter)
}

#' Get metadata of biospecimen & clinical data files with defined filter
#'
#' @param sDir String of directory for temporary files.
#' @param sArchive String of archive type: "legacy" or "".
#' @param sFieldMeta Character vector of colnames in metadata.
#' @param nEntityCount Number of entity (row) in the metadata: "-1" means all.
#' @param sEndpoint String of endpoint: "files".
#' @return A \code{data.frame} of metadata.
#' @example
#' dMetaDataClin <- MetaDataClin(sDir = ".")
MetaDataClin <- function(sDir = ".",
                         sArchive = "legacy",
                         sFieldMeta = "",
                         nEntityCount = (-1),
                         sEndpoint = "files") {
    library(rjson)
    if (nEntityCount == (-1)) {
        nEntityCount <- EntityCount(sArchive, sEndpoint)
    }
    if (sFieldMeta == "") {
        sFieldMeta <- FieldsMeta()
    } else {
        sFieldList <- FieldsList(sArchive, sEndpoint)
        print("fields names: checking ...")
        stopifnot(all(sFieldMeta %in% sFieldList))
        print("fields names: checking done!")
    }
    print("metadata file: preparing ...")
    sOutfile <- paste(sDir, "/tmp_metadata_", sEndpoint, ".tsv", sep = "")
    sUrl <- paste('"https://gdc-api.nci.nih.gov/', sArchive,
                  ifelse(sArchive == '', '', '/'), sEndpoint, '"',
                  sep = '')
    sOption <- paste("-o ", sOutfile,
                     " --silent --show-error --request POST",
                     " --header Content-Type:application/json --data @",
                     sDir, "/tmp_metadata.json",
                     sep = "")
    sArgument <- paste(sOption, sUrl)
    lFilter1 <- list()
    lFilter1$access <- list(op = "=",
                            content = list(field = "access", value = "open"))
    lFilter1$data_format <- list(op = "in",
                                 content = list(field = "data_format",
                                                value = "Biotab"))
    lFilterAll <- list(op = "and",
                       content = list(lFilter1$access,
                                      lFilter1$data_format))
    lPayload <- list(filters = lFilterAll,
                     format = "TSV",
                     sort = "file_id",
                     from = 1,
                     size = nEntityCount,
                     fields = paste(sFieldMeta, collapse = ","))
    cat(toJSON(lPayload), file = paste(sDir, "/tmp_metadata.json", sep = ""))
    stdOut <- system2("curl", sArgument, stdout = TRUE)
    if (!is.null(attr(stdOut, "status"))) {
        print("error (download): check the proxy")
    }
    stopifnot(is.null(attr(stdOut, "status")))
    if ("data" %in% dir()) {file.remove("data")}
    dMetaData <- tryCatch(read.csv(sOutfile,
                                   sep = "\t",
                                   row.names = NULL,
                                   as.is = TRUE,
                                   na.strings = "",
                                   check.names = FALSE),
                          error = function(e){print(e$message); return(NULL)})
    if (!is.null(dMetaData)) {  # > 0 lines available in input
        if (nrow(dMetaData) == 0) {
            dMetaData <- NULL
        } else {
            write.table(dMetaData,
                        file = paste(sDir,
                                     "/tmp_metadata_clinical__",
                                     sArchive,
                                     "__",
                                     sEndpoint,
                                     ".tsv",
                                     sep = ""),
                        quote = FALSE,
                        sep = "\t",
                        col.names = NA,
                        row.names = TRUE)
        }
    }
    print("metadata file: preparing done!")
    return(dMetaData)
}

#' Get metadata of somatic mutation data files with defined filter
#'
#' @param sCancer String of cancer type.
#' @param s1Assay String of assay platform, used in \code{Filter}.
#' @param sSampleTypeId Character vector of sample_type_id: "01", ..., "14", etc.
#' @param sDir String of directory for temporary files.
#' @param sArchive String of archive type: "legacy" or "".
#' @param sFieldMeta Character vector of colnames in metadata.
#' @param nEntityCount Number of entity (row) in the metadata: "-1" means all.
#' @param sEndpoint String of endpoint: "files".
#' @return A \code{data.frame} of metadata.
MetaDataSoma <- function(sCancer = "BRCA",
                         s1Assay,
                         sSampleTypeId = sprintf("%02d", c(seq(14), 20, 40, 50, 60, 61)),
                         sDir = ".",
                         sArchive = "legacy",
                         sFieldMeta = "",
                         nEntityCount = (-1),
                         sEndpoint = "files") {
    library(rjson)
    if (nEntityCount == (-1)) {
        nEntityCount <- EntityCount(sArchive, sEndpoint)
    }
    if (sFieldMeta == "") {
        sFieldMeta <- FieldsMeta()
    } else {
        sFieldList <- FieldsList(sArchive, sEndpoint)
        print("fields names: checking ...")
        stopifnot(all(sFieldMeta %in% sFieldList))
        print("fields names: checking done!")
    }
    lFilter <- Filter(s1Assay)
    print("metadata file: preparing ...")
    sOutfile <- paste(sDir, "/tmp_metadata_", sEndpoint, ".tsv", sep = "")
    sUrl <- paste('"https://gdc-api.nci.nih.gov/', sArchive,
                  ifelse(sArchive == '', '', '/'), sEndpoint, '"',
                  sep = '')
    sOption <- paste("-o ", sOutfile, " --silent --show-error --request POST ",
                     "--header Content-Type:application/json --data @",
                     sDir, "/tmp_metadata.json", sep = "")
    sArgument <- paste(sOption, sUrl)
    lFilter1 <- list()
    lFilter1$access <- list(op = "=" ,
                            content = list(field = "access",
                                           value = "open"))
    lFilter1$data_format <- list(op = "in",
                                 content = list(field = "data_format",
                                                value = "MAF"))
    lFilter1$project_id <- list(op = "in",
                                content = list(field = "cases.project.project_id",
                                               value = paste("TCGA", sCancer,
                                                             sep = "-")))
    lFilter1$data_category <- list(op = "in",
                                   content = list(field = "data_category",
                                                  value = lFilter$data_category))
    lFilter1$data_type <- list(op = "in",
                               content = list(field = "data_type",
                                              value = lFilter$data_type))
    lFilter1$experimental_strategy <- list(op = "in",
                                           content = list(field = "experimental_strategy",
                                                          value = lFilter$experimental_strategy))
    lFilter1$platform <- list(op = "in",
                              content = list(field = "platform",
                                             value = lFilter$platform))
    lFilterAll <- list(op = "and",
                       content = list(lFilter1$access,
                                      lFilter1$data_format,
                                      lFilter1$project_id,
                                      lFilter1$data_category,
                                      lFilter1$data_type,
                                      lFilter1$experimental_strategy,
                                      lFilter1$platform))
    lPayload <- list(filters = lFilterAll,
                     format = "TSV",
                     sort = "file_id",
                     from = 1,
                     size = nEntityCount,
                     fields = paste(sFieldMeta, collapse = ","))
    cat(toJSON(lPayload), file = paste(sDir, "/tmp_metadata.json", sep = ""))
    stdOut <- system2("curl", sArgument, stdout = TRUE)
    if (!is.null(attr(stdOut, "status"))) {
        print("error (download): check the proxy")
    }
    stopifnot(is.null(attr(stdOut, "status")))
    if ("data" %in% dir()) {file.remove("data")}
    dMetaData <- tryCatch(read.csv(sOutfile,
                                   sep = "\t",
                                   row.names = NULL,
                                   as.is = TRUE,
                                   na.strings = "",
                                   check.names = FALSE),
                          error = function(e){print(e$message); return(NULL)})
    if (!is.null(dMetaData)) {  # > 0 lines available in input
        sColNames <- c("archive.file_name",
                       "cases.0.project.project_id",
                       "data_category",
                       "data_type",
                       "experimental_strategy",
                       "file_id",
                       "file_name",
                       "file_size",
                       "platform",
                       "updated_datetime")
        dMetaData <- dMetaData[, sColNames]
        rownames(dMetaData) <- dMetaData[, "file_id"]
        write.table(dMetaData,
                    file = paste(sDir,
                                 "/tmp_metadata_byType__",
                                 sArchive,
                                 "__",
                                 sEndpoint,
                                 "__",
                                 paste(sCancer, collapse = "_"),
                                 "__",
                                 s1Assay,
                                 ".tsv",
                                 sep = ""),
                    quote = FALSE,
                    sep = "\t",
                    col.names = NA,
                    row.names = TRUE)
    }
    print("metadata file: preparing done!")
    return(dMetaData)
}


#' Get metadata of spefified assay platform files with defined filter
#'
#' @param sCancer String of cancer type.
#' @param s1Assay String of assay platform, used in \code{Filter}.
#' @param sSampleTypeId Character vector of sample_type_id: "01", ..., "14", etc.
#' @param sDir String of directory for temporary files.
#' @param sArchive String of archive type: "legacy" or "".
#' @param sFieldMeta Character vector of colnames in metadata.
#' @param nEntityCount Number of entity (row) in the metadata: "-1" means all.
#' @param sEndpoint String of endpoint: "files".
#' @return A \code{data.frame} of metadata.
#' @example
#' dMetaData <- MetaData(sCancer = "BRCA",
#'                       s1Assay = "gene_RNAseq",
#'                       sSampleTypeId = sprintf("%02d", c(seq(14), 20, 40, 50, 60, 61)),
#'                       sDir = ".",
#'                       sArchive = "legacy",
#'                       sFieldMeta = "",
#'                       nEntityCount = (-1),
#'                       sEndpoint = "files")
MetaData <- function(sCancer = "BRCA",
                     s1Assay,
                     sSampleTypeId = sprintf("%02d", c(seq(14), 20, 40, 50, 60, 61)),
                     sDir = ".",
                     sArchive = "legacy",
                     sFieldMeta = "",
                     nEntityCount = (-1),
                     sEndpoint = "files") {
    library(rjson)
    if (nEntityCount == (-1)) {
        nEntityCount <- EntityCount(sArchive, sEndpoint)
    }
    if (sFieldMeta == "") {
        sFieldMeta <- FieldsMeta()
    } else {
        sFieldList <- FieldsList(sArchive, sEndpoint)
        print("fields names: checking ...")
        stopifnot(all(sFieldMeta %in% sFieldList))
        print("fields names: checking done!")
    }
    lFilter <- Filter(s1Assay)
    print("metadata file: preparing ...")
    sOutfile <- paste(sDir, "/tmp_metadata_", sEndpoint, ".tsv", sep = "")
    sUrl <- paste('"https://gdc-api.nci.nih.gov/', sArchive,
                  ifelse(sArchive == '', '', '/'), sEndpoint, '"',
                  sep = '')
    sOption <- paste("-o ", sOutfile,
                     " --silent --show-error --request POST",
                     " --header Content-Type:application/json --data @",
                     sDir, "/tmp_metadata.json",
                     sep = "")
    sArgument <- paste(sOption, sUrl)
    lFilter1 <- list()
    lFilter1$access <- list(op = "=" ,
                            content = list(field = "access",
                                           value = "open"))
    lFilter1$data_format <- list(op = "in",
                                 content = list(field = "data_format",
                                                value = "TXT"))
    lFilter1$project_id <- list(op = "in",
                                content = list(field = "cases.project.project_id",
                                               value = paste("TCGA", sCancer,
                                                             sep = "-")))
    lFilter1$data_category <- list(op = "in",
                                   content = list(field = "data_category",
                                                  value = lFilter$data_category))
    lFilter1$data_type <- list(op = "in",
                               content = list(field = "data_type",
                                              value = lFilter$data_type))
    lFilter1$experimental_strategy <- list(op = "in",
                                           content = list(field = "experimental_strategy",
                                                          value = lFilter$experimental_strategy))
    lFilter1$platform <- list(op = "in",
                              content = list(field = "platform",
                                             value = lFilter$platform))
    names(sSampleTypeId) <- NULL  # avoid names added into tmp_metadata.json
    sFieldSampleTypeId <- "cases.samples.sample_type_id"
    lFilterSampleTypeId <- list(op = "in",
                                content = list(field = sFieldSampleTypeId,
                                               value = sSampleTypeId))
    lFilterAll <- list(op = "and",
                       content = list(lFilter1$access,
                                      lFilter1$data_format,
                                      lFilter1$project_id,
                                      lFilter1$data_category,
                                      lFilter1$data_type,
                                      lFilter1$experimental_strategy,
                                      lFilter1$platform,
                                      lFilterSampleTypeId))
    lPayload <- list(filters = lFilterAll,
                     format = "TSV",
                     sort = "file_id",
                     from = 1,
                     size = nEntityCount,
                     fields = paste(sFieldMeta, collapse = ","))
    cat(toJSON(lPayload), file = paste(sDir, "/tmp_metadata.json", sep = ""))
    stdOut <- system2("curl", sArgument, stdout = TRUE)
    if (!is.null(attr(stdOut, "status"))) {
        print("error (download): check the proxy")
    }
    stopifnot(is.null(attr(stdOut, "status")))
    if ("data" %in% dir()) {file.remove("data")}
    dMetaData <- tryCatch(read.csv(sOutfile,
                                   sep = "\t",
                                   row.names = NULL,
                                   as.is = TRUE,
                                   na.strings = "",
                                   colClasses = "character",
                                   check.names = FALSE),
                          error = function(e){print(e$message); return(NULL)})
    if (!is.null(dMetaData)) {  # > 0 lines available in input
        if (length(grep("\\.1\\.", colnames(dMetaData))) > 0) {
            dMetaData <- MetaCut(dMetaData, sColPattern = "\\.1\\.")
        }
        dMetaCut <- dMetaData[grep(lFilter$file_name, dMetaData$file_name),
                              sort(colnames(dMetaData))]  # filter with filename
        if (nrow(dMetaCut) == 0) {
            dMetaData <- NULL
        } else {
            dMetaCut <- dMetaCut[order(dMetaCut$file_name), ]  # sort before split !!
            # duplicated files with same barcode
            if (any(table(dMetaCut[, lFilter$submitter_id]) != 1)) {
                bRow <- ArchiveNewestInGroups(dMetaCut$archive.file_name,
                                              dMetaCut$file_name)
                # not lFilter$submitter_id because of one barcode to multi files
                dMetaCut <- dMetaCut[bRow, ]
            }
            # duplicated files with same barcode
            stopifnot(all(table(dMetaCut[, lFilter$submitter_id]) == 1))
            rownames(dMetaCut) <- dMetaCut[, lFilter$submitter_id]
            # sort by barcode
            dMetaData <- dMetaCut[order(rownames(dMetaCut)), order(colnames(dMetaCut))]
            stopifnot(all(rownames(dMetaData) == dMetaData[lFilter$submitter_id]))
            write.table(dMetaData,
                        file = paste(sDir,
                                     "/tmp_metadata_byType__",
                                     sArchive,
                                     "__",
                                     sEndpoint,
                                     "__",
                                     paste(sCancer, collapse = "_"),
                                     "__",
                                     s1Assay,
                                     ".tsv",
                                     sep = ""),
                        quote = FALSE,
                        sep = "\t",
                        col.names = NA,
                        row.names = TRUE)
        }
    }
    print("metadata file: preparing done!")
    return(dMetaData)
}

#' Download files by \code{file_id}
#'
#' @param sFileId2Bar Character vector of barcode with file_id as the name.
#' @param sDir String of directory for temporary files.
#' @param sArchive String of archive type: "legacy" or "".
#' @return Character vector of downloaded filename with file_id as the name.
#' @example
#' sFileNameById <- FileNameById(dMetaData$file_id, sDir = ".")
FileNameById <- function(sFileId2Bar,
                         sDir = ".",
                         sArchive = "legacy") {
    if (length(sFileId2Bar) == 0) {  # no files available for this filter
        sFileNameById <- NULL
    } else {
        stopifnot(length(sFileId2Bar) > 0)  # no files available for this filter
        library(rjson)
        sTar <- paste(sDir, "/gdc_download_", TimeNow(), ".tar.gz", sep = "")
        sUrl <- paste("https://gdc-api.nci.nih.gov/",
                      sArchive,
                      ifelse(sArchive == "", "", "/"),
                      "data",
                      sep = "")
        sOption <- paste("-o ", sTar,
                         " --silent --show-error --request POST ",
                         "--header Content-Type:application/json --data @",
                         sDir, "/tmp_id.json",
                         sep = "")
        sArgument <- paste(sOption, sUrl)
        cat(toJSON(list(ids = names(sFileId2Bar))),
            file = paste(sDir, "/tmp_id.json", sep = ""))
        print("*.tar.gz file: downloading & unzipping ...")
        sErr <- "error"
        while (sErr != 0) {
            stdOut <- system2("curl", sArgument, stdout = TRUE)
            if (!is.null(attr(stdOut, "status"))) {
                print("error (download): check the proxy")
            }
            stopifnot(is.null(attr(stdOut, "status")))
            sUntar <- strsplit(sTar, split = "\\.")[[1]][1]
            sErr <- tryCatch(untar(sTar, exdir = sUntar, tar = "internal"),
                             error = function(e){return(e$message)})
            if (length(sFileId2Bar) == 1) {  # only 1 unzipped file in *.tar.gz
                file.rename(from = sTar,
                            to = paste(sUntar, "/", sFileId2Bar[1], ".tsv", sep = ""))
                dManifest <- data.frame(filename = paste(sFileId2Bar[1],
                                                         ".tsv",
                                                         sep = ""))
                write.table(dManifest,
                            file = paste(sUntar, "/MANIFEST.txt", sep = ""),
                            quote = FALSE,
                            sep = "\t",
                            col.names = TRUE,
                            row.names = FALSE)
                sErr <- 0
            }
        }
        dManifest <- read.csv(paste(sUntar, "MANIFEST.txt", sep = "/"),
                              sep = "\t",
                              row.names = NULL,
                              as.is = TRUE,
                              na.strings = "",
                              check.names = FALSE)
        sFileNameById <- paste(sUntar, dManifest$filename, sep = "/")
        names(sFileNameById) <- sFileId2Bar[dManifest$id]
    }
    if ("data" %in% dir()) {file.remove("data")}
    print("*.tar.gz file: downloading & unzipping done!")
    return(sFileNameById[order(names(sFileNameById))])  # sort by file_id
}


#  =============================================================================
#  merge functions, NOT used directly by user
#  =============================================================================

#' Merge copy number variations ("cna") files, distributed by \code{Merge}
#'
#' @param sFileNameById Character vector of filename (named with file_id).
#' @return A \code{data.frame} of merged table.
MergeCopy <- function(sFileNameById) {
    sColNames <- c("Chromosome", "Start", "End", "Num_Probes", "Segment_Mean")
    ld <- ColumnsFromFiles(sFileNameById,
                           sColNames,
                           sSortByCol = "",
                           sSkipLines = 0,
                           sNa = "NA")
    dMerged <- do.call(rbind, ld)
    dMerged <- as.data.frame(as.matrix(dMerged), stringsAsFactors = FALSE)
    return(dMerged)
}

#' Merge methyloation ("methy") files, distributed by \code{Merge}
#'
#' @param sFileNameById Character vector of filename (named with file_id).
#' @return A \code{data.frame} of merged table.
MergeMethy <- function(sFileNameById) {
    sValue <- c("Beta_value")
    sProbe <- c("Composite Element REF",
                "Gene_Symbol",
                "Chromosome",
                "Genomic_Coordinate")
    for (n in seq(length(sFileNameById))) {
        d <- ColumnsFromFile(sFileNameById[n],
                             names(sFileNameById)[n],
                             c(sProbe, sValue),
                             sSortByCol = sProbe[1],
                             sSkipLines = 1,
                             sNa = "NA")
        if (n == 1) {
            sName <- sort(rownames(d))
            stopifnot(length(sName) == length(unique(sName)))
            dProbe <- d[sName, sProbe]
            m <- matrix(nrow = nrow(d), ncol = length(sFileNameById))
            m[, n] <- d[sName, sValue]
        } else {
            stopifnot(all(sort(rownames(d)) == sName))
            m[, n] <- d[sName, sValue]
        }
    }
    dMerged <- cbind(dProbe, m)
    colnames(dMerged) <- c("CpG",
                           "Gene_Symbol",
                           "Chromosome",
                           "Genomic_Coordinate",
                           names(sFileNameById))
    return(dMerged)
}

#' Merge microRNA ("mir") files, distributed by \code{Merge}
#'
#' @param sFileNameById Character vector of filename (named with file_id).
#' @return A \code{data.frame} of merged table.
MergeMir <- function(sFileNameById) {
    sValue <- c("read_count", "reads_per_million_miRNA_mapped")
    sProbe <- c("miRNA_ID")
    for (n in seq(length(sFileNameById))) {
        d <- ColumnsFromFile(sFileNameById[n],
                             names(sFileNameById)[n],
                             c(sProbe, sValue),
                             sSortByCol = sProbe[1],
                             sSkipLines = 0,
                             sNa = "NA")
        if (n == 1) {
            sName <- sort(rownames(d))
            stopifnot(length(sName) == length(unique(sName)))
            dProbe <- d[sName, sProbe]
            m <- matrix(nrow = nrow(d), ncol = length(sFileNameById) * 2)
            m[, 2 * n - 1] <- d[sName, sValue[1]]
            m[, 2 * n  ] <- d[sName, sValue[2]]
        } else {
            stopifnot(all(sort(rownames(d)) == sName))
            m[, 2 * n - 1] <- d[sName, sValue[1]]
            m[, 2 * n  ] <- d[sName, sValue[2]]
        }
    }
    dMerged <- cbind(dProbe, m)
    colnames(dMerged) <- c(sProbe,
                           paste(rep(names(sFileNameById),
                                     each = 2),
                                 rep(c("read_count",
                                       "reads_per_million_miRNA_mapped"),
                                     length(sFileNameById)),
                                 sep = ":"))
    return(dMerged)
}

#' Merge microRNA isoform ("mirIsoform") files, distributed by \code{Merge}
#'
#' @param sFileNameById Character vector of filename (named with file_id).
#' @return A \code{data.frame} of merged table.
MergeMirIso <- function(sFileNameById) {
    sValue <- c("read_count", "reads_per_million_miRNA_mapped")
    sProbe <- c("isoform_coords", "miRNA_ID", "miRNA_region", "cross-mapped")
    ld <- ColumnsFromFiles(sFileNameById,
                           c(sProbe, sValue),
                           sSortByCol = sProbe[1],
                           sSkipLines = 0,
                           sNa = "NA")
    lsRowname <- lapply(ld, function(x){paste(x[, "isoform_coords"],
                                              x[, "miRNA_ID"],
                                              x[, "miRNA_region"],
                                              x[, "cross-mapped"],
                                              sep = "|")})
    sMirIsoId <- sort(unique(unlist(lsRowname)))
    for (n in seq(length(ld))) {
        rownames(ld[[n]]) <- lsRowname[[n]]
    }
    m <- matrix(nrow = length(sMirIsoId),
                ncol = length(sFileNameById) * 2)
    for (n in seq(length(sFileNameById))) {
        m[, 2 * n - 1] <- ld[[n]][sMirIsoId, "read_count"]
        m[, 2 * n    ] <- ld[[n]][sMirIsoId, "reads_per_million_miRNA_mapped"]
    }
    dMerged <- cbind(do.call(rbind, strsplit(sMirIsoId, split = "\\|")), m)
    colnames(dMerged) <- c(sProbe,
                           paste(rep(names(sFileNameById),
                                     each = 2),
                                 rep(c("read_count",
                                       "reads_per_million_miRNA_mapped"),
                                     length(sFileNameById)),
                                 sep = ":"))
    return(dMerged)
}

#' Merge gene array ("gene_Array") files, distributed by \code{Merge}
#'
#' @param sFileNameById Character vector of filename (named with file_id).
#' @return A \code{data.frame} of merged table.
MergeGeneArray <- function(sFileNameById) {
    sValue <- c("log2 lowess normalized (cy5/cy3) collapsed by gene symbol")
    sProbe <- c("Composite Element REF")
    for (n in seq(length(sFileNameById))) {
        d <- ColumnsFromFile(sFileNameById[n],
                             names(sFileNameById)[n],
                             c(sProbe, sValue),
                             sSortByCol = sProbe[1],
                             sSkipLines = 1,
                             sNa = "NA")
        if (n == 1) {
            sName <- sort(rownames(d))
            stopifnot(length(sName) == length(unique(sName)))
            dProbe <- d[sName, sProbe]
            m <- matrix(nrow = nrow(d), ncol = length(sFileNameById))
            m[, n] <- d[sName, sValue[1]]
        } else {
            stopifnot(all(sort(rownames(d)) == sName))
            m[, n] <- d[sName, sValue[1]]
        }
    }
    dMerged <- cbind(dProbe, m)
    #dMerged <- rbind(c("gene_id", rep(names(sFileNameById), each = 1)),
    #                 dMerged)
    colnames(dMerged) <- c("gene_id", rep(names(sFileNameById), each = 1))
    return(dMerged)
}

#' Merge gene ("gene.normalized_RNAseq") files, distributed by \code{Merge}
#'
#' @param sFileNameById Character vector of filename (named with file_id).
#' @return A \code{data.frame} of merged table.
MergeGeneRnaSeqNorm <- function(sFileNameById) {
    sValue <- c("normalized_count")
    sProbe <- c("gene_id")
    for (n in seq(length(sFileNameById))) {
        d <- ColumnsFromFile(sFileNameById[n],
                             names(sFileNameById)[n],
                             c(sProbe, sValue),
                             sSortByCol = sProbe[1],
                             sSkipLines = 0,
                             sNa = "NA")
        if (n == 1) {
            sName <- sort(rownames(d))
            stopifnot(length(sName) == length(unique(sName)))
            dProbe <- d[sName, sProbe]
            m <- matrix(nrow = nrow(d), ncol = length(sFileNameById))
            m[, n] <- d[sName, sValue[1]]
        } else {
            stopifnot(all(sort(rownames(d)) == sName))
            m[, n] <- d[sName, sValue[1]]
        }
    }
    mMerged <- cbind(dProbe, m)
    colnames(mMerged) <- c(sProbe, rep(names(sFileNameById), each = 1))
    return(mMerged)
}

#' Merge unnormalized gene ("gene_RNAseq") files, distributed by \code{Merge}
#'
#' @param sFileNameById Character vector of filename (named with file_id).
#' @return A \code{data.frame} of merged table.
MergeGeneRnaSeqUnnorm <- function(sFileNameById) {
    sValue <- c("raw_count", "scaled_estimate")
    sProbe <- c("gene_id")
    for (n in seq(length(sFileNameById))) {
        d <- ColumnsFromFile(sFileNameById[n],
                             names(sFileNameById)[n],
                             c(sProbe, sValue),
                             sSortByCol = sProbe[1],
                             sSkipLines = 0,
                             sNa = "NA")
        if (n == 1) {
            sName <- sort(rownames(d))
            stopifnot(length(sName) == length(unique(sName)))
            dProbe <- d[sName, sProbe]
            m <- matrix(nrow = nrow(d), ncol = length(sFileNameById) * 2)
            m[, 2 * n - 1] <- d[sName, sValue[1]]
            m[, 2 * n    ] <- d[sName, sValue[2]]
        } else {
            stopifnot(all(sort(rownames(d)) == sName))
            m[, 2 * n - 1] <- d[sName, sValue[1]]
            m[, 2 * n    ] <- d[sName, sValue[2]]
        }
    }
    dMerged <- cbind(dProbe, m)
    colnames(dMerged) <- c(sProbe,
                           paste(rep(names(sFileNameById), each = 2),
                                 rep(c("raw_count",
                                       "scaled_estimate"),
                                     length(sFileNameById)),
                                 sep = ":"))
    return(dMerged)
}

#' Merge gene ("isoform.normalized_RNAseq") files, distributed by \code{Merge}
#'
#' @param sFileNameById Character vector of filename (named with file_id).
#' @return A \code{data.frame} of merged table.
MergeGeneIsoRnaSeqNorm <- function(sFileNameById) {
    sValue <- c("normalized_count")
    sProbe <- c("isoform_id")
    for (n in seq(length(sFileNameById))) {
        d <- ColumnsFromFile(sFileNameById[n],
                             names(sFileNameById)[n],
                             c(sProbe, sValue),
                             sSortByCol = sProbe[1],
                             sSkipLines = 0,
                             sNa = "NA")
        if (n == 1) {
            sName <- sort(rownames(d))
            stopifnot(length(sName) == length(unique(sName)))
            dProbe <- d[sName, sProbe]
            m <- matrix(nrow = nrow(d), ncol = length(sFileNameById))
            m[, n] <- d[sName, sValue[1]]
        } else {
            stopifnot(all(sort(rownames(d)) == sName))
            m[, n] <- d[sName, sValue[1]]
        }
    }
    dMerged <- cbind(dProbe, m)
    colnames(dMerged) <- c(sProbe, names(sFileNameById))
    return(dMerged)
}

#' Merge unnormalized gene ("isoform_RNAseq") files, distributed by \code{Merge}
#'
#' @param sFileNameById Character vector of filename (named with file_id).
#' @return A \code{data.frame} of merged table.
MergeGeneIsoRnaSeqUnnorm <- function(sFileNameById) {
    sValue <- c("raw_count", "scaled_estimate")
    sProbe <- c("isoform_id")
    for (n in seq(length(sFileNameById))) {
        d <- ColumnsFromFile(sFileNameById[n],
                             names(sFileNameById)[n],
                             c(sProbe, sValue),
                             sSortByCol = sProbe[1],
                             sSkipLines = 0,
                             sNa = "NA")
        if (n == 1) {
            sName <- sort(rownames(d))
            stopifnot(length(sName) == length(unique(sName)))
            dProbe <- d[sName, sProbe]
            m <- matrix(nrow = nrow(d), ncol = length(sFileNameById) * 2)
            m[, 2 * n - 1] <- d[sName, sValue[1]]
            m[, 2 * n  ] <- d[sName, sValue[2]]
        } else {
            stopifnot(all(sort(rownames(d)) == sName))
            m[, 2 * n - 1] <- d[sName, sValue[1]]
            m[, 2 * n  ] <- d[sName, sValue[2]]
        }
    }
    dMerged <- cbind(dProbe, m)
    colnames(dMerged) <- c(sProbe,
                           paste(rep(names(sFileNameById), each = 2),
                                 rep(c("raw_count", "scaled_estimate"), length(sFileNameById)),
                                 sep = ":"))
    return(dMerged)
}

#' Merge exon ("exon") files, distributed by \code{Merge}
#'
#' @param sFileNameById Character vector of filename (named with file_id).
#' @return A \code{data.frame} of merged table.
MergeExon <- function(sFileNameById) {
    sValue <- c("RPKM")
    sProbe <- c("exon")
    for (n in seq(length(sFileNameById))) {
        d <- ColumnsFromFile(sFileNameById[n],
                             names(sFileNameById)[n],
                             c(sProbe,
                               sValue),
                             sSortByCol = sProbe[1],
                             sSkipLines = 0,
                             sNa = "NA")
        if (n == 1) {
            sName <- sort(rownames(d))
            stopifnot(length(sName) == length(unique(sName)))
            dProbe <- d[sName, sProbe]
            m <- matrix(nrow = nrow(d), ncol = length(sFileNameById))
            m[, n] <- d[sName, sValue[1]]
        } else {
            stopifnot(all(sort(rownames(d)) == sName))
            m[, n] <- d[sName, sValue[1]]
        }
    }
    dMerged <- cbind(dProbe, m)
    colnames(dMerged) <- c(sProbe, names(sFileNameById))
    return(dMerged)
}

#' Merge exon ("exonJunction") files, distributed by \code{Merge}
#'
#' @param sFileNameById Character vector of filename (named with file_id).
#' @return A \code{data.frame} of merged table.
MergeExonJunction <- function(sFileNameById) {  # duplicated rows
    sValue <- c("raw_counts")
    sProbe <- c("junction")
    for (n in seq(length(sFileNameById))) {
        d <- ColumnsFromFile(sFileNameById[n],
                             names(sFileNameById)[n],
                             c(sProbe, sValue),
                             sSortByCol = sProbe[1],
                             sSkipLines = 0,
                             sNa = "NA")
        if (n == 1) {
            sName <- sort(rownames(d))
            stopifnot(length(sName) == length(unique(sName)))
            dProbe <- d[sName, sProbe]
            m <- matrix(nrow = nrow(d), ncol = length(sFileNameById))
            m[, n] <- d[sName, sValue[1]]
        } else {
            stopifnot(all(sort(rownames(d)) == sName))
            m[, n] <- d[sName, sValue[1]]
        }
    }
    dMerged <- cbind(dProbe, m)
    colnames(dMerged) <- c("exonJunction", names(sFileNameById))
    return(dMerged)
}

#' Merge protein ("protein_RPPA") files, distributed by \code{Merge}
#'
#' @param sFileNameById Character vector of filename (named with file_id).
#' @return A \code{data.frame} of merged table.
MergeProtein <- function(sFileNameById) {  # duplicated rows
    sPr <- c("ABL1|c-Abl",
             "ACACA ACACB|ACC_pS79",
             "ACACA|ACC1",
             "ACVRL1|ACVRL1",
             "ADAR|ADAR1",
             "AKT1 AKT2 AKT3|Akt",
             "AKT1 AKT2 AKT3|Akt_pS473",
             "AKT1 AKT2 AKT3|Akt_pT308",
             "AKT1S1|PRAS40_pT246",
             "ANXA1|Annexin-1",
             "ANXA7|Annexin_VII",
             "ARAF|A-Raf",
             "ARAF|A-Raf_pS299",
             "ARID1A|ARID1A",
             "AR|AR",
             "ASNS|ASNS",
             "ATM|ATM",
             "ATP5A1|Oxphos-complex-V_subunitb",
             "AXL|Axl",
             "BAD|Bad_pS112",
             "BAK1|Bak",
             "BAP1|Bap1-c-4",
             "BAX|Bax",
             "BCL2A1|Bcl2A1",
             "BCL2L11|Bim",
             "BCL2L1|Bcl-xL",
             "BCL2|Bcl-2",
             "BECN1|Beclin",
             "BID|Bid",
             "BIRC2|cIAP",
             "BRAF|B-Raf",
             "BRAF|B-Raf_pS445",
             "BRCA2|BRCA2",
             "BRD4|BRD4",
             "CA9|CA9",
             "CASP3|Caspase-3",
             "CASP7|Caspase-7_cleavedD198",
             "CASP8|Caspase-8",
             "CASP9|Caspase-9",
             "CAV1|Caveolin-1",
             "CCNB1|Cyclin_B1",
             "CCND1|Cyclin_D1",
             "CCNE1|Cyclin_E1",
             "CCNE2|Cyclin_E2",
             "CD274|CD274",
             "CD274|PD-L1",
             "CDH1|E-Cadherin",
             "CDH2|N-Cadherin",
             "CDH3|P-Cadherin",
             "CDK1|CDK1",
             "CDK1|CDK1_pY15",
             "CDKN1A|p21",
             "CDKN1B|p27",
             "CDKN1B|p27_pT157",
             "CDKN1B|p27_pT198",
             "CDKN2A|p16_INK4a",
             "CHEK1|Chk1",
             "CHEK1|Chk1_pS296",
             "CHEK1|Chk1_pS345",
             "CHEK2|Chk2",
             "CHEK2|Chk2_pT68",
             "CHGA|Chromogranin-A-N-term",
             "CLDN7|Claudin-7",
             "COG3|COG3",
             "COL6A1|Collagen_VI",
             "COPS5|JAB1",
             "CTLA4|CTLA4",
             "CTNNA1|alpha-Catenin",
             "CTNNB1|beta-Catenin",
             "DIABLO|Smac",
             "DIRAS3|DIRAS3",
             "DPP4|CD26",
             "DUSP4|DUSP4",
             "DVL3|Dvl3",
             "E2F1|E2F1",
             "EEF2K|eEF2K",
             "EEF2|eEF2",
             "EGFR|EGFR",
             "EGFR|EGFR_pY1068",
             "EGFR|EGFR_pY1173",
             "EIF4EBP1|4E-BP1",
             "EIF4EBP1|4E-BP1_pS65",
             "EIF4EBP1|4E-BP1_pT37_T46",
             "EIF4EBP1|4E-BP1_pT70",
             "EIF4E|eIF4E",
             "EIF4G1|eIF4G",
             "ENY2|ENY2",
             "EPPK1|EPPK1",
             "ERBB2|HER2",
             "ERBB2|HER2_pY1248",
             "ERBB3|HER3",
             "ERBB3|HER3_pY1289",
             "ERCC1|ERCC1",
             "ERCC5|ERCC5",
             "ERRFI1|MIG-6",
             "ESR1|ER-alpha",
             "ESR1|ER-alpha_pS118",
             "ETS1|ETS-1",
             "EZH2|EZH2",
             "FASN|FASN",
             "FN1|Fibronectin",
             "FOXM1|FoxM1",
             "FOXO3|FOXO3a",
             "FOXO3|FOXO3a_pS318_S321",
             "G6PD|G6PD",
             "GAB2|GAB2",
             "GAPDH|GAPDH",
             "GATA3|GATA3",
             "GATA6|GATA6",
             "GSK3A GSK3B|GSK3-alpha-beta",
             "GSK3A GSK3B|GSK3-alpha-beta_pS21_S9",
             "GSK3A GSK3B|GSK3_pS9",
             "GYG1|GYG-Glycogenin1",
             "GYS1|GYS",
             "GYS1|GYS_pS641",
             "HIF1A|HIF-1_alpha",
             "HSPA1A|HSP70",
             "IGFBP2|IGFBP2",
             "IGF1R|IGF1R_pY1135_Y1136",
             "INPP4B|INPP4B",
             "IRF1|IRF-1",
             "IRS1|IRS1",
             "ITGA2|CD49b",
             "JAK2|Jak2",
             "JUN|c-Jun_pS73",
             "KAT2A|GCN5L2",
             "KDR|VEGFR2",
             "KEAP1|KEAP1",
             "KIT|c-Kit",
             "KRT5|CK5",
             "LCK|Lck",
             "LCN2|LCN2a",
             "LDHA|LDHA",
             "LDHB|LDHB",
             "MACC1|MACC1",
             "MAP2K1|MEK1",
             "MAP2K1|MEK1_pS217_S221",
             "MAPK1 MAPK3|MAPK_pT202_Y204",
             "MAPK14|p38_MAPK",
             "MAPK14|p38_pT180_Y182",
             "MAPK1|ERK2",
             "MAPK8|JNK_pT183_pY185",
             "MAPK9|JNK2",
             "MET|c-Met",
             "MET|c-Met_pY1235",
             "MRE11|Mre11",
             "MS4A1|CD20",
             "MSH2|MSH2",
             "MSH6|MSH6",
             "MTCO2|Mitochondria",
             "MTOR|mTOR",
             "MTOR|mTOR_pS2448",
             "MYC|c-Myc",
             "MYH11|MYH11",
             "MYH9|Myosin-IIa",
             "MYH9|Myosin-IIa_pS1943",
             "NAPSA|Napsin-A",
             "NDRG1|NDRG1_pT346",
             "NF2|NF2",
             "NFE2L2|Nrf2",
             "NFKB1|NF-kB-p65_pS536",
             "NKX2-1|TTF1",
             "NOTCH1|Notch1",
             "NRAS|N-Ras",
             "NRG1|Heregulin",
             "PARK7|DJ-1",
             "PARP1|PARP-Ab-3",
             "PARP1|PARP1",
             "PARP1|PARP_cleaved",
             "PCNA|PCNA",
             "PDCD1|PDCD1",
             "PDCD4|PDCD4",
             "PDK1|PDK1",
             "PDK1|PDK1_pS241",
             "PEA15|PEA15",
             "PEA15|PEA15_pS116",
             "PECAM1|CD31",
             "PGR|PR",
             "PIK3CA|PI3K-p110-alpha",
             "PIK3R1|PI3K-p85",
             "PKM2|PKM2",
             "PRDX1|PRDX1",
             "PREX1|PREX1",
             "PRKAA1|AMPK_alpha",
             "PRKAA1|AMPK_pT172",
             "PRKCA|PKC-alpha",
             "PRKCA|PKC-alpha_pS657",
             "PRKCB|PKC-pan_BetaII_pS660",
             "PRKCD|PKC-delta_pS664",
             "PTEN|PTEN",
             "PTPN11|SHP-2_pY542",
             "PXN|Paxillin",
             "PYGB|PYGB",
             "PYGB|PYGB-AB2",
             "PYGL|PYGL",
             "PYGM|PYGM",
             "RAB11A RAB11B|Rab11",
             "RAB25|Rab25",
             "RAD50|Rad50",
             "RAD51|Rad51",
             "RAF1|C-Raf",
             "RAF1|C-Raf_pS338",
             "RB1|Rb",
             "RB1|Rb_pS807_S811",
             "RBM15|RBM15",
             "RET|Ret_pY905",
             "RICTOR|Rictor",
             "RICTOR|Rictor_pT1135",
             "RPS6KA1|p90RSK",
             "RPS6KA1|p90RSK_pT359_S363",
             "RPS6KB1|p70S6K",
             "RPS6KB1|p70S6K_pT389",
             "RPS6|S6",
             "RPS6|S6_pS235_S236",
             "RPS6|S6_pS240_S244",
             "RPTOR|Raptor",
             "SCD|SCD",
             "SDHB|Complex-II_subunit30",
             "SERPINE1|PAI-1",
             "SETD2|SETD2",
             "SHC1|Shc_pY317",
             "SLC1A5|SLC1A5",
             "SMAD1|Smad1",
             "SMAD3|Smad3",
             "SMAD4|Smad4",
             "SNAI1|Snail",
             "SQSTM1|p62-LCK-ligand",
             "SRC|Src",
             "SRC|Src_pY416",
             "SRC|Src_pY527",
             "SRSF1|SF2",
             "STAT3|STAT3_pY705",
             "STAT5A|STAT5-alpha",
             "STK11|LKB1",
             "STMN1|Stathmin",
             "SYK|Syk",
             "SYP|Synaptophysin",
             "TFRC|TFRC",
             "TGM2|Transglutaminase",
             "TIGAR|TIGAR",
             "TP53BP1|53BP1",
             "TP53|p53",
             "TP63|p63",
             "TSC1|TSC1",
             "TSC2|Tuberin",
             "TSC2|Tuberin_pT1462",
             "TUBA1B|Acetyl-a-Tubulin-Lys40",
             "TYMS|Thymidilate-Synthase",
             "WWTR1|TAZ",
             "XBP1|XBP1",
             "XRCC1|XRCC1",
             "XRCC5|Ku80",
             "YAP1|YAP",
             "YAP1|YAP_pS127",
             "YBX1|YB-1",
             "YBX1|YB-1_pS102",
             "YWHAB|14-3-3_beta",
             "YWHAE|14-3-3_epsilon",
             "YWHAZ|14-3-3_zeta" )
    names(sPr) <- sapply(strsplit(sPr, split = "\\|"), function(x){x[2]})
    sColNames <- c("Composite Element REF", "Protein Expression")
    ld  <- ColumnsFromFiles(sFileNameById,
                            sColNames,
                            sSortByCol = "Composite Element REF",
                            sSkipLines = 1,
                            sNa = "NA")
    lnAbValue <- lapply(ld,
                        function(x){
                            ProbeValue(x,
                                       sColProbe = "Composite Element REF",
                                       sColValue = "Protein Expression",
                                       nStrip = 0)})
    sAb <- sort(unique(unlist(lapply(lnAbValue, names))))
    sPr4Ab <- sPr[StripEnd(sUnstripped = sAb, nStrip = 4,
                           sStripEnd = "right")]
    sPr2Ab <- paste(sapply(strsplit(sPr4Ab, split = "\\|"),
                           function(x){x[1]}),
                    sAb, sep = "|")
    names(sPr2Ab) <- sAb
    m <- matrix(nrow = length(sPr2Ab), ncol = length(sFileNameById))
    rownames(m) <- sPr2Ab
    colnames(m) <- names(sFileNameById)
    for (s in names(sFileNameById)) {
        m[, s] <- lnAbValue[[s]][sAb]
    }
    dMerged <- cbind(protein = sPr2Ab, m)
    return(dMerged)
}

#' Main function of \code{Merge}, distribute filenames by assay platform
#'
#' @param sFileNameById Character vector of filename (named with file_id).
#' @param s1Assay String of assay platform, used in \code{Filter}.
#' @return A \code{data.frame} of merged table.
Merge <- function(sFileNameById,
                  s1Assay) {
    print("merging files: merging unzipped data files ...")
    if (is.null(sFileNameById)) {
        print("merging files: no file satisfies the lFilter!")
        dMerged <- NULL
        return(dMerged)
    } else if (s1Assay %in% c("cna_cnv.hg18",
                              "cna_cnv.hg19",
                              "cna_nocnv.hg18",
                              "cna_nocnv.hg19")) {
        dMerged <- CheckId(MergeCopy(sFileNameById), s1Assay)
    } else if (s1Assay %in% c("exonJunction_RNAseq")) {
        dMerged <- CheckId(MergeExonJunction(sFileNameById), s1Assay)
    } else if (s1Assay %in% c("exon_RNAseq")) {
        dMerged <- CheckId(MergeExon(sFileNameById), s1Assay)
    } else if (s1Assay %in% c("gene_Array")) {
        dMerged <- CheckId(MergeGeneArray(sFileNameById), s1Assay)
    } else if (s1Assay %in% c("gene.normalized_RNAseq")) {
        dMerged <- CheckId(MergeGeneRnaSeqNorm(sFileNameById), s1Assay)
    } else if (s1Assay %in% c("gene_RNAseq")) {
        dMerged <- CheckId(MergeGeneRnaSeqUnnorm(sFileNameById), s1Assay)
    } else if (s1Assay %in% c("isoform.normalized_RNAseq")) {
        dMerged <- CheckId(MergeGeneIsoRnaSeqNorm(sFileNameById), s1Assay)
    } else if (s1Assay %in% c("isoform_RNAseq")) {
        dMerged <- CheckId(MergeGeneIsoRnaSeqUnnorm(sFileNameById), s1Assay)
    } else if (s1Assay %in% c("methylation_27",
                              "methylation_450")) {
        dMerged <- CheckId(MergeMethy(sFileNameById), s1Assay)
    } else if (s1Assay %in% c("mir_GA.hg18",
                              "mir_GA.hg19",
                              "mir_GA.hg19.mirbase20",
                              "mir_HiSeq.hg18",
                              "mir_HiSeq.hg19",
                              "mir_HiSeq.hg19.mirbase20")) {
        dMerged <- CheckId(MergeMir(sFileNameById), s1Assay)
    } else if (s1Assay %in% c("mirIsoform_GA.hg18",
                              "mirIsoform_GA.hg19",
                              "mirIsoform_GA.hg19.mirbase20",
                              "mirIsoform_HiSeq.hg18",
                              "mirIsoform_HiSeq.hg19",
                              "mirIsoform_HiSeq.hg19.mirbase20")) {
        dMerged <- CheckId(MergeMirIso(sFileNameById), s1Assay)
    } else if (s1Assay %in% c("protein_RPPA")) {
        dMerged <- CheckId(MergeProtein(sFileNameById), s1Assay)
    }
    print("merging files: merging unzipped data files done!")
    return(dMerged)
}


#  =============================================================================
#  adapting functions, internal <> interface, NOT used directly by user
#  =============================================================================

#' Check the user specified parameters
#'
#' @param sCancer String of cancer type.
#' @param sAssay Character vector of assay platform.
#' @param sampleTypeName Character vector of name for sample_type_id.
#' @param sAssayGroup String of assay platform goup.
#' @return List of checked parameters.
CheckParam <- function(sCancer,
                       sAssay,
                       sampleTypeName,
                       sAssayGroup) {
    sCancerAll <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC",
                    "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML",
                    "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD",
                    "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT",
                    "THCA", "THYM", "UCEC", "UCS", "UVM")
    sSampleTypeIdAll <- sprintf("%02d", c(seq(14), 20, 40, 50, 60, 61))
    names(sSampleTypeIdAll) <-
        c("TP",   # 01, 'Primary Tumor'
          "TR",   # 02, 'Recurrent Tumor'
          "TB",   # 03, 'Primary Blood Derived Cancer - Peripheral Blood'
          "TRBM", # 04, 'Recurrent Blood Derived Cancer - Bone Marrow'
          "TAP",  # 05, 'Additional - New Primary'
          "TM",   # 06, 'Metastatic'
          "TAM",  # 07, 'Additional Metastatic'
          "THOC", # 08, 'Human Tumor Original Cells'
          "TBM",  # 09, 'Primary Blood Derived Cancer - Bone Marrow'
          "NB",   # 10, 'Blood Derived Normal'
          "NT",   # 11, 'Solid Tissue Normal'
          "NBC",  # 12, 'Buccal Cell Normal'
          "NEBV", # 13, 'EBV Immortalized Normal'
          "NBM",  # 14, 'Bone Marrow Normal'
          "CELLC",# 20, 'Control Analyte'
          "TRB",  # 40, 'Recurrent Blood Derived Cancer - Peripheral Blood'
          "CELL", # 50, 'Cell Lines'
          "XP",   # 60, 'Primary Xenograft Tissue'
          "XCL")  # 61, 'Cell Line Derived Xenograft Tissue'
    lsAssayGroup <- list(cna = c("cna_cnv.hg18",
                                 "cna_cnv.hg19",
                                 "cna_nocnv.hg18",
                                 "cna_nocnv.hg19"),
                         gene = c("gene_Array",
                                  "gene.normalized_RNAseq",
                                  "gene_RNAseq",
                                  "isoform.normalized_RNAseq",
                                  "isoform_RNAseq",
                                  "exon_RNAseq",
                                  "exonJunction_RNAseq"),
                         methy = c("methylation_27",
                                   "methylation_450"),
                         mir = c("mir_GA.hg18",
                                 "mir_GA.hg19",
                                 "mir_GA.hg19.mirbase20",
                                 "mir_HiSeq.hg18",
                                 "mir_HiSeq.hg19",
                                 "mir_HiSeq.hg19.mirbase20"),
                         mirIsoform = c("mirIsoform_GA.hg18",
                                        "mirIsoform_GA.hg19",
                                        "mirIsoform_GA.hg19.mirbase20",
                                        "mirIsoform_HiSeq.hg18",
                                        "mirIsoform_HiSeq.hg19",
                                        "mirIsoform_HiSeq.hg19.mirbase20"),
                         protein = c("protein_RPPA"),
                         somatic = c("somaticMutation_DNAseq"),
                         itraq = c("glycoproteome_iTRAQ",
                                   "phosphoproteome_iTRAQ",
                                   "proteome_iTRAQ"))
    sAssaySub <- lsAssayGroup[[sAssayGroup]]
    if (is.null(sCancer)) {
        sCancer <- sCancerAll
    } else if (!all(sCancer %in% sCancerAll)) {
        print(c("cancerType should be 'NULL' (all) or one of: ", sCancerAll))
        stopifnot(!all(sCancer %in% sCancerAll))
    }
    if (is.null(sAssay)) {
        sAssay <- sAssaySub
    } else if (!all(sAssay %in% sAssaySub)) {
        print(c("assayPlatform should be 'NULL' (all) or one of: ", sAssaySub))
        stopifnot(!all(sAssay %in% sAssaySub))
    }
    if (is.null(sampleTypeName)) {
        sSampleTypeId <- sSampleTypeIdAll
    } else if (!all(sampleTypeName %in% names(sSampleTypeIdAll))) {
        print(paste("tissueType should be 'NULL' (all) or one of:",
                    paste(sSampleTypeIdAll, collapse = ","),
                    "01) TP = 'Primary Tumor';",
                    "02) TR = 'Recurrent Tumor';",
                    "03) TB = 'Primary Blood Derived Cancer - Peripheral Blood';",
                    "04) TRBM = 'Recurrent Blood Derived Cancer - Bone Marrow';",
                    "05) TAP = 'Additional - New Primary';",
                    "06) TM = 'Metastatic';",
                    "07) TAM = 'Additional Metastatic';",
                    "08) THOC = 'Human Tumor Original Cells';",
                    "09) TBM = 'Primary Blood Derived Cancer - Bone Marrow';",
                    "10) NB = 'Blood Derived Normal';",
                    "11) NT = 'Solid Tissue Normal';",
                    "12) NBC = 'Buccal Cell Normal';",
                    "13) NEBV = 'EBV Immortalized Normal';",
                    "14) NBM = 'Bone Marrow Normal';",
                    "20) CELLC = 'Control Analyte';",
                    "40) TRB = 'Recurrent Blood Derived Cancer - Peripheral Blood';",
                    "50) CELL = 'Cell Lines';",
                    "60) XP = 'Primary Xenograft Tissue';",
                    "61) XCL = 'Cell Line Derived Xenograft Tissue'.",
                    sep = " "))
        stopifnot(all(sampleTypeName %in% names(sSampleTypeIdAll)))
    } else {
        sSampleTypeId <- sSampleTypeIdAll[sampleTypeName]
    }
    return(list(sCancer = sCancer, sAssay = sAssay,
                sSampleTypeId = sSampleTypeId))
}

#' Pipe of metadata, download and merge (general function)
#'
#' @param sCancer String of cancer type.
#' @param s1Assay String of assay platform, used in \code{Filter}.
#' @param sSampleTypeId Character vector of sample_type_id: "01", ..., "14", etc.
#' @param sBarcode Character vector of barcode, to specify patients.
#' @param sDir String of directory for temporary files.
#' @param sArchive String of archive type: "legacy" or "".
#' @param sFieldMeta Character vector of colnames in metadata.
#' @param nEntityCount Number of entity (row) in the metadata: "-1" means all.
#' @param sEndpoint String of endpoint: "files".
#' @return A \code{data.frame} of merged data.
Pipe <- function(sCancer,
                 s1Assay,
                 sSampleTypeId,
                 sBarcode = NULL,
                 sDir = ".",
                 sArchive = "legacy",
                 sFieldMeta = "",
                 nEntityCount = (-1),
                 sEndpoint = "files") {
    sTissueType <- c("TP", "TR", "TB", "TRBM", "TAP", "TM", "TAM", "THOC", "TBM",
                     "NB", "NT", "NBC", "NEBV", "NBM", "CELLC", "TRB", "CELL", "XP", "XCL")
    names(sTissueType) <- sprintf("%02d", c(seq(14), 20, 40, 50, 60, 61))
    dMeta <- MetaData(sCancer = sCancer,
                      s1Assay = s1Assay,
                      sSampleTypeId = sSampleTypeId,
                      sDir = sDir,
                      sArchive = sArchive,
                      sFieldMeta = sFieldMeta,
                      nEntityCount = nEntityCount,
                      sEndpoint = sEndpoint)
    if (!is.null(dMeta)) {
        sFileId2Bar <- rownames(dMeta)
        names(sFileId2Bar) <- dMeta$file_id
        if (!is.null(sBarcode)) {
            sFileId2Bar <- sFileId2Bar[ifelse(substr(sFileId2Bar, 1, 12) %in%
                                              substr(sBarcode, 1, 12), TRUE, FALSE)]
        }
        if (length(sFileId2Bar) > 0) {
            sFileNameById <- FileNameById(sFileId2Bar,
                                          sDir = sDir,
                                          sArchive = sArchive)
            dPiped <- Merge(sFileNameById, s1Assay)
        } else {
            dPiped <- NULL
        }
    } else {
        # print(paste("metadata = NULL, when cancerType = ",
        #             paste(sCancer, collapse = "|"),
        #             " & assayPlatform = ", s1Assay, " & tissueType = ",
        #             paste(sTissueType[sSampleTypeId], collapse = "|"),
        #             sep = ""))
        dPiped <- NULL
    }
    return(dPiped)
}

#' Pipe of metadata, download and merge (with batch download)
#'
#' @param sCancer String of cancer type.
#' @param s1Assay String of assay platform, used in \code{Filter}.
#' @param sSampleTypeId Character vector of sample_type_id: "01", ..., "14", etc.
#' @param sBarcode Character vector of barcode, to specify patients.
#' @param sDir String of directory for temporary files.
#' @param sArchive String of archive type: "legacy" or "".
#' @param sFieldMeta Character vector of colnames in metadata.
#' @param nEntityCount Number of entity (row) in the metadata: "-1" means all.
#' @param sEndpoint String of endpoint: "files".
#' @return A \code{data.frame} of merged data.
PipeBatch <- function(sCancer,
                      s1Assay,
                      sSampleTypeId,
                      sBarcode = NULL,
                      sDir = ".",
                      sArchive = "legacy",
                      sFieldMeta = "",
                      nEntityCount = (-1),
                      sEndpoint = "files") {
    sTissueType <- c("TP", "TR", "TB", "TRBM", "TAP", "TM", "TAM", "THOC", "TBM",
                     "NB", "NT", "NBC", "NEBV", "NBM", "CELLC", "TRB", "CELL", "XP", "XCL")
    names(sTissueType) <- sprintf("%02d", c(seq(14), 20, 40, 50, 60, 61))
    dMeta <- MetaData(sCancer = sCancer,
                      s1Assay = s1Assay,
                      sSampleTypeId = sSampleTypeId,
                      sDir = sDir,
                      sArchive = sArchive,
                      sFieldMeta = sFieldMeta,
                      nEntityCount = nEntityCount,
                      sEndpoint = sEndpoint)
    if (!is.null(dMeta)) {
        sFileId2Bar <- rownames(dMeta)
        names(sFileId2Bar) <- dMeta$file_id
        if (!is.null(sBarcode)) {
            sFileId2Bar <- sFileId2Bar[ifelse(substr(sFileId2Bar, 1, 12) %in%
                                              substr(sBarcode, 1, 12), TRUE, FALSE)]
        }
        if (length(sFileId2Bar) > 0) {
            sFileNameById <- vector()  # batch download
            n1Batch <- 50
            for (n in seq(ceiling(length(sFileId2Bar) / n1Batch)) - 1) {
                nStart <- n * n1Batch + 1
                nStop <- ifelse(n == (ceiling(length(sFileId2Bar)/n1Batch) - 1),
                                length(sFileId2Bar),
                                (n + 1) * n1Batch)
                sBatch <- FileNameById(sFileId2Bar[nStart : nStop],
                                       sDir = sDir,
                                       sArchive = sArchive)
                sFileNameById <- c(sFileNameById, sBatch)
            }  # batch download
            dPiped <- Merge(sFileNameById, s1Assay)
        } else {
            dPiped <- NULL
        }
    } else {
        # print(paste("metadata = NULL, when cancerType = ",
        #             paste(sCancer, collapse = "|"),
        #             " & assayPlatform = ", s1Assay, " & tissueType = ",
        #             paste(sTissueType[sSampleTypeId], collapse = "|"),
        #             sep = ""))
        dPiped <- NULL
    }
    return(dPiped)
}

#' Pipe of metadata, download and merge (for mir data with 705 probes)
#'
#' @param sCancer String of cancer type.
#' @param s1Assay String of assay platform, used in \code{Filter}.
#' @param sSampleTypeId Character vector of sample_type_id: "01", ..., "14", etc.
#' @param sBarcode Character vector of barcode, to specify patients.
#' @param sDir String of directory for temporary files.
#' @param sArchive String of archive type: "legacy" or "".
#' @param sFieldMeta Character vector of colnames in metadata.
#' @param nEntityCount Number of entity (row) in the metadata: "-1" means all.
#' @param sEndpoint String of endpoint: "files".
#' @return A \code{data.frame} of merged data.
PipeMirLt23k <- function(sCancer,
                         s1Assay,
                         sSampleTypeId,
                         sBarcode = NULL,
                         sDir = ".",
                         sArchive = "legacy",
                         sFieldMeta = "",
                         nEntityCount = (-1),
                         sEndpoint = "files") {
    sTissueType <- c("TP", "TR", "TB", "TRBM", "TAP", "TM", "TAM", "THOC", "TBM",
                     "NB", "NT", "NBC", "NEBV", "NBM", "CELLC", "TRB", "CELL", "XP", "XCL")
    names(sTissueType) <- sprintf("%02d", c(seq(14), 20, 40, 50, 60, 61))
    dMeta <- MetaData(sCancer = sCancer,
                      s1Assay = s1Assay,
                      sSampleTypeId = sSampleTypeId,
                      sDir = sDir,
                      sArchive = sArchive,
                      sFieldMeta = sFieldMeta,
                      nEntityCount = nEntityCount,
                      sEndpoint = sEndpoint)
    if (!is.null(dMeta)) {
        sFileId2Bar <- rownames(dMeta)
        names(sFileId2Bar) <- dMeta$file_id
        nLt23k <- as.numeric(dMeta$file_size) < 23000  # file_size < 23000
        sFileId2Bar <- sFileId2Bar[nLt23k]
        if (!is.null(sBarcode)) {
            sFileId2Bar <- sFileId2Bar[ifelse(substr(sFileId2Bar, 1, 12) %in%
                                              substr(sBarcode, 1, 12), TRUE, FALSE)]
        }
        if (length(sFileId2Bar) > 0) {
            sFileNameById <- FileNameById(sFileId2Bar,
                                          sDir = sDir,
                                          sArchive = sArchive)
            dPiped <- Merge(sFileNameById, s1Assay)
        } else {
            dPiped <- NULL
        }
    } else {
        # print(paste("metadata = NULL, when cancerType = ",
        #             paste(sCancer, collapse = "|"),
        #             " & assayPlatform = ", s1Assay, " & tissueType = ",
        #             paste(sTissueType[sSampleTypeId], collapse = "|"),
        #             sep = ""))
        dPiped <- NULL
    }
    return(dPiped)
}

#' Pipe of metadata, download and merge (for mir data with more than 705 probes)
#'
#' @param sCancer String of cancer type.
#' @param s1Assay String of assay platform, used in \code{Filter}.
#' @param sSampleTypeId Character vector of sample_type_id: "01", ..., "14", etc.
#' @param sBarcode Character vector of barcode, to specify patients.
#' @param sDir String of directory for temporary files.
#' @param sArchive String of archive type: "legacy" or "".
#' @param sFieldMeta Character vector of colnames in metadata.
#' @param nEntityCount Number of entity (row) in the metadata: "-1" means all.
#' @param sEndpoint String of endpoint: "files".
#' @return A \code{data.frame} of merged data.
PipeMirGt23k <- function(sCancer,
                         s1Assay,
                         sSampleTypeId,
                         sBarcode = NULL,
                         sDir = ".",
                         sArchive = "legacy",
                         sFieldMeta = "",
                         nEntityCount = (-1),
                         sEndpoint = "files") {
    sTissueType <- c("TP", "TR", "TB", "TRBM", "TAP", "TM", "TAM", "THOC", "TBM",
                     "NB", "NT", "NBC", "NEBV", "NBM", "CELLC", "TRB", "CELL", "XP", "XCL")
    names(sTissueType) <- sprintf("%02d", c(seq(14), 20, 40, 50, 60, 61))
    dMeta <- MetaData(sCancer = sCancer,
                      s1Assay = s1Assay,
                      sSampleTypeId = sSampleTypeId,
                      sDir = sDir,
                      sArchive = sArchive,
                      sFieldMeta = sFieldMeta,
                      nEntityCount = nEntityCount,
                      sEndpoint = sEndpoint)
    if (!is.null(dMeta)) {
        sFileId2Bar <- rownames(dMeta)
        names(sFileId2Bar) <- dMeta$file_id
        nGt23k <- as.numeric(dMeta$file_size) > 23000  # file_size > 23000
        sFileId2Bar <- sFileId2Bar[nGt23k]
        if (!is.null(sBarcode)) {
            sFileId2Bar <- sFileId2Bar[ifelse(substr(sFileId2Bar, 1, 12) %in%
                                              substr(sBarcode, 1, 12), TRUE, FALSE)]
        }
        if (length(sFileId2Bar) > 0) {
            sFileNameById <- FileNameById(sFileId2Bar,
                                          sDir = sDir,
                                          sArchive = sArchive)
            dPiped <- Merge(sFileNameById, s1Assay)
        } else {
            dPiped <- NULL
        }
    } else {
        # print(paste("metadata = NULL, when cancerType = ",
        #             paste(sCancer, collapse = "|"),
        #             " & assayPlatform = ", s1Assay, " & tissueType = ",
        #             paste(sTissueType[sSampleTypeId], collapse = "|"),
        #             sep = ""))
        dPiped <- NULL
    }
    return(dPiped)
}

#' Pipe of metadata, download and merge (for somatic mutation data)
#'
#' @param sCancer String of cancer type
#' @param s1Assay String of assay platform, used in \code{Filter}.
#' @param sSampleTypeId Character vector of sample_type_id: "01", ..., "14", etc.
#' @param sBarcode Character vector of barcode, to specify patients.
#' @param sDir String of directory for temporary files.
#' @param sArchive String of archive type: "legacy" or "".
#' @param sFieldMeta Character vector of colnames in metadata.
#' @param nEntityCount Number of entity (row) in the metadata: "-1" means all.
#' @param sEndpoint String of endpoint: "files".
#' @return A \code{data.frame} of merged data.
PipeSomatic <- function(sCancer,
                        s1Assay,
                        sSampleTypeId,
                        sBarcode = NULL,
                        sDir = ".",
                        sArchive = "legacy",
                        sFieldMeta = "",
                        nEntityCount = (-1),
                        sEndpoint = "files") {
    sTissueType <- c("TP", "TR", "TB", "TRBM", "TAP", "TM", "TAM", "THOC", "TBM",
                     "NB", "NT", "NBC", "NEBV", "NBM", "CELLC", "TRB", "CELL", "XP", "XCL")
    names(sTissueType) <- sprintf("%02d", c(seq(14), 20, 40, 50, 60, 61))
    dMeta <- MetaDataSoma(sCancer = sCancer,
                          s1Assay = s1Assay,
                          sSampleTypeId = sSampleTypeId,
                          sDir = sDir,
                          sArchive = sArchive,
                          sFieldMeta = sFieldMeta,
                          nEntityCount = nEntityCount,
                          sEndpoint = sEndpoint)
    sFileId2Bar <- dMeta$file_name
    names(sFileId2Bar) <- dMeta$file_id
    sFileNameById <- FileNameById(sFileId2Bar,
                                  sDir = sDir,
                                  sArchive = sArchive)
    sColNames <- c("hugo_symbol",
                   "entrez_gene_id",
                   "center",
                   "ncbi_build",
                   "chromosome",
                   "start_position",
                   "end_position",
                   "strand",
                   "variant_classification",
                   "variant_type",
                   "reference_allele",
                   "tumor_seq_allele1",
                   "tumor_seq_allele2",
                   "dbsnp_rs",
                   "dbsnp_val_status",
                   "tumor_sample_barcode",
                   "matched_norm_sample_barcode",
                   "match_norm_seq_allele1",
                   "match_norm_seq_allele2",
                   "tumor_validation_allele1",
                   "tumor_validation_allele2",
                   "match_norm_validation_allele1",
                   "match_norm_validation_allele2",
                   "verification_status",
                   "validation_status",
                   "mutation_status",
                   "sequencing_phase",
                   "sequence_source",
                   "validation_method",
                   "score",
                   "bam_file",
                   "sequencer",
                   "tumor_sample_uuid",
                   "matched_norm_sample_uuid")
    names(sColNames) <- c("Hugo_Symbol",
                          "Entrez_Gene_Id",
                          "Center",
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
                          "Tumor_Sample_Barcode",
                          "Matched_Norm_Sample_Barcode",
                          "Match_Norm_Seq_Allele1",
                          "Match_Norm_Seq_Allele2",
                          "Tumor_Validation_Allele1",
                          "Tumor_Validation_Allele2",
                          "Match_Norm_Validation_Allele1",
                          "Match_Norm_Validation_Allele2",
                          "Verification_Status",
                          "Validation_Status",
                          "Mutation_Status",
                          "Sequencing_Phase",
                          "Sequence_Source",
                          "Validation_Method",
                          "Score",
                          "BAM_File",
                          "Sequencer",
                          "Tumor_Sample_UUID",
                          "Matched_Norm_Sample_UUID")
    ldPiped <- list()
    for (sMaf in sFileNameById) {
        dMaf <- read.csv(sMaf,
                         sep = "\t",
                         row.names = NULL,
                         as.is = TRUE,
                         na.strings = "",
                         comment.char = "#",
                         check.names = FALSE)
        if (!is.null(dMaf)) {
            colnames(dMaf) <- tolower(colnames(dMaf))
            dPiped <- cbind(dMaf[, sColNames])
            colnames(dPiped) <- names(sColNames)
            if (!is.null(sBarcode)) {
                bBar <- ifelse(substr(dPiped[, "Tumor_Sample_Barcode"], 1, 12) %in%
                               substr(sBarcode, 1, 12), TRUE, FALSE)
                if (any(bBar)) {
                    dPiped <- dPiped[bBar, , drop = FALSE]
                } else {
                    dPiped <- NULL
                }
            }
            if (!is.null(dPiped) & length(sSampleTypeId) < 14) {
                bSampleTypeId <- ifelse(substr(dPiped[, "Tumor_Sample_Barcode"],
                                               14, 15) %in% sSampleTypeId, TRUE, FALSE)
                if (any(bSampleTypeId)) {
                    dPiped <- dPiped[bSampleTypeId, , drop = FALSE]
                } else {
                    dPiped <- NULL
                }
            }
            if (!is.null(dPiped)) {
                bSymbol <- ifelse(dPiped[, "Hugo_Symbol"] %in% c("."), FALSE, TRUE)
                if (any(bSymbol)) {
                    dPiped <- dPiped[bSymbol, , drop = FALSE]
                }
            }
        } else {
            # print(paste("metadata = NULL, when cancerType = ",
            #             paste(sCancer, collapse = "|"),
            #             " & assayPlatform = ", s1Assay, " & tissueType = ",
            #             paste(sTissueType[sSampleTypeId], collapse = "|"),
            #             sep = ""))
            dPiped <- NULL
        }
        ldPiped[[sMaf]] <- CheckId(dPiped, s1Assay = "somaticMutation_DNAseq")
    }
    return(ldPiped)
}


#  =============================================================================
#  interface functions, directly used by user
#  =============================================================================

#' DownloadmiRNASeqData: get miRNASeq data, assayPlatform %in% c("mir_GA.hg18", "mir_GA.hg19", "mir_GA.hg19.mirbase20", "mir_HiSeq.hg18", "mir_HiSeq.hg19", "mir_HiSeq.hg19.mirbase20")
#'
#' @param cancerType (i.e. sCancer, vector of cancer type), length(sCancer)> = 1
#' @param assayPlatform (i.e. sAssay, vector of type), length(sAssay)> = 1, assayPlatform %in% c("mir_GA.hg18", "mir_GA.hg19", "mir_GA.hg19.mirbase20", "mir_HiSeq.hg18", "mir_HiSeq.hg19", "mir_HiSeq.hg19.mirbase20")
#' @param tissueType (i.e. sampleTypeName, vector of sample_type_name, could be transfered to sample_type_id): c("TP", "TR", ...) -> c(01, 02, ...)
#' @param saveFilePath (string of path to save the merged data): absolute or relative path
#' @param saveFileName (string of filename prefix)
#' @param inputPatientIDs (i.e. sBarcode, vector of barcode), find specified patients
#' @return sFileName (vector of path/filename), could be used by module B
DownloadmiRNASeqData <- function(cancerType = NULL,
                                 assayPlatform = NULL,
                                 tissueType = NULL,
                                 saveFilePath = ".",
                                 saveFileName = "",
                                 inputPatientIDs = NULL) {
    sArchive <- "legacy"; sFieldMeta <- ""; nEntityCount <- (-1); sEndpoint <- "files"
    options(warn = -1)
    if (saveFilePath != ".") {dir.create(saveFilePath, recursive = TRUE)}
    sDir <- paste("tmp", TimeNow(), sep = "_"); dir.create(sDir)
    #  if (saveFilePath != ".") {
    #    if (!dir.exists(saveFilePath)) {
    #      dir.create(saveFilePath, recursive = TRUE)
    #    }
    #  } # dir.exists # R >= 3.2
    #  sDir <- paste("tmp", TimeNow(), sep = "_")
    #  if (!dir.exists(sDir)) {dir.create(sDir)} # dir.exists # R >= 3.2
    l <- CheckParam(sCancer = cancerType,
                    sAssay = assayPlatform,
                    sampleTypeName = tissueType,
                    sAssayGroup = "mir")
    sFileName <- rep(NA, length(l$sAssay))
    names(sFileName) <- l$sAssay
    for (s1Assay in l$sAssay) {
        dLt23k <- PipeMirLt23k(sCancer = l$sCancer,
                               s1Assay = s1Assay,
                               sSampleTypeId = l$sSampleTypeId,
                               sBarcode = inputPatientIDs,
                               sDir = sDir,
                               sArchive = sArchive,
                               sFieldMeta = sFieldMeta,
                               nEntityCount = nEntityCount,
                               sEndpoint = sEndpoint)
        dGt23k <- PipeMirGt23k(sCancer = l$sCancer,
                               s1Assay = s1Assay,
                               sSampleTypeId = l$sSampleTypeId,
                               sBarcode = inputPatientIDs,
                               sDir = sDir,
                               sArchive = sArchive,
                               sFieldMeta = sFieldMeta,
                               nEntityCount = nEntityCount,
                               sEndpoint = sEndpoint)
        if (!is.null(dLt23k)) {
            sFileName705 <- paste(saveFilePath,
                                  ifelse(endsWith(saveFilePath, "/"), "", "/"),
                                  ifelse(saveFileName == "", "",
                                         paste(saveFileName, "__", sep = "")),
                                  paste(l$sCancer, collapse = "_"),
                                  "__",
                                  s1Assay,
                                  "__705mir__",
                                  ifelse(is.null(tissueType), "tissueTypeAll",
                                         paste(tissueType, collapse = "_")),
                                  "__",
                                  TimeNow(),
                                  sep = "")
            dPiped <- dLt23k
            save(dLt23k, file = paste(sFileName705, ".rda", sep = ""))
            #write.table(dLt23k,
            #            file = paste(sFileName705, ".txt", sep = ""),
            #            quote = FALSE,
            #            sep = "\t",
            #            col.names = TRUE,
            #            row.names = FALSE,
            #            na = "")
            rm(dLt23k); gc()  # clear the memory
            sFileName[paste(s1Assay, "_705", sep = "")] <- paste(sFileName705, ".txt", sep = "")
        }
        if (!is.null(dGt23k)) {
            sFilename <- paste(saveFilePath,
                               ifelse(endsWith(saveFilePath, "/"), "", "/"),
                               ifelse(saveFileName == "", "",
                                      paste(saveFileName,
                                            "__", sep = "")),
                               paste(l$sCancer, collapse = "_"),
                               "__",
                               s1Assay,
                               "__",
                               ifelse(is.null(tissueType), "tissueTypeAll",
                                      paste(tissueType, collapse = "_")),
                               "__",
                               TimeNow(),
                               sep = "")
            dPiped <- dGt23k
            save(dPiped, file = paste(sFilename, ".rda", sep = ""))
            #write.table(dGt23k,
            #            file = paste(sFilename, ".txt", sep = ""),
            #            quote = FALSE,
            #            sep = "\t",
            #            col.names = TRUE,
            #            row.names = FALSE,
            #            na = "")
            rm(dGt23k); gc()  # clear the memory
            sFileName[s1Assay] <- paste(sFilename, ".txt", sep = "")
        }
    }
    unlink(sDir, recursive = TRUE)
    options(warn = 0)
    return(sFileName)
}

#' DownloadmiRisoformData: get miRisoform data, assayPlatform %in% c("mirIsoform_GA.hg18", "mirIsoform_GA.hg19", "mirIsoform_HiSeq.hg18", "mirIsoform_HiSeq.hg19")
#'
#' @param cancerType (i.e. sCancer, vector of cancer type), length(sCancer)> = 1
#' @param assayPlatform (i.e. sAssay, vector of type), length(sAssay)> = 1, assayPlatform %in% c("mir_GA.hg18", "mir_GA.hg19", "mir_GA.hg19.mirbase20", "mir_HiSeq.hg18", "mir_HiSeq.hg19", "mir_HiSeq.hg19.mirbase20")
#' @param tissueType (i.e. sampleTypeName, vector of sample_type_name, could be transfered to sample_type_id): c("TP", "TR", ...) -> c(01, 02, ...)
#' @param saveFilePath (string of path to save the merged data): absolute or relative path
#' @param saveFileName (string of filename prefix)
#' @param inputPatientIDs (i.e. sBarcode, vector of barcode), find specified patients
#' @return sFileName (vector of path/filename), could be used by module B
DownloadmiRisoformData <- function(cancerType = NULL,
                                   assayPlatform = NULL,
                                   tissueType = NULL,
                                   saveFilePath = ".",
                                   saveFileName = "",
                                   inputPatientIDs = NULL) {
    sArchive <- "legacy"; sFieldMeta <- ""; nEntityCount <- (-1); sEndpoint <- "files"
    options(warn = -1)
    if (saveFilePath != ".") {dir.create(saveFilePath, recursive = TRUE)}
    sDir <- paste("tmp", TimeNow(), sep = "_"); dir.create(sDir)
    #  if (saveFilePath != ".") {
    #    if (!dir.exists(saveFilePath)) {
    #      dir.create(saveFilePath, recursive = TRUE)
    #    }
    #  } # dir.exists # R >= 3.2
    #  sDir <- paste("tmp", TimeNow(), sep = "_")
    #  if (!dir.exists(sDir)) {dir.create(sDir)} # dir.exists # R >= 3.2
    l <- CheckParam(sCancer = cancerType,
                    sAssay = assayPlatform,
                    sampleTypeName = tissueType,
                    sAssayGroup = "mirIsoform")
    sFileName <- rep(NA, length(l$sAssay))
    names(sFileName) <- l$sAssay
    for (s1Assay in l$sAssay) {
        dPiped <- Pipe(sCancer = l$sCancer,
                       s1Assay = s1Assay,
                       sSampleTypeId = l$sSampleTypeId,
                       sBarcode = inputPatientIDs,
                       sDir = sDir,
                       sArchive = sArchive,
                       sFieldMeta = sFieldMeta,
                       nEntityCount = nEntityCount,
                       sEndpoint = sEndpoint)
        if (!is.null(dPiped)) {
            sFilename <- paste(saveFilePath,
                               ifelse(endsWith(saveFilePath, "/"), "", "/"),
                               ifelse(saveFileName == "", "",
                                      paste(saveFileName, "__", sep = "")),
                               paste(l$sCancer, collapse = "_"),
                               "__",
                               s1Assay,
                               "__",
                               ifelse(is.null(tissueType), "tissueTypeAll",
                                      paste(tissueType, collapse = "_")),
                               "__",
                               TimeNow(),
                               sep = "")
            save(dPiped, file = paste(sFilename, ".rda", sep = ""))
            #write.table(dPiped,
            #            file = paste(sFilename, ".txt", sep = ""),
            #            quote = FALSE,
            #            sep = "\t",
            #            col.names = TRUE,
            #            row.names = FALSE,
            #            na = "")
            rm(dPiped); gc()  # clear the memory
            sFileName[s1Assay] <- paste(sFilename, ".txt", sep = "")
        }
    }
    unlink(sDir, recursive = TRUE)
    options(warn = 0)
    return(sFileName)
}

#' DownloadmiRNASeqDataIncludeIsoform: get miRNASeq data, assayPlatform %in% c("mir_GA.hg18", "mir_GA.hg19", "mir_GA.hg19.mirbase20", "mir_HiSeq.hg18", "mir_HiSeq.hg19", "mir_HiSeq.hg19.mirbase20", "mirIsoform_GA.hg18", "mirIsoform_GA.hg19", "mirIsoform_GA.hg19.mirbase20", "mirIsoform_HiSeq.hg18", "mirIsoform_HiSeq.hg19", "mirIsoform_HiSeq.hg19.mirbase20")
#'
#' @param cancerType (i.e. sCancer, vector of cancer type), length(sCancer)> = 1
#' @param assayPlatform (i.e. sAssay, vector of type), length(sAssay)> = 1, assayPlatform %in% c("mir_GA.hg18", "mir_GA.hg19", "mir_GA.hg19.mirbase20", "mir_HiSeq.hg18", "mir_HiSeq.hg19", "mir_HiSeq.hg19.mirbase20", "mirIsoform_GA.hg18", "mirIsoform_GA.hg19", "mirIsoform_GA.hg19.mirbase20", "mirIsoform_HiSeq.hg18", "mirIsoform_HiSeq.hg19", "mirIsoform_HiSeq.hg19.mirbase20")
#' @param tissueType (i.e. sampleTypeName, vector of sample_type_name, could be transfered to sample_type_id): c("TP", "TR", ...) -> c(01, 02, ...)
#' @param saveFilePath (string of path to save the merged data): absolute or relative path
#' @param saveFileName (string of filename prefix)
#' @param inputPatientIDs (i.e. sBarcode, vector of barcode), find specified patients
#' @return sFileName (vector of path/filename), could be used by module B
DownloadmiRNASeqDataIncludeIsoform <- function(cancerType = NULL,
                                               assayPlatform = NULL,
                                               tissueType = NULL,
                                               saveFilePath = ".",
                                               saveFileName = "",
                                               inputPatientIDs = NULL) {
    sMir <- c("mir_GA.hg18",
              "mir_GA.hg19",
              "mir_HiSeq.hg18",
              "mir_HiSeq.hg19")
    sMirIso <- c("mirIsoform_GA.hg18",
                 "mirIsoform_GA.hg19",
                 "mirIsoform_HiSeq.hg18",
                 "mirIsoform_HiSeq.hg19")
    assayPlatformMir <- assayPlatform[assayPlatform %in% sMir]
    assayPlatformIso <- assayPlatform[assayPlatform %in% sMirIso]
    sFileNameMir <- DownloadmiRNASeqData(cancerType,
                                         assayPlatformMir,
                                         tissueType,
                                         saveFilePath,
                                         saveFileName,
                                         inputPatientIDs)
    sFileNameMirIso <- DownloadmiRisoformData(cancerType,
                                              assayPlatformIso,
                                              tissueType,
                                              saveFilePath,
                                              saveFileName,
                                              inputPatientIDs)
    sFileName <- c(sFileNameMir, sFileNameMirIso)
    return(sFileName)
}

#' DownloadMethylationData: get methylation data, assayPlatform %in% c("methylation_27", "methylation_450")
#'
#' @param cancerType (i.e. sCancer, vector of cancer type), length(sCancer)> = 1
#' @param assayPlatform (i.e. sAssay, vector of type), length(sAssay)> = 1, assayPlatform %in% c("methylation_27", "methylation_450")
#' @param tissueType (i.e. sampleTypeName, vector of sample_type_name, could be transfered to sample_type_id): c("TP", "TR", ...) -> c(01, 02, ...)
#' @param saveFilePath (string of path to save the merged data): absolute or relative path
#' @param saveFileName (string of filename prefix)
#' @param inputPatientIDs (i.e. sBarcode, vector of barcode), find specified patients
#' @return sFileName (vector of path/filename), could be used by module B
DownloadMethylationData <- function(cancerType = NULL,
                                    assayPlatform = NULL,
                                    tissueType = NULL,
                                    saveFilePath = ".",
                                    saveFileName = "",
                                    inputPatientIDs = NULL) {
    sArchive <- "legacy"; sFieldMeta <- ""; nEntityCount <- (-1); sEndpoint <- "files"
    options(warn = -1)
    if (saveFilePath != ".") {dir.create(saveFilePath, recursive = TRUE)}
    sDir <- paste("tmp", TimeNow(), sep = "_"); dir.create(sDir)
    #  if (saveFilePath != ".") {
    #    if (!dir.exists(saveFilePath)) {
    #      dir.create(saveFilePath, recursive = TRUE)
    #    }
    #  } # dir.exists # R >= 3.2
    #  sDir <- paste("tmp", TimeNow(), sep = "_")
    #  if (!dir.exists(sDir)) {dir.create(sDir)} # dir.exists # R >= 3.2
    l <- CheckParam(sCancer = cancerType,
                    sAssay = assayPlatform,
                    sampleTypeName = tissueType,
                    sAssayGroup = "methy")
    sFileName <- rep(NA, length(l$sAssay))
    names(sFileName) <- l$sAssay
    for (s1Assay in l$sAssay) {
        dPiped <- PipeBatch(sCancer = l$sCancer,
                            s1Assay = s1Assay,
                            sSampleTypeId = l$sSampleTypeId,
                            sBarcode = inputPatientIDs,
                            sDir = sDir,
                            sArchive = sArchive,
                            sFieldMeta = sFieldMeta,
                            nEntityCount = nEntityCount,
                            sEndpoint = sEndpoint)
        if (!is.null(dPiped)) {
            sFilename <- paste(saveFilePath,
                               ifelse(endsWith(saveFilePath, "/"), "", "/"),
                               ifelse(saveFileName == "", "",
                                      paste(saveFileName, "__", sep = "")),
                               paste(l$sCancer, collapse = "_"),
                               "__",
                               s1Assay,
                               "__",
                               ifelse(is.null(tissueType), "tissueTypeAll",
                                      paste(tissueType, collapse = "_")),
                               "__",
                               TimeNow(),
                               sep = "")
            save(dPiped, file = paste(sFilename, ".rda", sep = ""))
            #write.table(dPiped,
            #            file = paste(sFilename, ".txt", sep = ""),
            #            quote = FALSE,
            #            sep = "\t",
            #            col.names = TRUE,
            #            row.names = FALSE,
            #            na = "")
            rm(dPiped); gc()  # clear the memory
            sFileName[s1Assay] <- paste(sFilename, ".txt", sep = "")
        }
    }
    unlink(sDir, recursive = TRUE)
    options(warn = 0)
    return(sFileName)
}

#' DownloadCNAData: get copy number data, assayPlatform %in% c("cna_cnv.hg18", "cna_cnv.hg19", "cna_nocnv.hg18", "cna_nocnv.hg19")
#'
#' @param cancerType (i.e. sCancer, vector of cancer type), length(sCancer)> = 1
#' @param assayPlatform (i.e. sAssay, vector of type), length(sAssay)> = 1, assayPlatform %in% c("cna_cnv.hg18", "cna_cnv.hg19", "cna_nocnv.hg18", "cna_nocnv.hg19")
#' @param tissueType (i.e. sampleTypeName, vector of sample_type_name, could be transfered to sample_type_id): c("TP", "TR", ...) -> c(01, 02, ...)
#' @param saveFilePath (string of path to save the merged data): absolute or relative path
#' @param saveFileName (string of filename prefix)
#' @param inputPatientIDs (i.e. sBarcode, vector of barcode), find specified patients
#' @return sFileName (vector of path/filename), could be used by module B
DownloadCNAData <- function(cancerType = NULL,
                            assayPlatform = NULL,
                            tissueType = NULL,
                            saveFilePath = ".",
                            saveFileName = "",
                            inputPatientIDs = NULL) {
    sArchive <- "legacy"; sFieldMeta <- ""; nEntityCount <- (-1); sEndpoint <- "files"
    options(warn = -1)
    if (saveFilePath != ".") {dir.create(saveFilePath, recursive = TRUE)}
    sDir <- paste("tmp", TimeNow(), sep = "_"); dir.create(sDir)
    #  if (saveFilePath != ".") {
    #    if (!dir.exists(saveFilePath)) {
    #      dir.create(saveFilePath, recursive = TRUE)
    #    }
    #  } # dir.exists # R >= 3.2
    #  sDir <- paste("tmp", TimeNow(), sep = "_")
    #  if (!dir.exists(sDir)) {dir.create(sDir)} # dir.exists # R >= 3.2
    l <- CheckParam(sCancer = cancerType,
                    sAssay = assayPlatform,
                    sampleTypeName = tissueType,
                    sAssayGroup = "cna")
    sFileName <- rep(NA, length(l$sAssay))
    names(sFileName) <- l$sAssay
    for (s1Assay in l$sAssay) {
        dPiped <- Pipe(sCancer = l$sCancer,
                       s1Assay = s1Assay,
                       sSampleTypeId = l$sSampleTypeId,
                       sBarcode = inputPatientIDs,
                       sDir = sDir,
                       sArchive = sArchive,
                       sFieldMeta = sFieldMeta,
                       nEntityCount = nEntityCount,
                       sEndpoint = sEndpoint)
        if (!is.null(dPiped)) {
            colnames(dPiped)[1] <- c("Sample")  # for Module_B
            sFilename <- paste(saveFilePath,
                               ifelse(endsWith(saveFilePath, "/"), "", "/"),
                               ifelse(saveFileName == "", "",
                                      paste(saveFileName, "__", sep = "")),
                               paste(l$sCancer, collapse = "_"),
                               "__",
                               s1Assay,
                               "__",
                               ifelse(is.null(tissueType), "tissueTypeAll",
                                      paste(tissueType, collapse = "_")),
                               "__",
                               TimeNow(),
                               sep = "")
            save(dPiped, file = paste(sFilename, ".rda", sep = ""))
            #write.table(dPiped,
            #            file = paste(sFilename, ".txt", sep = ""),
            #            quote = FALSE,
            #            sep = "\t",
            #            col.names = TRUE,
            #            row.names = FALSE,
            #            na = "")
            rm(dPiped); gc()  # clear the memory
            sFileName[s1Assay] <- paste(sFilename, ".txt", sep = "")
        }
    }
    unlink(sDir, recursive = TRUE)
    options(warn = 0)
    return(sFileName)
}

#' DownloadRNASeqData: get gene expression data, assayPlatform %in% c("gene_Array", "gene.normalized_RNAseq", "gene_RNAseq", "isoform.normalized_RNAseq", "isoform_RNAseq", "exon_RNAseq", "exonJunction_RNAseq")
#'
#' @param cancerType (i.e. sCancer, vector of cancer type), length(sCancer)> = 1
#' @param assayPlatform (i.e. sAssay, vector of type), length(sAssay)> = 1, assayPlatform %in% c("gene_Array", "gene.normalized_RNAseq", "gene_RNAseq", "isoform.normalized_RNAseq", "isoform_RNAseq", "exon_RNAseq", "exonJunction_RNAseq")
#' @param tissueType (i.e. sampleTypeName, vector of sample_type_name, could be transfered to sample_type_id): c("TP", "TR", ...) -> c(01, 02, ...)
#' @param saveFilePath (string of path to save the merged data): absolute or relative path
#' @param saveFileName (string of filename prefix)
#' @param inputPatientIDs (i.e. sBarcode, vector of barcode), find specified patients
#' @return sFileName (vector of path/filename), could be used by module B
DownloadRNASeqData <- function(cancerType = NULL,
                               assayPlatform = NULL,
                               tissueType = NULL,
                               saveFilePath = ".",
                               saveFileName = "",
                               inputPatientIDs = NULL) {
    sArchive <- "legacy"; sFieldMeta <- ""; nEntityCount <- (-1); sEndpoint <- "files"
    options(warn = -1)
    if (saveFilePath != ".") {dir.create(saveFilePath, recursive = TRUE)}
    sDir <- paste("tmp", TimeNow(), sep = "_"); dir.create(sDir)
    #  if (saveFilePath != ".") {
    #    if (!dir.exists(saveFilePath)) {
    #      dir.create(saveFilePath, recursive = TRUE)
    #    }
    #  } # dir.exists # R >= 3.2
    #  sDir <- paste("tmp", TimeNow(), sep = "_")
    #  if (!dir.exists(sDir)) {dir.create(sDir)} # dir.exists # R >= 3.2
    l <- CheckParam(sCancer = cancerType,
                    sAssay = assayPlatform,
                    sampleTypeName = tissueType,
                    sAssayGroup = "gene")
    sFileName <- rep(NA, length(l$sAssay))
    names(sFileName) <- l$sAssay
    for (s1Assay in l$sAssay) {
        if (s1Assay %in% c("isoform.normalized_RNAseq",
                           "isoform_RNAseq",
                           "exon_RNAseq",
                           "exonJunction_RNAseq")) {
            dPiped <- PipeBatch(sCancer = l$sCancer,
                                s1Assay = s1Assay,
                                sSampleTypeId = l$sSampleTypeId,
                                sBarcode = inputPatientIDs,
                                sDir = sDir,
                                sArchive = sArchive,
                                sFieldMeta = sFieldMeta,
                                nEntityCount = nEntityCount,
                                sEndpoint = sEndpoint)
        } else {
            dPiped <- Pipe(sCancer = l$sCancer,
                           s1Assay = s1Assay,
                           sSampleTypeId = l$sSampleTypeId,
                           sBarcode = inputPatientIDs,
                           sDir = sDir,
                           sArchive = sArchive,
                           sFieldMeta = sFieldMeta,
                           nEntityCount = nEntityCount,
                           sEndpoint = sEndpoint)
        }
        if (!is.null(dPiped)) {
            sFilename <- paste(saveFilePath,
                               ifelse(endsWith(saveFilePath, "/"), "", "/"),
                               ifelse(saveFileName == "", "",
                                      paste(saveFileName, "__", sep = "")),
                               paste(l$sCancer, collapse = "_"),
                               "__",
                               s1Assay,
                               "__",
                               ifelse(is.null(tissueType), "tissueTypeAll",
                                      paste(tissueType, collapse = "_")),
                               "__",
                               TimeNow(),
                               sep = "")
            save(dPiped, file = paste(sFilename, ".rda", sep = ""))
            #write.table(dPiped,
            #            file = paste(sFilename, ".txt", sep = ""),
            #            quote = FALSE,
            #            sep = "\t",
            #            col.names = TRUE,
            #            row.names = FALSE,
            #            na = "")
            rm(dPiped); gc()  # clear the memory
            sFileName[s1Assay] <- paste(sFilename, ".txt", sep = "")
        }
    }
    unlink(sDir, recursive = TRUE)
    options(warn = 0)
    return(sFileName)
}

#' DownloadRPPAData: get protein expression data, assayPlatform %in% c("protein_RPPA")
#'
#' @param cancerType (i.e. sCancer, vector of cancer type), length(sCancer)> = 1
#' @param assayPlatform (i.e. sAssay, vector of type), length(sAssay)> = 1, assayPlatform %in% c("protein_RPPA")
#' @param tissueType (i.e. sampleTypeName, vector of sample_type_name, could be transfered to sample_type_id): c("TP", "TR", ...) -> c(01, 02, ...)
#' @param saveFilePath (string of path to save the merged data): absolute or relative path
#' @param saveFileName (string of filename prefix)
#' @param inputPatientIDs (i.e. sBarcode, vector of barcode), find specified patients
#' @return sFileName (vector of path/filename), could be used by module B
DownloadRPPAData <- function(cancerType = NULL,
                             assayPlatform = NULL,
                             tissueType = NULL,
                             saveFilePath = ".",
                             saveFileName = "",
                             inputPatientIDs = NULL) {
    sArchive <- "legacy"; sFieldMeta <- ""; nEntityCount <- (-1); sEndpoint <- "files"
    options(warn = -1)
    if (saveFilePath != ".") {dir.create(saveFilePath, recursive = TRUE)}
    sDir <- paste("tmp", TimeNow(), sep = "_"); dir.create(sDir)
    #  if (saveFilePath != ".") {
    #    if (!dir.exists(saveFilePath)) {
    #      dir.create(saveFilePath, recursive = TRUE)
    #    }
    #  } # dir.exists # R >= 3.2
    #  sDir <- paste("tmp", TimeNow(), sep = "_")
    #  if (!dir.exists(sDir)) {dir.create(sDir)} # dir.exists # R >= 3.2
    l <- CheckParam(sCancer = cancerType,
                    sAssay = assayPlatform,
                    sampleTypeName = tissueType,
                    sAssayGroup = "protein")
    sFileName <- rep(NA, length(l$sAssay))
    names(sFileName) <- l$sAssay
    for (s1Assay in l$sAssay) {
        dPiped <- Pipe(sCancer = l$sCancer,
                       s1Assay = s1Assay,
                       sSampleTypeId = l$sSampleTypeId,
                       sBarcode = inputPatientIDs,
                       sDir = sDir,
                       sArchive = sArchive,
                       sFieldMeta = sFieldMeta,
                       nEntityCount = nEntityCount,
                       sEndpoint = sEndpoint)
        if (!is.null(dPiped)) {
            sFilename <- paste(saveFilePath,
                               ifelse(endsWith(saveFilePath, "/"), "", "/"),
                               ifelse(saveFileName == "", "",
                                      paste(saveFileName, "__", sep = "")),
                               paste(l$sCancer, collapse = "_"),
                               "__",
                               s1Assay,
                               "__",
                               ifelse(is.null(tissueType), "tissueTypeAll",
                                      paste(tissueType, collapse = "_")),
                               "__",
                               TimeNow(),
                               sep = "")
            save(dPiped, file = paste(sFilename, ".rda", sep = ""))
            #write.table(dPiped,
            #            file = paste(sFilename, ".txt", sep = ""),
            #            quote = FALSE,
            #            sep = "\t",
            #            col.names = TRUE,
            #            row.names = FALSE,
            #            na = "")
            rm(dPiped); gc()  # clear the memory
            sFileName[s1Assay] <- paste(sFilename, ".txt", sep = "")
        }
    }
    unlink(sDir, recursive = TRUE)
    options(warn = 0)
    return(sFileName)
}

#' DownloadSomaticMutationData: get somatic mutation, assayPlatform %in% c("somaticMutation_DNAseq")
#'
#' @param cancerType (i.e. sCancer, vector of cancer type), length(sCancer)> = 1
#' @param assayPlatform (i.e. sAssay, vector of type), length(sAssay)> = 1, assayPlatform %in% c("somaticMutation_DNAseq")
#' @param tissueType (i.e. sampleTypeName, vector of sample_type_name, could be transfered to sample_type_id): c("TP", "TR", ...) -> c(01, 02, ...)
#' @param saveFilePath (string of path to save the merged data): absolute or relative path
#' @param saveFileName (string of filename prefix)
#' @param inputPatientIDs (i.e. sBarcode, vector of barcode), find specified patients
#' @return sFileName (vector of path/filename), could be used by module B
#' @examples
#' v <- DownloadSomaticMutationData(cancerType = "BRCA", assayPlatform = NULL, tissueType = NULL, saveFilePath = ".", saveFileName = "", inputPatientIDs = NULL)
DownloadSomaticMutationData <- function(cancerType = NULL,
                                        assayPlatform = NULL,
                                        tissueType = NULL,
                                        saveFilePath = ".",
                                        saveFileName = "",
                                        inputPatientIDs = NULL) {
    sArchive <- "legacy"; sFieldMeta <- ""; nEntityCount <- (-1); sEndpoint <- "files"
    options(warn = -1)
    if (saveFilePath != ".") {dir.create(saveFilePath, recursive = TRUE)}
    sDir <- paste("tmp", TimeNow(), sep = "_"); dir.create(sDir)
    #  if (saveFilePath != ".") {
    #    if (!dir.exists(saveFilePath)) {
    #      dir.create(saveFilePath, recursive = TRUE)
    #    }
    #  } # dir.exists # R >= 3.2
    #  sDir <- paste("tmp", TimeNow(), sep = "_")
    #  if (!dir.exists(sDir)) {dir.create(sDir)} # dir.exists # R >= 3.2
    l <- CheckParam(sCancer = cancerType,
                    sAssay = assayPlatform,
                    sampleTypeName = tissueType,
                    sAssayGroup = "somatic")
    sFileName <- NULL
    for (s1Cancer in l$sCancer) {
        ldPiped <- PipeSomatic(sCancer = s1Cancer,
                               s1Assay = l$sAssay,
                               sSampleTypeId = l$sSampleTypeId,
                               sBarcode = inputPatientIDs,
                               sDir = sDir,
                               sArchive = sArchive,
                               sFieldMeta = sFieldMeta,
                               nEntityCount = nEntityCount,
                               sEndpoint = sEndpoint)
        for (s1Pathname in names(ldPiped)) {
            dPiped <- ldPiped[[s1Pathname]]
            if (!is.null(dPiped)) {
                sFilename <- paste(saveFilePath,
                                   ifelse(endsWith(saveFilePath, "/"), "", "/"),
                                   ifelse(saveFileName == "", "",
                                          paste(saveFileName, "__", sep = "")),
                                   paste(s1Cancer, collapse = "_"),
                                   "__",
                                   l$sAssay,
                                   "__",
                                   ifelse(is.null(tissueType), "tissueTypeAll",
                                          paste(tissueType, collapse = "_")),
                                   "__",
                                   TimeNow(),
                                   "__",
                                   rev(strsplit(s1Pathname, split = "/")[[1]])[1],
                                   sep = "")
                save(dPiped, file = paste(sFilename, ".rda", sep = ""))
                #write.table(dPiped,
                #            file = paste(sFilename, ".txt", sep = ""),
                #            quote = FALSE,
                #            sep = "\t",
                #            col.names = TRUE,
                #            row.names = FALSE,
                #            na = "")
                sFileName <- c(sFileName, paste(sFilename, ".txt", sep = ""))
            }
            if (.Platform$OS.type == "windows") {
                dir.create(paste(saveFilePath,
                                 ifelse(endsWith(saveFilePath, "/"), "", "/"),
                                 "originalSomaticMutationFiles",
                                 sep = ""),
                           recursive = TRUE)
                bFileRename <-file.rename(from = s1Pathname,
                                          to = paste(saveFilePath,
                                                     ifelse(endsWith(saveFilePath, "/"), "", "/"),
                                                     "originalSomaticMutationFiles/",
                                                     rev(strsplit(s1Pathname,
                                                                  split = "/")[[1]])[1],
                                                     sep = ""))
                stopifnot(all(bFileRename))
            }
        }
    }
    rm(ldPiped); gc()  # clear the memory
    unlink(sDir, recursive = TRUE)
    options(warn = 0)
    return(sFileName)
}

#' DownloadCPTACData: get CPTAC data, assayPlatform %in% c("glycoproteome_iTRAQ", "phosphoproteome_iTRAQ", "proteome_iTRAQ")
#'
#' @param cancerType (i.e. sCancer, vector of cancer type), length(sCancer)> = 1. Now only c("BRCA", "OV", "COAD", "READ"), "BRCA"->"Breast", "OV"->"OV", c("COAD", "READ")->"Colorectal"
#' @param assayPlatform (i.e. sAssay, vector of type), length(sAssay)> = 1, assayPlatform %in% c("glycoproteome_iTRAQ", "phosphoproteome_iTRAQ", "proteome_iTRAQ")
#' @param tissueType (i.e. sampleTypeName, vector of sample_type_name, could be transfered to sample_type_id): c("TP", "TR", ...) -> c(01, 02, ...). Now only "TP"(01)
#' @param saveFilePath (string of path to save the merged data): absolute or relative path
#' @param saveFileName (string of filename prefix)
#' @param inputPatientIDs (i.e. sBarcode, vector of barcode), find specified patients
#' @return sFileName (vector of path/filename), could be used by module B
#' @examples
#' v <- DownloadCPTACData(cancerType = NULL, assayPlatform = NULL, tissueType = NULL, saveFilePath = ".", saveFileName = "", inputPatientIDs = NULL)
DownloadCPTACData <- function(cancerType = NULL,
                              assayPlatform = NULL,
                              tissueType = NULL,
                              saveFilePath = ".",
                              saveFileName = "",
                              inputPatientIDs = NULL) {
    sBarcode <- inputPatientIDs
    options(warn = -1)
    if (saveFilePath != ".") {dir.create(saveFilePath, recursive = TRUE)}
    sDir <- paste("tmp", TimeNow(), sep = "_"); dir.create(sDir)
    #  if (saveFilePath != ".") {
    #    if (!dir.exists(saveFilePath)) {
    #      dir.create(saveFilePath, recursive = TRUE)
    #    }
    #  } # dir.exists # R >= 3.2
    #  sDir <- paste("tmp", TimeNow(), sep = "_")
    #  if (!dir.exists(sDir)) {dir.create(sDir)} # dir.exists # R >= 3.2
    l <- CheckParam(sCancer = cancerType,
                    sAssay = assayPlatform,
                    sampleTypeName = tissueType,
                    sAssayGroup = "itraq")
    sCancerCptac <- c("BRCA", "OV", "COAD", "READ")
    sCancer <- intersect(l$sCancer, sCancerCptac)
    if (is.null(sCancer)) {
        print(c("cancerType should be 'NULL' (for all cancerType) or one of: ",
                sCancerCptac))
        stopifnot(is.null(sCancer))
    }
    sUrlPre <- "https://cptc-xfer.uis.georgetown.edu/publicData/Phase_II_Data/"
    print("CPTAC files  : downloading ...")
    sFileName <- NULL
    for (s1Cancer in sCancer) {
        if (s1Cancer %in% c("BRCA")) {
            sPathname <- c("TCGA_Breast_Cancer/TCGA_Breast_BI_Proteome_CDAP_Protein_Report.r3/TCGA_Breast_BI_Proteome.itraq.tsv",
                           ##"TCGA_Breast_Cancer/TCGA_Breast_BI_Proteome_CDAP_Protein_Report.r2/TCGA_Breast_BI_Proteome_CDAP.r2.itraq.tsv",
                           # "TCGA_Breast_Cancer/TCGA_Breast_BI_Proteome_CDAP_Protein_Report.r2/TCGA_Breast_BI_Proteome_CDAP.r2.peptides.tsv",
                           # "TCGA_Breast_Cancer/TCGA_Breast_BI_Phosphoproteome_CDAP_Protein_Report.r3/TCGA_Breast_BI_Phosphoproteome.peptides.tsv",
                           # "TCGA_Breast_Cancer/TCGA_Breast_BI_Phosphoproteome_CDAP_Protein_Report.r3/TCGA_Breast_BI_Phosphoproteome.phosphopeptide.itraq.tsv",
                           "TCGA_Breast_Cancer/TCGA_Breast_BI_Phosphoproteome_CDAP_Protein_Report.r4/TCGA_Breast_BI_Phosphoproteome.phosphosite.itraq.tsv")
            ##"TCGA_Breast_Cancer/TCGA_Breast_BI_Phosphoproteome_CDAP_Protein_Report.r3/TCGA_Breast_BI_Phosphoproteome.phosphosite.itraq.tsv")
        } else if (s1Cancer %in% c("COAD", "READ")) {
            sPathname <- c(# "TCGA_Colorectal_Cancer/TCGA_Colon_VU_Proteome_CDAP_Protein_Report.r2/TCGA_Colon_VU_Proteome_CDAP.r2.peptides.tsv",
                           # "TCGA_Colorectal_Cancer/TCGA_Colon_VU_Proteome_CDAP_Protein_Report.r2/TCGA_Colon_VU_Proteome_CDAP.r2.precursor_area.tsv",
                           "TCGA_Colorectal_Cancer/TCGA_Colon_VU_Proteome_CDAP_Protein_Report.r2/TCGA_Colon_VU_Proteome_CDAP.r2.spectral_counts.tsv")
        } else if (s1Cancer %in% c("OV")) {
            sPathname <- c("TCGA_Ovarian_Cancer/TCGA_Ovarian_JHU_Proteome_CDAP_Protein_Report.r3/TCGA_Ovarian_JHU_Proteome.itraq.tsv",
                           ##"TCGA_Ovarian_Cancer/TCGA_Ovarian_JHU_Proteome_CDAP_Protein_Report.r2/TCGA_Ovarian_JHU_Proteome_CDAP.r2.itraq.tsv",
                           # "TCGA_Ovarian_Cancer/TCGA_Ovarian_JHU_Proteome_CDAP_Protein_Report.r2/TCGA_Ovarian_JHU_Proteome_CDAP.r2.peptides.tsv",
                           # "TCGA_Ovarian_Cancer/TCGA_Ovarian_JHU_Glycoproteome_CDAP_Protein_Report.r3/TCGA_Ovarian_JHU_Glycoproteome.glycopeptide.itraq.tsv",
                           "TCGA_Ovarian_Cancer/TCGA_Ovarian_JHU_Glycoproteome_CDAP_Protein_Report.r4/TCGA_Ovarian_JHU_Glycoproteome.glycosite.itraq.tsv",
                           ##"TCGA_Ovarian_Cancer/TCGA_Ovarian_JHU_Glycoproteome_CDAP_Protein_Report.r3/TCGA_Ovarian_JHU_Glycoproteome.glycosite.itraq.tsv",
                           # "TCGA_Ovarian_Cancer/TCGA_Ovarian_JHU_Glycoproteome_CDAP_Protein_Report.r3/TCGA_Ovarian_JHU_Glycoproteome.peptides.tsv",
                           "TCGA_Ovarian_Cancer/TCGA_Ovarian_PNNL_Proteome_CDAP_Protein_Report.r3/TCGA_Ovarian_PNNL_Proteome.itraq.tsv",
                           ##"TCGA_Ovarian_Cancer/TCGA_Ovarian_PNNL_Proteome_CDAP_Protein_Report.r2/TCGA_Ovarian_PNNL_Proteome_CDAP.r2.itraq.tsv",
                           # "TCGA_Ovarian_Cancer/TCGA_Ovarian_PNNL_Proteome_CDAP_Protein_Report.r2/TCGA_Ovarian_PNNL_Proteome_CDAP.r2.peptides.tsv",
                           # "TCGA_Ovarian_Cancer/TCGA_Ovarian_PNNL_Phosphoproteome_CDAP_Protein_Report.r3/TCGA_Ovarian_PNNL_Phosphoproteome.peptides.tsv",
                           # "TCGA_Ovarian_Cancer/TCGA_Ovarian_PNNL_Phosphoproteome_CDAP_Protein_Report.r3/TCGA_Ovarian_PNNL_Phosphoproteome.phosphopeptide.itraq.tsv",
                           "TCGA_Ovarian_Cancer/TCGA_Ovarian_PNNL_Phosphoproteome_CDAP_Protein_Report.r4/TCGA_Ovarian_PNNL_Phosphoproteome.phosphosite.itraq.tsv")
            ##"TCGA_Ovarian_Cancer/TCGA_Ovarian_PNNL_Phosphoproteome_CDAP_Protein_Report.r3/TCGA_Ovarian_PNNL_Phosphoproteome.phosphosite.itraq.tsv")
        }
        for (s1Assay in l$sAssay) {
            lFilter <- Filter(s1Assay)
            for (s1Pathname in sPathname[grep(lFilter$file_name,sPathname)]) {
                sFilename <- rev(strsplit(s1Pathname, split = "/")[[1]])[1]
                sUrl <- paste(sUrlPre, s1Pathname, sep = "")
                sOutfile <- paste(sDir, "/", sFilename, sep = "")
                sOption <- paste("--silent --show-error -o", sOutfile)
                sArgument <- paste(sOption, sUrl)
                stdOut <- system2("curl", sArgument, stdout = TRUE)
                if (!is.null(attr(stdOut, "status"))) {
                    print("error (download): check the proxy")
                }
                stopifnot(is.null(attr(stdOut, "status")))
                # sPathName <- paste(saveFilePath, ifelse(endsWith(saveFilePath, "/"), "", "/"), sFilename, sep = "")
                # bFileRename <- file.rename(from = sOutfile, to = sPathName)
                # stopifnot(all(bFileRename))
                # sFilename <- rev(strsplit(sPathName, split = "/")[[1]])[1]
                d <- read.csv(sOutfile,
                              sep = "\t",
                              skip = 0,
                              row.names = NULL,
                              as.is = TRUE,
                              na.strings = "",
                              check.names = FALSE)
                if (s1Assay == "proteome_iTRAQ") {
                    if (s1Cancer %in% c("COAD", "READ")) {
                        sRowNot <- c("Total")
                        sColNot <- c("Total Spectral Counts",
                                     "Total Unshared Spectral Counts")
                        sColDes <- c("Gene",
                                     "Description",
                                     "Organism",
                                     "Chromosome",
                                     "Locus")
                    } else {
                        sRowNot <- c("Mean", "Median", "StdDev")
                        sColNot <- NULL
                        sColDes <- c("Gene",
                                     "Description",
                                     "Organism",
                                     "Chromosome",
                                     "Locus")
                    }
                } else if (s1Assay == "phosphoproteome_iTRAQ") {
                    sRowNot <- NULL
                    sColNot <- NULL
                    sColDes <- c("Phosphosite","Peptide", "Gene", "Organism")
                } else if (s1Assay == "glycoproteome_iTRAQ") {
                    sRowNot <- NULL
                    sColNot <- NULL
                    sColDes <- c("Glycosite","Peptide", "Gene", "Organism")
                }
                d <- d[!(d[,sColDes[1]] %in% sRowNot), !(colnames(d) %in% sColNot)]
                mInfo <- as.matrix(d[, sColDes])
                mData <- as.matrix(d[, setdiff(colnames(d), sColDes)])
                for (n in grep("^X[0-9A-Za-z][0-9A-Za-z]\\.", colnames(mData))) {
                    colnames(mData)[n] <- sub("^X", "", colnames(mData)[n])
                }
                for (n in grep("^[0-9A-Za-z][0-9A-Za-z]-", colnames(mData))) {
                    colnames(mData)[n] <- paste("TCGA-",
                                                #gsub("\\.", "-", colnames(mData)[n]),
                                                colnames(mData)[n],
                                                sep = "")
                }
                for (n in grep("-[0-9][0-9]* (Log|Unshared|Spectral)", colnames(mData))) {
                    substr(colnames(mData)[n], 17, 17) <- ":"
                }
                for (n in grep("^OVARIAN-CONTROL\\.", colnames(mData))) {
                    colnames(mData)[n] <- gsub("\\.", ":", colnames(mData)[n])
                    #colnames(mData)[n] <- colnames(mData)[n]
                    substr(colnames(mData)[n], 16, 16) <- ":"
                }
                if (!is.null(sBarcode)) {
                    bCol <- ifelse(substr(colnames(mData), 1, 12) %in%
                                   substr(sBarcode, 1, 12), TRUE, FALSE)
                    if (!any(bCol)) {
                        print("CPTAC files  : no inputPatientIDs found!")
                    }
                    mData <- mData[, bCol, drop = FALSE]
                }
                if (length(l$sSampleTypeId) < 14) {
                    bSampleTypeId <- ifelse(substr(colnames(mData), 14, 15) %in%
                                            l$sSampleTypeId, TRUE, FALSE)
                    if (!any(bSampleTypeId)) {
                        print("CPTAC files  : no tissueType found!")
                    }
                    mData <- mData[, bSampleTypeId, drop = FALSE]
                }
                if (s1Cancer == "COAD") {
                    sPatient <-
                        c("TCGA-A6-3807", "TCGA-A6-3808", "TCGA-A6-3810",
                          "TCGA-AA-3518", "TCGA-AA-3525", "TCGA-AA-3526",
                          "TCGA-AA-3529", "TCGA-AA-3531", "TCGA-AA-3534",
                          "TCGA-AA-3552", "TCGA-AA-3554", "TCGA-AA-3558",
                          "TCGA-AA-3561", "TCGA-AA-3664", "TCGA-AA-3666",
                          "TCGA-AA-3672", "TCGA-AA-3684", "TCGA-AA-3695",
                          "TCGA-AA-3710", "TCGA-AA-3715", "TCGA-AA-3818",
                          "TCGA-AA-3848", "TCGA-AA-3864", "TCGA-AA-3986",
                          "TCGA-AA-3989", "TCGA-AA-A004", "TCGA-AA-A00A",
                          "TCGA-AA-A00E", "TCGA-AA-A00F", "TCGA-AA-A00J",
                          "TCGA-AA-A00K", "TCGA-AA-A00N", "TCGA-AA-A00O",
                          "TCGA-AA-A00R", "TCGA-AA-A00U", "TCGA-AA-A010",
                          "TCGA-AA-A017", "TCGA-AA-A01C", "TCGA-AA-A01D",
                          "TCGA-AA-A01F", "TCGA-AA-A01I", "TCGA-AA-A01K",
                          "TCGA-AA-A01P", "TCGA-AA-A01R", "TCGA-AA-A01S",
                          "TCGA-AA-A01T", "TCGA-AA-A01V", "TCGA-AA-A01X",
                          "TCGA-AA-A01Z", "TCGA-AA-A022", "TCGA-AA-A024",
                          "TCGA-AA-A029", "TCGA-AA-A02E", "TCGA-AA-A02H",
                          "TCGA-AA-A02J", "TCGA-AA-A02O", "TCGA-AA-A02R",
                          "TCGA-AA-A02Y", "TCGA-AA-A03F", "TCGA-AA-A03J")
                    bPatient <- substr(colnames(mData), 1, 12) %in% sPatient
                    mData <- mData[, bPatient, drop = FALSE]
                }
                if (s1Cancer == "READ") {
                    sPatient <-
                        c("TCGA-AF-2691", "TCGA-AF-2692", "TCGA-AF-3400",
                          "TCGA-AF-3913", "TCGA-AG-3574", "TCGA-AG-3580",
                          "TCGA-AG-3584", "TCGA-AG-3593", "TCGA-AG-3594",
                          "TCGA-AG-4007", "TCGA-AG-A002", "TCGA-AG-A008",
                          "TCGA-AG-A00C", "TCGA-AG-A00H", "TCGA-AG-A00Y",
                          "TCGA-AG-A011", "TCGA-AG-A014", "TCGA-AG-A015",
                          "TCGA-AG-A016", "TCGA-AG-A01J", "TCGA-AG-A01L",
                          "TCGA-AG-A01N", "TCGA-AG-A01W", "TCGA-AG-A01Y",
                          "TCGA-AG-A020", "TCGA-AG-A026", "TCGA-AG-A02N",
                          "TCGA-AG-A02X", "TCGA-AG-A032", "TCGA-AG-A036")
                    bPatient <- substr(colnames(mData), 1, 12) %in% sPatient
                    mData <- mData[, bPatient, drop = FALSE]
                }
                sFileName1 <- paste(saveFilePath,
                                    ifelse(endsWith(saveFilePath, "/"), "", "/"),
                                    ifelse(saveFileName == "", "",
                                           paste(saveFileName, "__", sep = "")),
                                    paste(s1Cancer, collapse = "_"),
                                    "__",
                                    s1Assay,
                                    "__",
                                    ifelse(is.null(tissueType), "tissueTypeAll",
                                           paste(tissueType, collapse = "_")),
                                    "__",
                                    # ifelse(s1Cancer == "OV",
                                    strsplit(sFilename, split = "_")[[1]][3],
                                    #      ""),
                                    "__",
                                    TimeNow(),
                                    ".txt",
                                    sep = "")
                dPiped <- cbind(mInfo, mData)
                save(dPiped, file = paste(sFileName1, ".rda", sep = ""))
                #write.table(dPiped,
                #            file = sFileName1,
                #            quote = FALSE,
                #            sep = "\t",
                #            col.names = TRUE,
                #            row.names = FALSE,
                #            na = "")
                rm(mInfo,mData); gc()  # clear the memory
                sFileName <- c(sFileName, sFileName1)
                # #
                # vnUnshared <- grep("Unshared", colnames(mData))
                # mDataShUn <- mData[, vnUnshared - 1]  # barcode-Log-Ratio
                # mDataUnsh <- mData[, vnUnshared    ]  # barcode-Unshared-Log-Ratio
                # colnames(mDataShUn) <- unlist(strsplit(colnames(mDataShUn),
                #                                        split = "-Log-Ratio"))
                # colnames(mDataUnsh) <- unlist(strsplit(colnames(mDataUnsh),
                #                                        split = "-Unshared-Log-Ratio"))
                # substr(colnames(mDataShUn), 17, 17) <- ":"
                # substr(colnames(mDataUnsh), 17, 17) <- ":"
                # mDataShUn <- mDataShUn[, order(colnames(mDataShUn))]
                # mDataUnsh <- mDataUnsh[, order(colnames(mDataUnsh))]
                # stopifnot(all(colnames(mDataUnsh) == colnames(mDataShUn)))
                # if (!is.null(sBarcode)) {
                #   bCol <- ifelse(substr(colnames(mDataShUn), 1, 12) %in%
                #                   substr(sBarcode, 1, 12), TRUE, FALSE)
                #   if (!any(bCol)) {
                #     print("CPTAC files  : no inputPatientIDs found!")
                #   }
                #   mDataShUn <- mDataShUn[, bCol, drop = FALSE]
                #   mDataUnsh <- mDataUnsh[, bCol, drop = FALSE]
                # }
                # if (length(l$sSampleTypeId) < 14) {
                #   bSampleTypeId <- ifelse(substr(colnames(mDataShUn), 14, 15) %in%
                #                            l$sSampleTypeId, TRUE, FALSE)
                #   if (!any(bSampleTypeId)) {
                #     print("CPTAC files  : no tissueType found!")
                #   }
                #   mDataShUn <- mDataShUn[, bSampleTypeId, drop = FALSE]
                #   mDataUnsh <- mDataUnsh[, bSampleTypeId, drop = FALSE]
                # }
                # sFileNameShUn <- paste(saveFilePath,
                #                        ifelse(endsWith(saveFilePath, "/"), "", "/"),
                #                        ifelse(saveFileName == "", "",
                #                               paste(saveFileName, "__", sep = "")),
                #                        paste(s1Cancer, collapse = "_"),
                #                        "__",
                #                        l$sAssay,
                #                        "__",
                #                        strsplit(sFilename, split = "\\.")[[1]][1],
                #                        "__LogRatio",
                #                        "__",
                #                        TimeNow(),
                #                        ".txt",
                #                        sep = "")
                # sFileNameUnsh <- paste(saveFilePath,
                #                        ifelse(endsWith(saveFilePath, "/"), "", "/"),
                #                        ifelse(saveFileName == "", "",
                #                               paste(saveFileName, "__", sep = "")),
                #                        paste(s1Cancer, collapse = "_"),
                #                        "__",
                #                        l$sAssay,
                #                        "__",
                #                        strsplit(sFilename, split = "\\.")[[1]][1],
                #                        "__LogRatio_Unshared",
                #                        "__",
                #                        TimeNow(),
                #                        ".txt",
                #                        sep = "")
                ## write.table(cbind(mInfo, mDataShUn),
                ##             file = sFileNameShUn,
                ##             quote = FALSE,
                ##             sep = "\t",
                ##             col.names = TRUE,
                ##             row.names = FALSE,
                ##             na = "")
                ## write.table(cbind(mInfo, mDataUnsh),
                ##             file = sFileNameUnsh,
                ##             quote = FALSE,
                ##             sep = "\t",
                ##             col.names = TRUE,
                ##             row.names = FALSE,
                ##             na = "")
                # sFileName <- c(sFileName, sFileNameShUn, sFileNameUnsh)
                # #
        }}
    }
    unlink(sDir, recursive = TRUE)
    options(warn = 0)
    print("CPTAC files  : downloading done!")
    return(sFileName)
}

#' DownloadBiospecimenClinicalData: get biospecimen and clinical data
#'
#' @param cancerType String indicating the specified cancer type
#' for which data should be downloaded.
#' Its value can be one of the cancer type abbreviations:
#' \code{c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA",
#' "GBM", HNSC, KICH, KIRC, KIRP, LAML, LGG, LIHC, LUAD, LUSC, MESO, OV, PAAD, PCPG, PRAD, READ, SARC, SKCM, STAD, TGCT, THCA, THYM, UCEC, UCS, UVM. Please refer to TCGA (https://tcga-data.nci.nih.gov/docs/publications/tcga/) for information about cancer type. The cancer type abbreviation table (Table 1) shows the full cancer type name.
#' @param saveFilePath (string of path to save the merged data): absolute or relative path
#' @param saveFileName (string of filename prefix)
#' @return sFileName (vector of path/filename), could be used by module B
DownloadBiospecimenClinicalData <- function(cancerType = NULL,
                                            saveFilePath =
                                                "./BiospecimenClinicalData",
                                            saveFileName = "") {
    sArchive <- "legacy"; sFieldMeta <- ""; nEntityCount <- (-1); sEndpoint <- "files"
    options(warn = -1)
    if (saveFilePath != ".") {dir.create(saveFilePath, recursive = TRUE)}
    sDir <- paste("tmp", TimeNow(), sep = "_"); dir.create(sDir)
    #  if (saveFilePath != ".") {
    #    if (!dir.exists(saveFilePath)) {
    #      dir.create(saveFilePath, recursive = TRUE)
    #    }
    #  } # dir.exists # R >= 3.2
    #  sDir <- paste("tmp", TimeNow(), sep = "_")
    #  if (!dir.exists(sDir)) {dir.create(sDir)} # dir.exists # R >= 3.2
    d <- MetaDataClin(sDir = sDir,
                      sArchive = "legacy",
                      sFieldMeta = "",
                      nEntityCount = (-1),
                      sEndpoint = "files")
    sId2Name <- d$file_name
    stopifnot(length(sId2Name) == length(unique(sId2Name))) # duplicated file_name
    names(sId2Name) <- d$file_id
    sCancerAll <- c("ACC" , "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC",
                    "ESCA", "GBM" , "HNSC", "KICH", "KIRC", "KIRP", "LAML",
                    "LGG" , "LIHC", "LUAD", "LUSC", "MESO", "OV"  , "PAAD",
                    "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT",
                    "THCA", "THYM", "UCEC", "UCS" , "UVM")
    if (is.null(cancerType)) {
        cancerType <- sCancerAll
    } else if (!all(cancerType %in% sCancerAll)) {
        print(c("cancerType should be 'NULL' (for all cancerType) or one of: ",
                sCancerAll))
        stopifnot(!all(cancerType %in% sCancerAll))
    } else {
        bCancer <- ifelse(
                          toupper(
                                  sapply(
                                         strsplit(
                                                  sapply(
                                                         strsplit(sId2Name,
                                                                  split = "\\."),
                                                         function(x){rev(x)[2]}),
                                                  split = "_"),
                                         function(y){rev(y)[1]})
                                  ) %in% cancerType,
                          TRUE, FALSE)
        sId2Name <- sId2Name[bCancer]  # cancerType
    }
    sFileNameById <- FileNameById(sId2Name,
                                  sDir = sDir,
                                  sArchive = "legacy")
    bFileRename <- file.rename(from = sFileNameById,
                               to = paste(saveFilePath,
                                          sapply(strsplit(sFileNameById,
                                                          split = "/"),
                                                 function(x){rev(x)[1]}),
                                          sep = "/"))
    stopifnot(all(bFileRename))
    if (saveFileName != "") {
        bFileRename <- FALSE
        bFileRename <- file.rename(from = paste(saveFilePath,
                                                sId2Name,
                                                sep = "/"),
                                   to = paste(saveFilePath,
                                              paste(saveFileName,
                                                    sId2Name,
                                                    sep = "__"),
                                              sep = "/"))
        stopifnot(all(bFileRename))
    }
    unlink(sDir, recursive = TRUE)
    sFileName <- paste(saveFilePath, dir(saveFilePath), sep = "/")
    options(warn = 0)
    return(sFileName)
}

#' Unique 1st column of matrix with sum the value of rows with same rowhead
#'
#' @param m0 character matrix with numbers in columns 2:end
#' @return matrix with 1st column of m0 as rownames
UniqueRow <- function(m0) {
    lm0 <- lapply(split(m0[, -1], m0[, 1]), as.numeric)
    lm1 <- lapply(lm0, matrix, ncol = ncol(m0) - 1)
    m <- do.call(rbind, lapply(lm1, colSums))
    colnames(m) <- colnames(m0)[-1]
    return(m)
}

#' Check gene symbol according to NCBI
#'
#' @param mMerged vector of "Symbol|EntrezGeneID", eg. TP53|7157.
#' @return character vector of "Symbol|EntrezGeneID" with Symbol corrected.
CheckId <- function(mMerged, s1Assay) {
    load("./SupportingFiles/lGeneInfo.rda")
    load("./SupportingFiles/dGeneHg19.rda")
    SymbolToId <- function(sSymbol) {
        #load(file = "./SupportingFiles/lGeneInfo.rda")
        lsId2Synonyms <- strsplit(lGeneInfo$sId2Synonyms, split = "\\|")
        lId <- lapply(sSymbol,
                      function(x){
                          if (x %in% lGeneInfo$sId2SymbolOfficial) {
                              names(lGeneInfo$sId2SymbolOfficial)[lGeneInfo$sId2SymbolOfficial %in% x]
                          } else if (x %in% lGeneInfo$sId2SymbolMerged) {
                              names(lGeneInfo$sId2SymbolMerged)[lGeneInfo$sId2SymbolMerged %in% x]
                          } else {
                              names(lGeneInfo$sId2Synonym)[sapply(lsId2Synonyms, function(y){x %in% y})]
                          }
                      })
        return(lId)
    }
    if (s1Assay %in% c("cna_cnv.hg18",
                       "cna_cnv.hg19",
                       "cna_nocnv.hg18",
                       "cna_nocnv.hg19")) {
        dChecked <- mMerged
    } else if (s1Assay == "somaticMutation_DNAseq") {
        mMerged[, 2] = as.character(mMerged[, 2])
        lId0 <- ifelse(mMerged[, 2] != 0,
                       ifelse(mMerged[, 2] %in% names(lGeneInfo$sGeneHist),
                              lGeneInfo$sGeneHist[mMerged[, 2]],
                              mMerged[, 2]),
                       SymbolToId(mMerged[, 1]))
        b0 <- sapply(lId0, function(x){!(length(x) < 1 | all(is.na(x)))})
        lId1 <- lId0[b0]
        mMerged <- mMerged[b0, ]
        #lChr <- lapply(lId1, function(x){lGeneInfo$sId2Chr[x]})
        lId2 <- lapply(seq(nrow(mMerged)),
                       function(x) {
                           #v = vector()
                           #for (x in seq(nrow(mMerged))) {
                           if (length(lId1[[x]]) > 1) {
                               #v <- c(v,x)
                               lId1[[x]][which(lGeneInfo$sId2Chr[x] == mMerged[x, "Chromosome"])]
                           } else {
                               lId1[[x]]
                       }})
        lId3 <- lapply(seq(nrow(mMerged)),
                       function(x) {
                           if (length(lId2[[x]]) > 1) {
                               lId2[[x]][which(paste("chr",
                                                     mMerged[x, "Chromosome"],
                                                     sep = "") == as.character(dGeneHg19[x, "chr"]) &
                                               mMerged[x, "Start_Position"] >= dGeneHg19[x, "start"] &
                                               mMerged[x,   "End_Position"] <= dGeneHg19[x,   "end"]
                                           )]
                           } else {
                               lId2[[x]]
                       }})
        b3 <- sapply(lId3, length) == 1
        lId <- lId3[b3]
        mMerged <- mMerged[b3, ]
        stopifnot(unique(sapply(lId, length)) == 1)
        sId <- unlist(lId)
        mMerged[, 2] <- sId
        mMerged[, 1] <- paste(lGeneInfo$sId2SymbolMerged[sId], sId, sep = "|")
        dChecked <- mMerged
    } else if (s1Assay %in% c("exonJunction_RNAseq")) {
        dChecked <- mMerged
    } else if (s1Assay %in% c("exon_RNAseq")) {
        dChecked <- mMerged
    } else if (s1Assay %in% c("gene_Array")) {
        dChecked <- mMerged
    } else if (s1Assay %in% c("gene.normalized_RNAseq")) {
        sMerged <- mMerged[, 1]
        lMerged <- strsplit(sMerged, split = "\\|")
        sId0 <- sapply(lMerged, function(x){x[2]})
        sId1 <- ifelse(sId0 %in% names(lGeneInfo$sGeneHist),
                       lGeneInfo$sGeneHist[sId0],
                       sId0)
        bId <- sId1 %in% names(lGeneInfo$sId2SymbolMerged)
        sId <- sId1[bId]
        sSymbolId <- paste(lGeneInfo$sId2SymbolMerged[sId],
                           sId,
                           sep = "|")
        mChecked0 <- cbind(sSymbolId,
                           mMerged[bId, 2:ncol(mMerged)])[order(sSymbolId), ]
        colnames(mChecked0) <- c("geneSymbolId",
                                 colnames(mMerged)[2:ncol(mMerged)])
        mChecked <- UniqueRow(mChecked0)
        dChecked <- as.data.frame(cbind("geneSymbolId" = rownames(mChecked),
                                        mChecked))
    } else if (s1Assay %in% c("gene_RNAseq")) {
        dChecked <- mMerged
    } else if (s1Assay %in% c("isoform.normalized_RNAseq")) {
        dChecked <- mMerged
    } else if (s1Assay %in% c("isoform_RNAseq")) {
        dChecked <- mMerged
    } else if (s1Assay %in% c("methylation_27",
                              "methylation_450")) {
        load("./SupportingFiles/sMethyProbe2geneSymbolId.rda")
        mMerged[, "Gene_Symbol"] <- sMethyProbe2geneSymbolId[mMerged[, "CpG"]]
        mMerged[, "Gene_Symbol"] <- sMethyProbe2geneSymbolId[mMerged[, "CpG"]]
        dChecked <- mMerged
    } else if (s1Assay %in% c("mir_GA.hg18",
                              "mir_GA.hg19",
                              "mir_GA.hg19.mirbase20",
                              "mir_HiSeq.hg18",
                              "mir_HiSeq.hg19",
                              "mir_HiSeq.hg19.mirbase20")) {
        # discontinued: c("hsa-mir-3118-5", "hsa-mir-3118-6", "hsa-mir-3669", "hsa-mir-3673")
        # twice: c("100302179", "100500894", "100616134", "102464837", "102465906", "102466225")
        sMerged <- mMerged[, 1]
        lMerged <- lapply(strsplit(sMerged, split = "-"),
                          function(v){c(paste(v, collapse = "-"),
                                        paste(v[-1], collapse = "-"))})
        sIdMir <- names(lGeneInfo$sId2SymbolMerged[lGeneInfo$sId2Type == "ncRNA" &
                        substr(lGeneInfo$sId2SymbolMerged, 1, 3) == "MIR"])
        lsId2Synonyms <- strsplit(lGeneInfo$sId2Synonyms[sIdMir],
                                  split = "\\|")
        lbId2Synonyms <- lapply(lMerged,
                                function(v){sapply(lsId2Synonyms,
                                                   function(v1){any(v1 %in% v)})})
        lId2Designations <- strsplit(lGeneInfo$sId2Designations[sIdMir],
                                     split = "\\|")
        lbId2Designations <- lapply(lMerged,
                                    function(v){sapply(lId2Designations,
                                                       function(v1){any(v1 %in% v)})})
        sId1 <- sapply(seq(lMerged),
                       function(n){ifelse(any(lbId2Designations[[n]]),
                                          sIdMir[lbId2Designations[[n]]],
                                          sIdMir[lbId2Synonyms[[n]]])})
        bId <- !is.na(sId1)
        sId <- sId1[bId]
        sSymbolId <- paste(lGeneInfo$sId2SymbolMerged[sId],
                           sId,
                           sep = "|")
        mChecked0 <- cbind(sSymbolId,
                           mMerged[bId, 2:ncol(mMerged)])[order(sSymbolId), ]
        colnames(mChecked0) <- c("geneSymbolId",
                                 colnames(mMerged)[2:ncol(mMerged)])
        mChecked <- UniqueRow(mChecked0)
        dChecked <- as.data.frame(cbind("geneSymbolId" = rownames(mChecked),
                                        mChecked))
    } else if (s1Assay %in% c("mirIsoform_GA.hg18",
                              "mirIsoform_GA.hg19",
                              "mirIsoform_GA.hg19.mirbase20",
                              "mirIsoform_HiSeq.hg18",
                              "mirIsoform_HiSeq.hg19",
                              "mirIsoform_HiSeq.hg19.mirbase20")) {
        dChecked <- mMerged
    } else if (s1Assay %in% c("protein_RPPA")) {
        mAb <- mMerged[, -1]
        mode(mAb) <- "numeric"
        rownames(mAb) <- sapply(strsplit(mMerged[, 1], split = "\\|"), function(v){v[2]})
        lAbGene <- sapply(strsplit(mMerged[, 1], split = "\\|"), function(v){strsplit(v[1], split = " ")})
        sIdProtein <- names(lGeneInfo$sId2SymbolMerged[lGeneInfo$sId2Type ==
                            "protein-coding"])
        lbId2SymbolOfficial <-
            lapply(lAbGene,
                   function(v){sapply(lGeneInfo$sId2SymbolOfficial[sIdProtein],
                                      function(s){any(s %in% v)})})
        lsId2SymbolOfficial <-
            lapply(lbId2SymbolOfficial,
                   function(v){names(lGeneInfo$sId2SymbolOfficial[sIdProtein])[v]})
        names(lsId2SymbolOfficial) <- rownames(mAb)
        sId2Ab <- unlist(sapply(names(lsId2SymbolOfficial),
                                function(s){paste(lsId2SymbolOfficial[[s]],
                                                  s,
                                                  sep = "|")}))
        sId <- sapply(strsplit(sId2Ab, split = "\\|"), function(v){v[1]})
        sAb <- sapply(strsplit(sId2Ab, split = "\\|"), function(v){v[2]})
        sSymbolId <- paste(lGeneInfo$sId2SymbolMerged[sId],
                           sId,
                           sAb,
                           sep = "|")
        mChecked <- mAb[sAb,]
        rownames(mChecked) <- sSymbolId
        mChecked <- mChecked[order(rownames(mChecked)), ]
        dChecked <- as.data.frame(cbind("geneSymbolIdTypeNameAntibody" = rownames(mChecked),
                                        mChecked))
    } else if (s1Assay %in% c("proteome_iTRAQ")) {
        dChecked <- mMerged
    }
    return(dChecked)
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

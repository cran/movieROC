\name{HCC}
\alias{HCC}
\docType{data}
\title{
Hepatocellular carcinoma data
}
\description{
This dataset is derived from gene expression arrays of tumor and adjacent non-tumor tissues of 62 Taiwanese cases of hepatocellular carcinoma. The complete dataset was
deposited in \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi}{NCBI’s Gene Expression Omnibus (GEO)} and it is available through series accession number GSE37988. This dataset contains 948 from the 27,578 autosomal CpG sites screened.
}
\usage{data("HCC")}
\format{
  A data frame with 124 observations on 952 variables. First 4 variables are \code{tissue} (identification number for the tissue; from 1 to 62), \code{sex} (female or male), \code{age} (age in years of the patient), and \code{tumor} (status of the tissue; nontumor or tumor). The following 948 from \code{cg03409548} to \code{cg20240860} are numeric variables containing the relative gene expression intensities of the corresponding gene.
}
\source{
NCBI’s Gene Expression Omnibus (GEO) - Series accession number GSE37988 (\url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37988})
}
\examples{
data(HCC)
str(HCC) # Data structure
table(HCC$tumor) # Number of non-tumor and tumor tissues

# Histograms of gene 18384097 expression intensities for non-tumor and tumor tissues
oldpar <- par(mfrow = c(2,1))
x <- subset(HCC, tumor == "nontumor")$cg18384097
y <- subset(HCC, tumor == "tumor")$cg18384097
hist(x, xlim = c(0,1), main = "Gene 18384097 in non-tumor tissues")
hist(y, xlim = c(0,1), main = "Gene 18384097 in tumor tissues")
par(oldpar)
}
\keyword{datasets}

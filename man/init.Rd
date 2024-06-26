% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/init.R
\name{init}
\alias{init}
\title{Create a dataset.}
\usage{
init(
  rna,
  atac,
  segmentation,
  rna_normalisation_factors = rna \%>\% auto_normalisation_factor(),
  atac_normalisation_factors = atac \%>\% auto_normalisation_factor(),
  rna_likelihood = "NB",
  atac_likelihood = "NB",
  reference_genome = "GRCh38",
  description = "(R)CONGAS+ model",
  smooth = FALSE,
  multiome = FALSE
)
}
\arguments{
\item{rna}{A tibble with single-cell RNA data.}

\item{atac}{A tibble with single-cell ATAC data.}

\item{segmentation}{A tibble with the input segmentation.}

\item{rna_normalisation_factors}{The RNA tibble with the input per-cell normalisation factors.
By default these are computed by function \code{auto_normalisation_factor}.}

\item{atac_normalisation_factors}{The ATAC tibble with the input per-cell normalisation factors.
By default these are computed by function \code{auto_normalisation_factor}.}

\item{rna_likelihood}{Type of likelihood used for RNA data (\code{"G"} for Gaussian and
\code{""NB} for Negative Binomial). The RNA default is \code{"G"}.}

\item{atac_likelihood}{Type of likelihood used for ATAC data, with default \code{"NB"}.}

\item{reference_genome}{Either \code{"GRCh38"} or \code{"hg19"}.}

\item{description}{A model in-words description.}

\item{smooth}{If yes, input segments are smootheed by joining per chromosome segments that
have the same ploidy.}

\item{mutiome}{Default to FALSE. Flag indicating whether the RNA and ATAC observations are the result of a matched RNA-ATAC sequencing assay such as 10x multiome assay.
(i.e., there is a 1:1 correspondence between barcodes of the two modalities.)}
}
\value{
An object of class \code{rcongasplus}
}
\description{
This function creates a dataset (an object of class \code{rcongasplus}) by assembling multiple single-cell input measurements
(ATAC and/or RNA data modalities), the input segmentation (from bulk DNA sequencing),
and the per-cell normalisation factors for the data.

All input data are passed as tibbles; the input formats are as follows:
\itemize{
\item for single-cell ATAC/RNA data, the \code{cell} identifier, the genomic coordinates
(\code{chr}, \code{from}, \code{to}) which refer either to an ATAC peak, or an RNA gene
identifier, and a \code{value} reporting the reads mapped.
\item for the input segmentation, the genomic coordinates
(\code{chr}, \code{from}, \code{to}) which refer to the segment, and the number of
\code{copies} (i.e., DNA ploidy) of the segment.
\item for normalization factors the \code{cell} identifier, the actual \code{normalisation_factor}
and the \code{modality} to wihch the factor refers to
}

This function receives also other parameters - e.g., the models likelihoods - which
will determine the overall behaviour of the underlying model, and how data are preared for inference.
\itemize{
\item A Negative Binomial likelihood (\code{"NB"}), which works directly from raw counts data
\item A Gaussian likelihood (\code{"G"}), which requires a z-score transformation of the data. This consists
in :
\itemize{
\item scaling raw counts by the input normalization factors;
\item computing z-scores per cell;
\item summing up z-scores per segment;
\item computing z-scores per segment;
\item center the z-scores mean to the input ploidy.
}
}
}
\examples{
data("example_input")

# For instance, RNA data
example_input$x_rna \%>\% print

# .. or ATAC data
example_input$x_atac \%>\% print

# .. and segmentation
example_input$x_segmentation \%>\% print

# .. and normalisation factors can be computed (default)
example_input$x_rna \%>\% auto_normalisation_factor()

x = init(
  rna = example_input$x_rna,
  atac = example_input$x_atac,
  segmentation = example_input$x_segmentation,
  rna_likelihood = "G",
  atac_likelihood = 'NB',
  description = 'My model')

print(x)
}

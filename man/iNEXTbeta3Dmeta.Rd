% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Mainfunction.R
\name{iNEXTbeta3Dmeta}
\alias{iNEXTbeta3Dmeta}
\title{Estimates the difference of standardized 3D diversity with common sample coverage (for alpha, beta, gamma diversity, and four classes of dissimilarity measures) between two treatments, and fit a fixed- or random-effects model to perform meta analysis}
\usage{
iNEXTbeta3Dmeta(
  data,
  model = "FE",
  diversity = "TD",
  order.q = 0,
  datatype = "abundance",
  level = NULL,
  nboot,
  treatment_order,
  conf = 0.95,
  PDtree,
  PDreftime = NULL,
  PDtype = "meanPD",
  FDdistM,
  FDtype = "AUC",
  FDtau = NULL,
  FDcut_number = 30
)
}
\arguments{
\item{data}{(a) For \code{datatype = "abundance"}, data can be input as a data.frame. The data frame has study/site and treatment as the first two columns, followed by columns for species names. Here an assemblage refers to a combination of study/site and treatment. \cr
(b) For \code{datatype = "incidence_raw"}, data can be input as a data.frame. The data frame has study/site, treatment, patch as the first three columns, followed by columns for species names.}

\item{model}{Selection of model type: 'FE' = Fixed-effects model, 'RE' = Random-effects model.}

\item{diversity}{Selection of diversity type: 'TD' = Taxonomic diversity, 'PD' = Phylogenetic diversity, and 'FD' = Functional diversity.}

\item{order.q}{A numerical value specifying the diversity order, Default is \code{q = 0, 1, 2}.}

\item{datatype}{Data type of input data: individual-based abundance data (\code{datatype = "abundance"}) or species by sampling-units incidence data (\code{datatype = "incidence_raw"}) with all entries being 0 (non-detection) or 1 (detection).}

\item{level}{A numerical value between 0 and 1 specifying the sample coverage level used for computing standardized diversity and dissimilarity. By default(\code{level = NULL}), the function automatically calculates standardized 3D gamma, alpha, and beta diversities, along with four dissimilarity indices, up to the minimum coverage achieved by doubling the reference sample size across all site and treatment combinations.}

\item{nboot}{A positive integer specifying the number of bootstrap replications when assessing sampling uncertainty for estimating standardized beta3D diversity and the associated confidence intervals. Default is 10. If more accurate results are required, set \code{nboot = 100} (or \code{nboot = 200}).}

\item{treatment_order}{A character vector for the names of treatment. The difference of standardized 3D diversity will be computed as diversity of the first treatment minus the diversity of second treatment.}

\item{conf}{A positive number < 1 specifying the level of confidence interval. Default is 0.95.}

\item{PDtree}{(required only when \code{diversity = "PD"}), a phylogenetic tree in Newick format for all observed species in the pooled data.}

\item{PDreftime}{(required only when \code{diversity = "PD"}), a numerical value specifying reference times for PD. Default is NULL (i.e., the age of the root of PDtree).}

\item{PDtype}{(argument only for \code{diversity = "PD"}), select PD type: \code{PDtype = "PD"} (effective total branch length) or \code{PDtype = "meanPD"} (effective number of equally divergent lineages). Default is \code{PDtype = "meanPD"}, where \code{meanPD = PD/tree depth}.}

\item{FDdistM}{(required only when \code{diversity = "FD"}), a species pairwise distance matrix for all species in the pooled data.}

\item{FDtype}{(required only when \code{diversity = "FD"}), select FD type: \code{FDtype = "tau_value"} for FD under a specified threshold value, or \code{FDtype = "AUC"} (area under the curve of tau-profile) for an overall FD which integrates all threshold values between zero and one. Default is \code{FDtype = "AUC"}.}

\item{FDtau}{(required only when \code{diversity = "FD"} and \code{FDtype = "tau_value"}), a numerical value between 0 and 1 specifying the tau value (threshold level) that will be used to compute FD. If \code{FDtau = "NULL"} (default), then threshold is set to be the mean distance between any two individuals randomly selected from the pooled data (i.e., quadratic entropy).}

\item{FDcut_number}{(required only when \code{diversity = "FD"} and \code{FDtype = "AUC"}) a numeric number to cut [0, 1] interval into equal-spaced sub-intervals to obtain the AUC value by integrating the tau-profile. Equivalently, the number of tau values that will be considered to compute the integrated AUC value. Default is 30. A larger value can be set to obtain more accurate AUC value.}
}
\value{
The function returns a list containing seven components: Gamma, Alpha, Beta, 1-C, 1-U, 1-V, and 1-S. Each component consists of two parts. The first part is a dataframe that includes a column named Site (or the user-defined name in the input data), representing the study or site identity. It also contains a column Difference, which indicates the difference in diversity between the two treatments,
along with SE for the standard error of the difference, LCL and UCL for the lower and upper confidence limits, and Order.q for the diversity order q. The type of diversity measure is recorded in the Diversity column, which can be TD, PD, or FD. Additionally, the dataframe includes two columns showing the estimated diversity values of the two treatments for each site,
and a Weight column representing the weight assigned to each site for the fixed or random effect model. The second part of each component is a summary table reporting meta-analytic statistics under the fixed or random effect model, including Cochran's Q statistic (Q_val), the degrees of freedom (df_val), the associated p-value (p_val) for the heterogeneity test, the heterogeneity percentage (I2_val), and the estimated between-site variance (tau2_val).
}
\description{
\code{iNEXTbeta3Dmeta} is a function that estimates the difference of standardized 3D (taxonomic, phylogenetic and functional) beta diversity between two treatments (e.g., enhanced vs. control), and perform meta analysis (fixed or random effect model) for several studies/sites.
}
\examples{

## Taxonomic diversity for abundance data

# Coverage-based standardized TD
data("Spider_abundance_data")
output1_abu <- iNEXTbeta3Dmeta(data = Spider_abundance_data, model = "RE", diversity = "TD",
                               order.q = c(0, 1, 2), datatype = "abundance", level = NULL,
                               nboot = 10, treatment_order = c("Enhanced", "Control"), conf = 0.95)

output1_abu

## Phylogenetic diversity for abundance data

# Coverage-based standardized PD
data("Spider_abundance_data")
data("Spider_tree")
output2_abu <- iNEXTbeta3Dmeta(data = Spider_abundance_data, model = "RE", diversity = "PD",
                               order.q = c(0, 1, 2), datatype = "abundance", level = NULL,
                               nboot = 10, treatment_order = c("Enhanced", "Control"), conf = 0.95,
                               PDtree = Spider_tree, PDreftime = NULL, PDtype = "meanPD")
output2_abu

## Functional diversity for abundance data

# Coverage-based standardized FD
data("Spider_abundance_data")
data("Spider_distM")
output3_abu <- iNEXTbeta3Dmeta(data = Spider_abundance_data, model = "RE", diversity = "FD",
                               order.q = c(0, 1, 2), datatype = "abundance", level = NULL,
                               nboot = 10, treatment_order = c("Enhanced", "Control"), conf = 0.95,
                               FDdistM = Spider_distM, FDtype = "AUC", FDcut_number = 30)
output3_abu

## Taxonomic diversity for incidence data

# Coverage-based standardized TD
data("Bat_incidence_data")
output1_inc <- iNEXTbeta3Dmeta(data = Bat_incidence_data, model = "RE", diversity = "TD",
                               order.q = c(0, 1, 2), datatype = "incidence_raw", level = NULL,
                               nboot = 10, treatment_order = c("Enhanced", "Control"), conf = 0.95)
output1_inc

## Phylogenetic diversity for incidence data

# Coverage-based standardized PD
data("Bat_incidence_data")
data("Bat_tree")
output2_inc <- iNEXTbeta3Dmeta(data = Bat_incidence_data, model = "RE", diversity = "PD",
                               order.q = c(0, 1, 2), datatype = "incidence_raw", level = NULL,
                               nboot = 10, treatment_order = c("Enhanced", "Control"), conf = 0.95,
                               PDtree = Bat_tree, PDreftime = NULL, PDtype = "meanPD")
output2_inc

## Functional diversity for incidence data

# Coverage-based standardized FD
data("Bat_incidence_data")
data("Bat_distM")
output3_inc <- iNEXTbeta3Dmeta(data = Bat_incidence_data, model = "RE", diversity = "FD",
                               order.q = c(0, 1, 2), datatype = "incidence_raw", level = NULL,
                               nboot = 10, treatment_order = c("Enhanced", "Control"), conf = 0.95,
                               FDdistM = Bat_distM, FDtype = "AUC", FDcut_number = 30)
output3_inc

}

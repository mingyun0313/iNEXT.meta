\name{Spider_distM}
\alias{Spider_distM}
\docType{data}
\title{
  Species pairwise distance matrix for spider species
}
\description{
  This dataset includes a spider species pairwise symmetric distance matrix for 199 species listed in the \code{Spider_abundance_data} dataset. Each element in the matrix represents a Gower distance between two species computed from species traits. \cr
}
\usage{data("Spider_distM")}
\format{
  \code{Spider_distM} is a 199 x 199 (species by species) symmetric data.frame. Each element of the data.frame is between zero (for species pairs with identical traits) and one. \cr
}

\keyword{datasets}

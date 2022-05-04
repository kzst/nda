\name{CWTS_2020}
\alias{CWTS_2020}
\docType{data}
\title{
CWTS Leiden's 2020 dataset, where the data frame has 1176 observations of 42 variables.
}
\description{
Sample datasets for Network-based, non-parametric dimensionality reduction and analysis (NDA)

CWTS Leiden's 2020 dataset, where the data frame has 1176 observations of 42 variables.
}
\usage{data("CWTS_2020")}
\format{
  A data frame with 1176 observations of 42 variables.
}

\source{
CWTS Leiden Ranking 2020: https://www.leidenranking.com/ranking/2020/list
}

\examples{
data(CWTS_2020)

%%% Convert Unicodes to UTF-8 in the university names
% library(stringi)
% rownames(CWTS_2020)<-stri_unescape_unicode((gsub("<U\\+(....)>", "\\\\u\\1",(rownames(CWTS_2020)))))

}
\keyword{datasets}

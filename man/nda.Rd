\name{nda}
\alias{nda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Network-based, non-parametric dimensionality reduction and analysis
%%  ~~function to do ... ~~
}
\description{
Network-based dimensionality reduction of high-dimensional, low-sample-size datasets
}
\usage{
nda(data,cor_method=1,min_R=0,min_comm=2,Gamma=1,null_mode_type=4,mod_mode=6,min_evalue=0,min_communality=0,com_communalities=0,use_rotation=FALSE)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A numeric data frame}
  \item{cor_method}{Correlation method (optional). 
'1' Pearson's correlation (default),
'2' Spearman's correlation, '3' Kendall's correlation,
'4' Distance correlation}
  \item{min_R}{Minimal square correlation between indicators (default: 0)}
  \item{min_comm}{Minimal number of indicators in a community (default: 2)}
  \item{Gamma}{Gamma parameter in multiresolution null_modell (default: 1)}
  \item{null_modell_type}{'1' Differencial Newmann-Grivan's null modell, '2' The null model is the mean of square correlations between indicators, '3' The null modell is the specified minimal square correlation, '4' Newmann-Grivan's modell (default)}
  \item{min_evalue}{Minimal eigenvector centrality value (default: 0)}
  \item{min_communality}{Minimal communality value of indicators (default: 0)}
  \item{com_communalities}{Minimal common communalities (default: 0)}
  \item{use_rotation}{FALSE no rotation (default), TRUE varimax rotation}
}
\details{
NDA both works on low and high simple size datasets. If min_evalue=min_communality=com_communalities=0 than there is no feature selection.
}
\value{

\item{communality}{Communality estimates for each item. These are merely the sum of squared factor loadings for that item.}
\item{loadings}{A standard loading matrix of class “loadings"}
\item{uniqueness}{Uniqueness value of indicators.}
\item{factors}{Number of found factors}
\item{scores}{Estimates of the factor scores are reported}
\item{n.obs}{Number of observations specified or found}
\item{fn}{NDA}

\references{
%% ~put references to the literature/web site here ~
}
\author{
Prof. Zsolt T. Kosztyan, Department of Quantitative Methods, Institute of Management, Faculty of Business and Economics, University of Pannonia, Hungary

e-mail: kzst@gtk.uni-pannon.hu
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

data(swiss)
df<-swiss
p<-nda(df)
p

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):

\keyword{multivariate}
\keyword{nonparametric}
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{dimensionality reduction}
% \concept{ ~cpt2 }
% Use only one concept per line.
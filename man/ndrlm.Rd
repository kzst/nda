\name{ndrlm}
\alias{ndrlm}
\title{Genearlized Network-based Dimensionality Reduction and Regression (GNDR)
}
\description{
The main function of Generalized Network-based Dimensionality Reduction and Regression (GNDR) for supervised learning.
}
\usage{
ndrlm(Y,X,latents="in",dircon=FALSE,optimize=TRUE,
                target="adj.r.square",rel_weight=FALSE,
                cor_method=1,
                cor_type=1,min_comm=2,Gamma=1,
                null_model_type=4,mod_mode=1,use_rotation=FALSE,
                rotation="oblimin",pareto=FALSE,fit_weights=NULL,
                lower.bounds.x = c(rep(-100,ncol(X))),
                upper.bounds.x = c(rep(100,ncol(X))),
                lower.bounds.latentx = c(0,0,0,0),
                upper.bounds.latentx = c(0.6,0.6,0.6,0.3),
                lower.bounds.y = c(rep(-100,ncol(Y))),
                upper.bounds.y = c(rep(100,ncol(Y))),
                lower.bounds.latenty = c(0,0,0,0),
                upper.bounds.latenty = c(0.6,0.6,0.6,0.3),
                popsize = 20, generations = 30, cprob = 0.7, cdist = 5,
                mprob = 0.2, mdist=10, seed=NULL)

}
\arguments{
  \item{Y}{A numeric data frame of output variables}
  \item{X}{A numeric data frame of input variables}
  \item{latents}{The employs of latent variables: "in" employs latent-independent variables (default); "out" employs latent-dependent variables; "both" employs both latent-dependent and latent independent variables; "none" do not employs latent variable (= multiple regression)}
  \item{dircon}{Wether enable or disable direct connection between input and output variables (default=FALSE)}
  \item{optimize}{Optimization of fittings (default=TRUE)}
  \item{target}{Target performance measures. The possible target measure are "adj.r.square" = adjusted R  square (default), "r.sqauare" = R square, "MAE" = mean absolute error, "MAPE" = mean absolute percentage error, "MASE" = mean absolute scaled error ,"MSE"= mean square error,"RMSE" = root mean square error}
  \item{rel_weight}{Use relative weights. In this case, all weights should be non-negative. (default=FALSE)}
  \item{cor_method}{Correlation method (optional).
'1' Pearson's correlation (default),
'2' Spearman's correlation, '3' Kendall's correlation,
'4' Distance correlation}
  \item{cor_type}{Correlation type (optional).
'1' Bivariate correlation (default),
'2' partial correlation, '3' semi-partial correlation}
  \item{min_comm}{Minimal number of indicators per community (default: 2).}
  \item{Gamma}{Gamma parameter in multiresolution null modell (default: 1).}
  \item{null_model_type}{'1' Differential Newmann-Grivan's null model, '2' The null model is the mean of square correlations between indicators, '3' The null model is the specified minimal square correlation, '4' Newmann-Grivan's modell (default)}
  \item{mod_mode}{Community-based modularity calculation mode: '1' Louvain modularity (default), '2' Fast-greedy modularity, '3' Leading Eigen modularity, '4' Infomap modularity, '5' Walktrap modularity, '6' Leiden modularity}
  \item{use_rotation}{FALSE no rotation (default), TRUE the rotation is used.}
  \item{rotation}{"none", "varimax", "quartimax", "promax", "oblimin", "simplimax", and "cluster" are possible rotations/transformations of the solution. "oblimin" is the default, if use_rotation is TRUE.}
  \item{pareto}{in the case of multiple objectives TRUE (default value) provides pareto-optimal solution, while FALSE provides weighted mean of objective functions (see out_weights)}
  \item{fit_weights}{weights of fitting the output variables (weights of means of objectives)}
  \item{lower.bounds.x}{Lower bounds of weights of independent variables in GNDA}
  \item{upper.bounds.x}{Upper bounds of weights of independent variables in GNDA}
  \item{lower.bounds.latentx}{Lower bounds of hyper-parementers of GNDA for independent variables (values must be positive)}
  \item{upper.bounds.latentx}{Upper bounds of hyper-parementers of GNDA for independent variables (value must be lower than one)}
  \item{lower.bounds.y}{Lower bounds of weights of dependent variables in GNDA}
  \item{upper.bounds.y}{Upper bounds of weights of dependent variables in GNDA}
  \item{lower.bounds.latenty}{Lower bounds of hyper-parementers of GNDA for dependent variables (values must be positive)}
  \item{upper.bounds.latenty}{Upper bounds of hyper-parementers of GNDA for dependent variables (value must be lower than one)}
  \item{popsize}{size of population of NSGA-II for fitting betas (default=20)}
  \item{generations}{number of generations to breed of NSGA-II for fitting betas (default=30)}
  \item{cprob}{crossover probability of NSGA-II for fitting betas (default=0.7)}
  \item{cdist}{crossover distribution index of NSGA-II for fitting betas (default=5)}
  \item{mprob}{mutation probability of NSGA-II for fitting betas (default=0.2)}
  \item{mdist}{mutation distribution index of NSGA-II for fitting betas (default=10)}
  \item{seed}{default seed value (default=NULL, no seed)}
}
\details{
  NDRLM is a variable fitting with feature selection based on the tunes of GNDA method with NSGA-II algorithm for parameter fittings.
}
\value{
\item{fval}{Objective function for fitting}
\item{target}{Target performance measures. The possible target measure are "adj.r.square" = adjusted R  square (default), "r.sqauare" = R square, "MAE" = mean absolute error, "MAPE" = mean absolute percentage error, "MASE" = mean absolute scaled error ,"MSE"= mean square error,"RMSE" = root mean square error}
\item{hyperparams}{optimized hyperparameters}
\item{pareto}{in the case of multiple objectives TRUE provides pareto-optimal solution, while FALSE (default) provides weighted mean of objective functions (see out_weights)}
\item{Y}{A numeric data frame of output variables}
\item{X}{A numeric data frame of input variables}
\item{latents}{Latent model: "in", "out", "both", "none"}

\item{NDAin}{GNDA object, which is the result of model reduction and features selection in the case of employing latent-independent variables}
\item{NDAin_weight}{Weights of input variables (used in \code{\link{ndr}})}
\item{NDAin_min_evalue}{Optimized minimal eigenvector centrality value (used in \code{\link{ndr}})}
\item{NDAin_min_communality}{Optimized minimal communality value of indicators (used in \code{\link{ndr}})}
\item{NDAin_com_communalities}{Optimized
minimal common communalities (used in \code{\link{ndr}})}
\item{NDAin_min_R}{Optimized
minimal square correlation between indicators (used in \code{\link{ndr}})}
\item{NDAout}{GNDA object, which is the result of model reduction and features selection in the case of employing latent-dependent variables}
\item{NDAout_weight}{Weights of input variables (used in \code{\link{ndr}})}
\item{NDAout_min_evalue}{Optimized minimal eigenvector centrality value (used in \code{\link{ndr}})}
\item{NDAout_min_communality}{Optimized minimal communality value of indicators (used in \code{\link{ndr}})}
\item{NDAout_com_communalities}{Optimized
minimal common communalities (used in \code{\link{ndr}})}
\item{NDAout_min_R}{Optimized
minimal square correlation between indicators (used in \code{\link{ndr}})}
\item{fits}{List of linear regrassion models}
\item{otimized}{Wheter fittings are optimized or not}
\item{NSGA}{Outpot structure of NSGA-II optimization (list), if the optimization value is true (see in \code{\link[mco:nsga2]{mco::nsga2}})}
\item{extra_vars.X}{Logic variable. If direct connection (dircon=TRUE) is allowed not only the latent but the excluded input variables are analyized in the linear models as extra input variables.}
\item{extra_vars.Y}{Logic variable. If direct connection (dircon=TRUE) is allowed not only the latent but the excluded output variables are analyized in the linear models as extra input variables.}
\item{dircon_X}{The list of input variables which are directly connected to output variables.}
\item{dircon_Y}{The list of output variables which are directly connected to output variables.}
\item{seed}{applied seed value (default=NULL, no seed)}
\item{fn}{Function (regression) name: NDRLM}
\item{Call}{Callback function}
}

\references{
Kosztyan, Z. T., Kurbucz, M. T., & Katona, A. I. (2022). Network-based dimensionality reduction of high-dimensional, low-sample-size datasets. Knowledge-Based Systems, 109180. doi:10.1016/j.knosys.2022.109180
}

\author{
Zsolt T. Kosztyan*, Marcell T. Kurbucz, Attila I. Katona

e-mail*: kosztyan.zsolt@gtk.uni-pannon.hu
}

\seealso{\code{\link{ndr}}, \code{\link{plot}}, \code{\link{summary}}, \code{\link[mco:nsga2]{mco::nsga2}}.}

\examples{

# Using NDRLM without fitting optimization
X<-freeny.x
Y<-freeny.y
NDRLM<-ndrlm(Y,X,optimize=FALSE)
summary(NDRLM)
plot(NDRLM)

\dontrun{
# Using NDRLM with optimized fitting

NDRLM<-ndrlm(Y,X)
summary(NDRLM)

# Using Leiden's modularity for grouping variables

X<-freeny.x
Y<-freeny.y
NDRLM<-ndrlm(Y,X,mod_mode=6)
plot(NDRLM)

# Using relative weights

NDRLM<-ndrlm(Y,X,mod_mode=6,rel_weight=TRUE)
plot(NDRLM)

# Using Spearman's correlation

NDRLM<-ndrlm(Y,X,cor_method=2)
summary(NDRLM)

# Using greater population and generations

NDRLM<-ndrlm(Y,X,popsize=52,generations=40)
summary(NDRLM)

# No latent variables
NDRLM<-ndrlm(Y,X,latents="none")
plot(NDRLM)

# In-out model
library(lavaan)
df<-PoliticalDemocracy # Data of Political Democracy

dem<-PoliticalDemocracy[,c(1:8)]
ind60<-PoliticalDemocracy[,-c(1:8)]

NBSEM<-ndrlm(dem,ind60,latents = "both",seed = 2)
plot(NBSEM)
}


}


\keyword{multivariate}
\keyword{nonparametric}
\concept{dimensionality}
\concept{reduction}

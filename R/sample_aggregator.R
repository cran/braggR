#' Revealed Aggregator
#'
#' This function allows the user to compute the revealed aggregator from \emph{Satopää, V.A. (2021):
#' Regularized Aggregation of One-off Probability Predictions}.
#' The current version of the paper is available at \url{https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3769945}.
#' @param p Vector of \eqn{K \ge 2} forecasters' probability estimates of a future binary event.
#' These values represent probability predictions and must be strictly between 0 and 1.
#' @param p0 The forecasters' common prior. This represents a probability prediction based 
#' on some of the forecasters' common evidence and must be strictly between 0 and 1.
#' @param alpha,beta The shape and scale parameters of the prior beta distribution of the common prior.
#' If omitted, the sampler uses the fixed common prior given by \code{p0}.
#' However, if \code{alpha} and \code{beta} are provided, they must be strictly positive.
#' In this case, the common prior \code{p0} will be treated as a random variable and sampled along with the other model parameters.
#' @param a,b The parameters for the prior distribution of \eqn{(\rho, \gamma, \delta)} in  \emph{Satopää, V.A. (2021)}.
#' The default choice \code{a = 1/2} and \code{b = 1/2} gives the Jeffreys' independence prior.
#' If \eqn{p_0} is not equal to \eqn{0.5}, then \code{a = 1} and \code{b = 1/2} give the Jeffreys' prior.
#' @param num_sample The number of posterior samples to be drawn.
#' This does not take into account burnin and thinning.
#' @param burnin The number of the initial \code{num_sample} posterior draws that are discarded for burnin.
#' This value cannot exceed \code{num_sample}.
#' @param thin After \code{burnin} draws have been discarded, the final sample is formed by keeping every \code{thin}'th value.
#' To ensure that the final sample holds at least two draws, \code{thin} can be at most \code{(num_sample-burnin)/2}.
#' @param seed The seed value for random value generation.
#' @return A data frame with rows representing posterior draws of \eqn{(p*, \rho, \gamma, \delta, p0)}. The columns are:
#' \itemize{
##'  \item{\code{aggregate}: The posterior samples of the oracle aggregator \eqn{p*}.
##'  The average of these values gives the revealed aggregator \eqn{p''}.
##'  The 95\% interval of these values gives the 95\% credible interval of the oracle aggregator.}
##'  \item{\code{rho}: The posterior samples of the forecasters' shared evidence, \eqn{\rho}.}
##'  \item{\code{gamma}: The posterior samples of the forecasters' total evidence, \eqn{\gamma}.
##'  The difference \code{gamma}-\code{rho} gives the posterior samples of
##'  the forecasters' rational disagreement.}
##'  \item{\code{delta}: The posterior samples of the forecasters' total evidence plus noise, \eqn{\delta}.
##'  The difference \code{delta}-\code{gamma} gives the posterior samples of
##'  the forecasters' irrational disagreement.}
##'  \item{\code{p0}: The posterior samples of the forecasters' common prior.
##'  If a beta prior distribution is not specified via the arguments \code{alpha} and \code{beta},
##'  then all elements of this column are equal to the fixed common prior given by the \code{p0} argument.}
##' }
#' @examples
#' # Illustration on Scenario B in Satopää, V.A. (2021).
#' # Forecasters' probability predictions:
#' p = c(1/2, 5/16, 1/8, 1/4, 1/2)
#'
#' # Aggregate with a fixed common prior of 0.5.
#' # Sample the posterior distribution:
#' post_sample = sample_aggregator(p, p0 = 0.5, num_sample = 10^6, seed = 1)
#' # The posterior means of the model parameters:
#' colMeans(post_sample[,-1])
#' # The posterior mean of the oracle aggregator, a.k.a., the revealed aggregator:
#' mean(post_sample[,1])
#' # The 95% credible interval for the oracle aggregator:
#' quantile(post_sample[,1], c(0.025, 0.975))
#'
#'
#' # Aggregate based a uniform distribution on the common prior
#' # Recall that Beta(1,1) corresponds to the uniform distribution.
#' # Sample the posterior distribution:
#' post_sample = sample_aggregator(p, alpha = 1, beta = 1, num_sample = 10^6, seed = 1)
#' # The posterior means of the oracle aggregate and the model parameters:
#' colMeans(post_sample)
#' @export

#' @useDynLib braggR, .registration = TRUE
#' @importFrom Rcpp sourceCpp

sample_aggregator <- function(p, p0 = NULL,
                              alpha = NULL, beta = NULL,
                              a = 1/2, b = 1/2,
                              num_sample = 1000000, burnin = num_sample/2, thin = 1,
                              seed = 1) {
  set.seed(seed)
  beta_distr = !is.null(alpha) & !is.null(beta)

  ########################################
  # Check the input:
  # Either p0 or (alpha,beta) must be specified:
  if(is.null(p0) & (is.null(alpha) | is.null(beta))) stop("Either p0 or (alpha,beta) must be specified!")
  
  # All values of p and p0 must be non-NAs:
  if(any(is.na(c(p, p0)))) stop("Predictions and common prior cannot be NAs!")

  # All values of p must be within (0,1):
  if(any( (p < 0) | (p > 1))) stop("All predictions must be within (0,1)!")

  if(!beta_distr){
    # Check that provided p0 value is within (0,1):
    if((p0<0) | (p0>1)) stop("Common prior must be within (0,1)!")
  } else {
    # Check that alpha and beta are positive numbers:
    if(is.na(alpha) | is.na(beta)) stop("Alpha or beta cannot be NAs!")
    if((alpha<=0) | (beta<=0)) stop("Alpha and beta must be positive numbers!")
  }

  # num_sample, burnin, and thin must be positive numbers:
  if((num_sample <= 0) || is.na(num_sample) || is.null(num_sample)) stop("num_sample must be a positive number!")
  if((burnin <= 0) || is.na(burnin) || is.null(burnin)) stop("burnin must be a positive number!")
  if((thin <= 0) || is.na(thin) || is.null(thin)) stop("thin must be a positive number!")

  # burnin cannot be larger than num_sample:
  if(burnin > num_sample) stop("num_sample must be larger than burnin!")

  # At least two draws must remain after burnin and thinning:
  if(2*thin >= (num_sample-burnin)) stop("thin can be at most (num_sample-burnin)/2!")

  # Therefore should be at least 2 predictions:
  if(length(p) == 1) stop("Aggregation needs at least 2 predictions!")

  # Both a and b must be positive:
  if(a <= 0 || b <= 0) stop("a and b must be positive!")
  ########################################


  ########################################
  # Sample the posterior distribution:
  if(beta_distr){
    post_sample = .sample_parameters_cpp(p = p,
                                        alpha = alpha, beta = beta,
                                        a = a, b = b,
                                        num_sample = num_sample,
                                        seed = seed)
  } else {
    post_sample = .sample_parameters_cpp(p = p,
                                        p0 = p0,
                                        alpha = -1, beta = -1,
                                        a = a, b = b,
                                        num_sample = num_sample,
                                        seed = seed)
  }
  ########################################

  ########################################
  # Process the sampled values:
  post_sample = post_sample[-c(1:burnin),]                   ## burnin
  post_sample = post_sample[seq(1,nrow(post_sample), thin),] ## thinning
  post_sample = data.frame(post_sample)
  names(post_sample) = c("aggregate", "rho", "gamma", "delta", "p0")
  return(post_sample)
}

#' One dimensional Monte-Carlo integration
#'
#' @param g the known function to integrate over a given interval
#' @param a (finite) left endpoint of the given interval
#' @param b (finite) right endpoint of the given interval
#' @param naive boolean whether to naive estimate or efficient estimator, see details.
#' @param n number of variates to simulate
#' @param sig_lvl significance level of confidence intervals
#' @param ... additional arguments to pass to the function \code{g}; must be named
#'
#' @description {One dimensional integrals over an interval can be approximated
#' via sample means of (functions of) uniformly distributed random variables, via
#' the SLLN and error estimates can be obtained via CLT. This naive approach is
#' not ideal for multi-dimensional integrals, however. The integrand must be known.}
#' @details {The naive estimate computes the sample mean of \eqn{g(U)} over an IID sample
#' of uniformly distributed RVs, \eqn{U_1,...,U_n} while
#' the efficient estimator uses the sample mean of \eqn{0.5(g(U)+g(a+b-U))}, which
#' can be shown to have a lower variance than the aforementioned.}
#' @return data.frame containing
#' \itemize{
#' \item \code{estimate} the point estimate of the integral,
#' \item \code{lb} the lower bound of the confidence interval,
#' \item \code{ub} the upper bound of the confidence interval,
#' \item \code{std_error} the standard error of the point-estimate. }
#' @export monte_carlo_integrator
#' @examples
#' # Sophomore's dream: value to 10 decimals is 0.7834305107...
#' monte_carlo_integrator(function(x) x^x, 0, 1)
#' # Integral related to gamma function: exact value is Gamma(p+1), here p = 2
#' monte_carlo_integrator(function(x) log(1/x)^2, 0, 1)
monte_carlo_integrator <- function(g, a = 0, b = 1, naive = FALSE, n = 10^6, sig_lvl = 0.05, ...)
{
  if(!is.finite(a) || !is.finite(b))
  {
    stop("This integrator does not handle improper integrals over infinite regions. Endpoints of intervals must be finite.")
  }
  u <- stats::runif(n, min = a, max = b)
  if(!naive)
  {
    # Efficient estimator
    w <- 0.5*(g(u, ...)+g(b+a-u, ...))
  } else
  {
    # Naive estimator
    w <- g(u, ...)
  }
  x <- (b-a)*mean(w)
  # Sample mean is approximately normally distributed via CLT for large samples
  z_alpha <- stats::qnorm(1-sig_lvl/2)
  err <- stats::sd(w) * (b - a) / sqrt(n)
  lb <- x - z_alpha * err
  ub <- x + z_alpha * err
  # Gather results into data-frame
  results <- data.frame(estimate = x, lb = lb, ub = ub, std_error = err)
  return(results)
}

#' Multi-dimensional integration via Monte-Carlo methods
#'
#' @param g the known integrand, a function of a vector defined over a 'rectangular' region
#' @param lbs the left-end points of each interval per coordinate
#' @param ubs the right end-points of each interval per coordinate
#' @param n the number of variates to simulate of the IID uniform vector
#' @param sig_lvl the significance level for the confidence intervals
#' @param ... additional arguments to pass to the function \code{g}
#'
#' @description {A naive but efficient estimator for multi-dimensional integrals
#' via Monte-Carlo methods. The integrand must take a vector as an argument.}
#' @return data.frame containing
#' \itemize{
#' \item \code{estimate} the point estimate of the integral,
#' \item \code{lb} the lower bound of the confidence interval,
#' \item \code{ub} the upper bound of the confidence interval,
#' \item \code{std_error} the standard error of the point-estimate. }
#' @export multidim_integrator
#' @examples
#' # A simple example: integrating f(x,y,z)=xyz over [0,1]^3; exact value 1/8
#' multidim_integrator(function(x) x[1]*x[2]*x[3], lbs = c(0, 0, 0), ubs = c(1, 1, 1))
#'
#' # More involved example: approximating E(log(1+w dot R))
#' # where R is the vector of returns of N assets, w is the vector of weights
#' # Lower and upper bounds of returns, assuming Uniform distribution
#' lbs <- c(-1, -0.5, -0.25, 0)
#' ubs <- c(2, 0.5, 1, 1.5)
#' w <- c(0.1, 0.1, 0.1, 0.7)
#' N <- length(lbs) # Number of assets
#' # Mean's of returns from uniform model
#' mu <- 0.5*(ubs+lbs)
#' # Quadratic approximation requires the matrix E(R_iR_j) for all assets i,j
#' B <- matrix(0, N, N)
#' for(i in 1:N)
#' {
#'   for(j in 1:N)
#'   {
#'     B[i, j] <- mu[i]*mu[j]
#'   }
#' }
#' # Quadratic approximation to the expectation:
#' quad_approx <- as.numeric(t(mu)%*%w-0.5*t(w)%*%B%*%w)
#' # Function of a vector to pass to the integrator
#' g <- function(x, w)
#' {
#'   log(1+w%*%x)
#' }
#' # Since we want to approximate the expectation, we divided by the measure
#' # of the region we are integrating over.
#' multidim_integrator(g, lbs, ubs, w = w)/prod(ubs-lbs)
#' quad_approx
multidim_integrator <- function(g, lbs, ubs, n = 10^4, sig_lvl = 0.05, ...)
{
  # Dimensions come from intervals.
  N <- length(lbs)
  # Simulate n variates of vectors of IID uniform of dimension N
  u <- matrix(0, nrow = n, ncol = N)
  for(i in 1:n)
  {
    u[i, ] <- stats::runif(N, min = lbs, max = ubs)
  }
  # For each dimension and over every variate, transform u_i to b_i+a_i-u_i
  transformed_u <- t(apply(u, 1, function(x) ubs+lbs-x))
  # Compute g(u)
  gs <- apply(u, 1, function(x) g(x, ...))
  # Compute g(a+b-u)
  g_shifted <- apply(transformed_u, 1, function(x) g(x, ...))
  # Take the average of the two test statistics
  X <- 0.5*(gs+g_shifted)
  # Take the sample mean over the number of variates and multiply by measure of the region.
  estimate <- mean(X)*prod(ubs-lbs)
  # Sample mean is approximately normally distributed via CLT for large samples
  z_alpha <- stats::qnorm(1-sig_lvl/2)
  err <- stats::sd(X) * prod(ubs - lbs) / sqrt(n)
  lb <- estimate - z_alpha * err
  ub <- estimate + z_alpha * err
  # Gather results into data-frame
  results <- data.frame(estimate = estimate, lb = lb, ub = ub, std_error = err)
  return(results)

}

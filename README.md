
<!-- README.md is generated from README.Rmd. Please edit that file -->

# braggR

<!-- badges: start -->
<!-- badges: end -->

The goal of `braggR` is to provide easy access to the revealed
aggregator proposed in [Satopää
(2021)](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3769945).

## Installation

You can install the released version of `braggR` from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("braggR")
```

## Example

This section illustrates `braggR` on Scenario 2 in [Satopää
(2021)](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3769945).

``` r
library(braggR)
# Forecasters' probability predictions:
p = c(0.0625, 0.3125, 0.1250, 0.3125, 0.1875)

## Aggregate with a fixed base rate of 0.5.

# Sample the posterior distribution:
post_sample = sample_aggregator(p, p0 = 0.5, num_sample = 10^6, seed = 1)
# The posterior means of the model parameters:
colMeans(post_sample[,-1])
#>       rho     gamma     delta        p0 
#> 0.3821977 0.4742795 0.6561926 0.5000000
# The posterior mean of the level of rational disagreement:
mean(post_sample[,3]-post_sample[,2])
#> [1] 0.09208173
# The posterior mean of the level of irrational disagreement:
mean(post_sample[,4]-post_sample[,3])
#> [1] 0.1819131

# The revealed aggregator (a.k.a., the posterior mean of the oracle aggregator):
mean(post_sample[,1])
#> [1] 0.1405172
# The 95% credible interval of the oracle aggregator:
quantile(post_sample[,1], c(0.025, 0.975))
#>        2.5%       97.5% 
#> 0.001800206 0.284216903
```

This illustration aggregates the predictions in `p` by sampling the
posterior distribution `1,000,000` times. The base rate is fixed to
`p0 = 0.5`. By default, the level of burnin and thinning have been set
to `num_sample/2` and `1`, respectively. Therefore, in this case, out of
the `1,000,000` initially sampled values, the first `500,000` are
discarded for burnin. Given that thinning is equal to `1`, no more draws
are discarded. The final output `post_sample` then holds `500,000`draws
for the `aggregate` and the model parameters, `rho`, `gamma`, `delta`,
and `p0`. Given that `p0` was fixed to `0.5`, it is not sampled in this
case. Therefore all values in the final column of `post_sample` are
equal to `0.5.` The other quantities, however, show posterior
variability and can be summarized with the posterior mean. The first
column of `post_sample` represents the posterior sample of the oracle
aggregator. The average of these values is called the revealed
aggregator in [Satopää
(2021)](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3769945).
The final line shows the `95%` credible interval of the oracle
aggregator.

``` r
# Aggregate based on a prior beta(2,1) distribution on the base rate.
# Recall that Beta(1,1) corresponds to the uniform distribution.
# Beta(2,1) has mean alpha / (alpha + beta) = 2/3 and 
# variance alpha * beta / ((alpha+beta)^2*(alpha+beta+1)) = 1/18

# Sample the posterior distribution:
post_sample = sample_aggregator(p, alpha = 2, beta = 1, num_sample = 10^6, seed = 1)
# The posterior means of the oracle aggregator and the model parameters:
colMeans(post_sample)
#> aggregate       rho     gamma     delta        p0 
#> 0.1724935 0.5636953 0.6376554 0.9892552 0.6662238
```

This repeats the first illustration but, instead of fixing `p0` to
`0.5`, the base rate is now sampled from a `beta(2,1)` distribution. As
a result, the final column of `post_sample` shows posterior variability
and averages to a value close to the prior mean `2/3`.

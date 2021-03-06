---
title: "Scaling up abundance with confidence (intervals)"
author: "David J. Harris"
date: "August 14, 2015"
output: html_document
---

***Update**: This is all incorrect, since I misread Harte (2008).  The geometric distribution is only used when the total abundance of the species at the largest spatial scale is already known. I won't delete the repository in case it's useful to someone some day, but don't believe anything you read in it.


```r
set.seed(1) # set random number seed for reproducibility
```

## Parameter estimation

According to Harte et al. (2008), the distribution of abundance values, $n$, across plots of size $A$ follows a probability distribution that declines exponentially (assuming $A$ is less than half the total study area, which is the only case we typically care about). In statistical terminology, that means we want the geometric distribution, which describes the distribution of the number of failures in a sequence of Bernoulli trials before the first success occurs.  This distribution has a parameter, $p$ describing the probability of success in these trials.  Figure 1 shows two examples of probability distribution functions from this family, with $p=.8$  (in black) and $p=.5$ (in red).


```r
plot(0:10 - .05, dgeom(0:10, p = .8), type = "h", lwd = 3, xlab = "n", ylab = "frequency")
points(0:10 + .05, dgeom(0:10, p = .5), type = "h", lwd = 3, col = "red")
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2-1.png) 

Let's say we observe a site with $n=7$.  That outcome could have been generated from a range of possible geometric distributions, and we can assign a likelihood ($\mathrm{prob}(n|p)$) to them.

I'm pretty sure the likelihood is proportional to a beta distribution (because conjugate priors). Trial and error indicates that the distribution has $\mathrm{shape1} = 2$ and $\mathrm{shape2} = n + A$. While the maximum likelihood estimate for $p$ is $1/n$, the distribution is skewed so the posterior mean will be somewhat closer to $0.5$ (especially when $n$ is small).  Here, I rescale the units of $A$ so that it equals 1.


```r
Z = exp(4.276666) # Proportionality constant, found empirically

curve(dgeom(7, x), from = .0001, xlab = "p", ylab = "likelihood", lwd = 4, main = "scaled likelihood with one observed site")
curve(dbeta(x, shape1 = 2, shape2 = 7 + 1) / Z, add = TRUE, col = "red", lwd = 2, lty = 2)
legend("topright", c("dgeom", "dbeta"), lty = 1:2, lwd = c(4, 2), col = 1:2)
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png) 

If we'd observed $n$ at more than one site, then our estimate of $p$ would still follow a beta distribution, it would just have larger shape parameters (and less uncertainty). Trial and error seems to indicate that $\mathrm{shape1}=1+\sum_i{A_i}$ and $\mathrm{shape2} = 1 + \sum_i{n_i}$, where $i$ indexes the sites.  Here, I scaled $A$'s units so that the average area of the sampled sites equals one.


```r
Z = exp(8.291547) # Proportionality constant, found empirically

curve(dgeom(7, x) * dgeom(11, x) * Z, from = .0001, xlab = "p", ylab = "likelihood", lwd = 4, main = "scaled likelihood with two observed sites")
curve(dbeta(x, shape1 = 1 + 2, shape2 = 7 + 11 + 1), add = TRUE, col = "red", lwd = 2, lty = 2)
legend("topright", c("dgeom", "dbeta"), lty = 1:2, lwd = c(4, 2), col = 1:2)
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png) 

## Predicting one new site

This beta distribution describes our uncertainty about the "true" value of the parameter $p$.  I'm reasonably sure we can say $p=1/\bar{n}$, i.e. the reciprocal of the expected number of individuals in a plot of size $A$, but it would be worth double-checking.

If the parameter $p$ can take on different values from this distribution, then a wide range of outcomes are possible the next time we collect a sample from an area of size $A$.


```r
p_parameter_possibilities = rbeta(10000, 2, 7 + 2)

# 95% confidence interval for 1/p, which should be the expected abundance per plot 
#of size A
quantile(1 / p_parameter_possibilities, c(.025, .975))
```

```
##      2.5%     97.5% 
##  2.224526 41.316903
```

Based on this, our observation of 7 individuals in one area of size $A$ leads to a 95% confidence interval for expected abundance that ranges from 2.2 to 41.3.  This sounds reasonable: the geometric distribution has a long tail, so it's possible that the true mean is pretty small and 7 was an outlier. Likewise, the mode of a geometric distribution is always zero, so we'd expect to occasionally find a small number of individuals (like 7), even when the true mean is much higher.

The long tail also means that a newly-sampled site could have an absolutely enormous number of individuals, especially if the true mean density is much higher than we expected.


```r
sample_n_values = rgeom(10000, p_parameter_possibilities)
plot(table(sample_n_values))
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-1.png) 

This mixture distribution has a much longer tail than what we'd have gotten if we had assumed we knew $p$ exactly, because it accounts for the fact that $1/p$ could be huge. If we had just assumed $p = 1/n$ (the maximum likelihood estimate), we'd have expected a distribution like this instead:


```r
plot(table(rgeom(10000, 1/7)), xlim = range(sample_n_values))
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png) 

## Scaling up

Now let's scale up to multiple sites (e.g. 1000).  The METE model assumes that there's a single landscape-wide value of $p$, we just don't know what it is.  So for each possible value of $p$, we could sample 1000 sites of size $A=1$ and add up the $n$ values to get $n_0$.

With a fixed value of $p$, the sum of geometrically-distributed variables seems to be a negative binomial, with $\mathrm{size}=A_0/A$) and $\mathrm{prob}=p$ (which entails a mean/variance ratio equal to $p$).  I haven't really tested this, but it fits my understanding of the negative binomial versus the geometric (negative binomial generalizes the geometric distribution to describe time to *multiple* failures instead of just one) and seems to work correctly in the graph below.


```r
plot(density(replicate(100000, sum(rgeom(1000, p = 1/7)))), lwd = 4)
lines(0:10000, dnbinom(0:10000, size = 1000, prob = 1/7), col = "red", lwd = 2, lty = 2)
legend("topright", c("rgeom", "dnbinom"), lty = 1:2, lwd = c(4, 2), col = 1:2)
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-1.png) 


Again, if we average over our uncertainty, we get a much wider range of possibilities (note the range of the x axis).


```r
n_0_possibilities = sapply(
  rbeta(1000, 2, 7 + 2),
  function(p){
    sum(rgeom(1000, p))
  }
)

plot(density(n_0_possibilities))
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-1.png) 

While this is a uselessly wide distribution (95% confidence intervals between 1.2e+03 and 4.7e+04), it's also pretty sensible that we can't pin things down very well with only one quadrat.  With more samples, we could be more confident about the value of $p$, and we could cut the overdispersion down until things started to look like the previous figure when $A$ is large.

## Scaling up with spatial heterogeneity

If there's spatial heterogeneity, then we won't have the same $p$ at all sites. We'd thus expect an even broader range of possibilities when we scaled up than above (although we could keep the variance down by actually sampling multiple sites).  Depending on what kind of distribution we expect for $p$ across sites, we could do this in different ways. The most principled approach would be to fit a mixture of geometric distributions in the first step, instead of just one, and then carry the results through the rest of the process. Or we could estimate $p$ as a function of the environment (like fitting an SDM with abundance data and then taking the reciprocal of $n$ to get $p$ at all the unexplored areas).  Alternatively, we might be able to just increase the overdispersion of our negative binomial a bit further while, which should probably work fine in a lot of cases even if there isn't any particular reason to think the mixture distribution implied by spatial heterogeneity is still negative binomial.

## Remaining questions

* Did I make any boneheaded mistakes?
* What's the *joint* distribution of $n$ values for different species in a given area $A$, assuming $S_0$ is known/fixed and we've marginalized out $E$? This would make a useful prior for the kinds of multi-species distribution modeling I've been working on.
* I feel like there's a disconnect in my mind regarding when to expect geometric versus negative binomial distributions.  The math I'm using seems to indicate that one large quadrat would have a geometrically distributed $n$, but cutting it up into lots of tiny quadrats (small $A$) and summing them up would lead to negative binomially distributed $n$. I guess Harte (2008) says that the slope of the geometric distribution should change with the value of $A$, but it doesn't seem to allow for the same flexibility we get with the negative binomial.  So I'm still confused.

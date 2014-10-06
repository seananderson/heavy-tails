# Data selection

We followed the following data selection and quality control rules with the GPDD:

1. To remove populations with unreliable population indices that could be strongly confounded with economics and sampling effort, we removed all populations with a sampling protocol listed as `harvest` as well populations with the words `harvest` or `fur` in the cited reference title.

2. We removed all populations with uneven sampling intervals. I.e. we removed populations that didn't have a constant difference between the "decimal year begin" and “decimal year end” columns.

3. We removed all populations rated as $< 2$ in the GPDD quality assessment (on a scale of $1$ to $5$, with $1$ being the lowest quality data) [following @sibly2005; @ziebarth2010]

4. Populations with negative abundance values were assumed to be log values and transformed by taking $10$ to the power of the recorded abundance. In many cases this was noted for the population, but not in all cases. We inspected each of these \totalAssumedLog\ time series to make sure our assumption made sense.

5. We filled in all missing time steps with `NA` values and imputed single missing values with the geometric mean of the previous and following values. We chose a geometric mean to be linear on the log scale that the Gompertz and Ricker-logistic models were fit on.

6. We filled in single recorded values of zero with the lowest non-zero value in the time series [following @brook2006a]. This assumes that single values of zero result from abundance being low enough that censusing missed present individuals. We turned multiple zero values in a row into `NA` values. This implies that multiple zero values were either censusing errors or caused by emigration. Regardless, our population models were fit on a multiplicative (log) scale and so could not account for zero abundance.

7. We removed all populations with four or more identical values in a row since these suggest either recording error or extrapolation between two observations.

8. We removed all populations without at least four unique values [following @brook2006a].

9. We then wrote an algorithm to find the longest unbroken window of abundance (no `NA`s) with at least 20 time steps in each population time series. If there were any populations with multiple windows of identical length, we took the most recent window. This is a longer window than used in some previous analyses [e.g. @brook2006a], but since our model attempts to capture the shape of the distribution tails, our model requires more data.

10. Finally, we removed GPDD Main ID `20531`, which we noticed was a duplicate of `10139` (a heron population).

We provide a supplemental figure of all the time series included in our analysis and indicate which values were interpolated (\percImputedPops\% of populations had at least one point interpolated and only \percImputedPoints\% of the total observations were interpolated) (Fig. \ref{fig:all-ts}). Table S1 shows the final taxonomic breakdown and the number of populations with interpolated values.

# Details on the heavy-tailed Gompertz probability model

For the Gompertz model, our weakly-informative priors (Fig. \ref{fig:priors}) were:
$$\begin{aligned}
b &\sim \mathrm{Uniform}(-1, 2)\\
\lambda &\sim \mathrm{Normal}(0, 10^2)\\
\nu &\sim \mathrm{Truncated\mhyphen Exponential}(0.01, \mathrm{min.} = 2)\\
\sigma &\sim \mathrm{Half\mhyphen Cauchy} (0, 2.5).\end{aligned}$$

The $b$ parameter was bounded slightly past values that would imply a stationary time series so we can detect if the time series are non-stationary while being non-informative within the specified range. The estimates of $b$ were generally well within these bounds.

The prior on $\lambda$ is very weakly informative within the range of expected values for population growth and is similar to priors suggested by  \citet{gelman2008d} for intercepts of regression models.

The prior on $\nu$ is based on @fernandez1998; they chose a more informative exponential rate parameter of 0.1. 
We chose a less informative rate parameter of 0.01 and truncated the distribution at two, since at $\nu < 2$ the variance of the t distribution is undefined. This prior gives only a 7.7% probability $\nu < 10$ but constrains the sampling sufficiently to avoid wandering off towards infinity --- above approximately $\nu = 20$ or $\nu = 30$ the t distribution is so similar to the normal distribution (Fig. \ref{fig:didactic}) that time series of the length considered here are unlikely to be sufficiently informative about the precise value of $\nu$. In the scenario where the data are not heavy tailed (e.g. Fig. \ref{fig:didactic}e, h) the posterior will approximately match the prior (median $= 71$, mean $= 102$) and not be flagged as likely heavy tailed using the metrics we used in our paper (e.g. Pr($\nu < 10$) > 0.5).

The prior on $\sigma$ follows \citet{gelman2006c} and \citet{gelman2008d} and is based on our expected range of the value in nature from previous studies \citep[e.g.][]{connors2014}. In our testing, our parameter estimates were not substantially changed by switching to an uninformative uniform prior on $\sigma$, but the weakly informative prior substantially sped up chain convergence.

We fit our models with Stan 2.4.0 [@stan-manual2014], and R 3.1.1 [@r2014]. We began with four chains and 2000 iterations, discarding the first 1000 as warm up (i.e. 4000 total samples). If $\hat{R}$ (a measure of chain convergence) was greater than 1.05 for any parameter or the minimum effective sample size (a measure of the effective number of uncorrelated samples) for any parameter was less than 200, we doubled both the total iterations and warm up period and sampled from the model again. These thresholds are in excess of the minimums recommended by \citet{gelman2006a} of $\hat{R} < 1.1$ and effective sample size $> 100$ for reliable point estimates and confidence intervals. In the majority of cases our minimum thresholds were greatly exceeded. We continued this procedure up to 8000 iterations (16000 total samples) by which all chains were deemed to have sufficiently converged.

# Simulation testing the model

Throughout all of this — show that the model, if anything, under-predicts heavy tails but with sufficient data is unbiased.

2 parts: how many samples from the true population t distribution do you need to detect low nu? And, given that you have a set of deviations in which nu is detectable (effective nu is within 0.5 CV of true nu), can the more complex Gompertz still capture this?

First part: We drew from t distributions with different nu values and mean of 0, scale of 1. We started with 1600 samples and then fitted again at the first 800, 400, 200, 100, 50, 25. Each time we recorded the posterior samples of $\nu$.

Second part: To generate a series of process deviations with an effective nu approximately equal to the true nu, we generated process deviation sets repeatedly and estimated the mean, scale, and nu values each time. We recorded when the estimated nu was within 0.2 coefficient of variations of the true $\nu$ value and used this set of random seed values in our Gompertz simulation. We then fit our Gompertz models to the simulated datasets with all parameters (except $\nu$) set near the median values estimated in the GPDD.

# Alternative population models

## Ricker-logistic

We also fitted a Ricker-logistic model:
$$\begin{aligned}
x_t &= x_{t-1} + r_{\mathrm{max}}\left(1 - \frac{N_{t-1}}{K}\right) + \epsilon_t\\
\epsilon_t &\sim \mathrm{Student\mhyphen t}(\nu, 0, \sigma),
\end{aligned}$$
where $r_\mathrm{max}$ represents the maximum population growth rate that is obtained when $N$ (abundance) $= 0$. The parameter $K$ represents the carrying capacity and, as before, $x_t$ represents the $\ln$ transformed abundance at time $t$. The Ricker-logistic model assumes a linear decrease in population growth rate ($x_t - x_{t-1}$) with increases in abundance ($N_t$). In contrast, the Gompertz model assumes a linear decrease in population growth rate with increases in $\ln$ abundance ($x_t$) (REF).

To fit the Ricker-logistic models, we chose a prior on $K$ uniform between zero and twice the maximum observed abundance (@clark2010 chose uniform between zero and maximum observed, which is less conservative). We set the prior on $r_\mathrm{max}$ as uniform between 0 and 20 as in @clark2010. We used the same priors on $\nu$ and $\sigma$ as in the Gompertz model.

## Autocorrelated residuals

We considered a version of the Gompertz model in which an autoregressive parameter was fit to the process noise residuals:
$$\begin{aligned}
x_t &= \lambda + b x_{t-1} + \epsilon_t\\
\epsilon_t &\sim \mathrm{Student\mhyphen t}(\nu, \phi \epsilon_{t-1}, \sigma).
\end{aligned}$$

In addition to the parameters in the original Gompertz model, an additional parameter, $\phi$, is estimated that represents the relationship between of subsequent process noise residuals. Based on the results of previous analyses with the GPDD [e.g. @connors2014] and the chosen priors in previous analyses [e.g. @thorson2014a] and to greatly speed up chain convergence when running our model across all populations, we placed a weakly informative prior on $\phi$ that assumed the greatest probability density near zero with the reduced possibility of $\phi$ being near $-1$ or $1$. Specifically, we chose $\phi \sim \mathrm{Truncated\mhyphen Normal}(0, 1, \mathrm{min.} = -1, \mathrm{max.} = 1)$.

## Assumed observation error

We considered a version of the base Gompertz model that allowed for a specified level of observation error:
$$\begin{aligned}
U_t &= \lambda + b U_{t-1} + \epsilon_t\\
x_t &\sim \mathrm{Normal}(U_t, \sigma_\mathrm{obs})\\
\epsilon_t &\sim \mathrm{Student\mhyphen t}(\nu, 0, \sigma_\mathrm{proc}),
\end{aligned}$$

where $U$ represents the unobserved state vector, and $\sigma_\mathrm{obs}$ represents the standard deviation of observation error (on a log scale), which was set at $0.2$. TODO justify $0.2$.

# Modelling predictors of heavy-tailed dynamics

<!--Multilevel logistic regression model:-->
<!--$$\begin{aligned}-->
<!--\mathrm{Pr}(y_i = 1) &= -->
  <!--\mathrm{logit}^{-1}(\alpha + \alpha_{k[i]} + \alpha_{j[i]} + -->
    <!--X_i \beta + \epsilon_i)\\-->
<!--\alpha_j &\sim \mathrm{Normal}(0, \sigma^2_{\alpha[j]})\\-->
<!--\alpha_k &\sim \mathrm{Normal}(0, \sigma^2_{\alpha[k]})\\-->
<!--\epsilon &\sim \mathrm{Normal}(0, \sigma^2_{\epsilon}),-->
<!--\end{aligned}$$-->

We fit a multilevel beta regression model to the predicted probability of heavy tails, Pr($\nu < 10$), to investigate potential correlates of heavy-tailed dynamics.
The model was as follows:
$$\begin{aligned}
\mathrm{Pr}(\nu_i < 0.5) &\sim \mathrm{Beta}(A_i, B_i)\\
A_i &= \phi \mu_i\\
B_i &= \phi (1 - \mu_i)\\
\mu_i &= \mathrm{logit}^{-1}(\alpha 
+ \alpha_{\mathrm{class}[i]} 
+ \alpha_{\mathrm{order}[i]} 
+ \alpha_{\mathrm{species}[i]} 
+ X_i \beta)\\
\alpha_\mathrm{class} &\sim \mathrm{Normal}(0, \sigma^2_{\alpha[\mathrm{class}]})\\
\alpha_\mathrm{order} &\sim \mathrm{Normal}(0, \sigma^2_{\alpha[\mathrm{order}]})\\
\alpha_\mathrm{species} &\sim 
  \mathrm{Normal}(0, \sigma^2_{\alpha[\mathrm{species}]}),
\end{aligned}$$

\noindent
where $A$ and $B$ represent the beta distribution shape parameters, $\mu_i$ represents the predicted value for population $i$, $\phi$ represents the dispersion parameter, and $X_i$ represents a vector of predictors for population $i$ with associated $\beta$ parameters. The intercepts are allowed to vary from the overall intercept $\alpha$ by taxonomic classes ($\alpha_\mathrm{class}$), taxonomic orders ($\alpha_\mathrm{order}$), and species ($\alpha_\mathrm{species}$) with standard deviations $\sigma_{\alpha[\mathrm{class}]}$, $\sigma_{\alpha[\mathrm{order}]}$, and $\sigma_{\alpha[\mathrm{species}]}$. Where possible, we also allowed for error distributions around the predictors by incorporating the standard deviation of the posterior samples for the Gompertz parameters $\lambda$, $b$, and $\sigma$ around the mean point value. We log transformed $\sigma$, time-series length, and lifespan to match the way they are visually represented in Fig. \ref{fig:correlates}. All input variables were standardized by subtracting their mean and dividing by two standard deviations \citep{gelman2008c} to make their coefficients comparable in magnitude. 

We incorporated weakly informative priors into our model: $\mathrm{Cauchy}(0, 10)$ on the global intercept $\alpha$, $\mathrm{Half\mhyphen Cauchy}(0, 2.5)$ on all standard deviation parameters, $\mathrm{Half\mhyphen Cauchy}(0, 10)$ on the dispersion parameter $\phi$, and $\mathrm{Cauchy}(0, 2.5)$ on all other parameters \citep{gelman2006c,gelman2008d}. Compared to normal priors, the Cauchy priors concentrate more probability density at reasonably parameter values while allowing for higher probability density far into the tails, thereby allowing the data to dominate the posterior more strongly if it disagrees with the prior. Our conclusions were not qualitatively changed by using uninformative priors.

\renewcommand{\thetable}{S\arabic{table}}
\setcounter{table}{0}

\include{stat-table.tex}

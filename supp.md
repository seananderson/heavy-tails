# Data selection

To TODO …

1. To remove populations with unreliable population indices that could be strongly confounded with economics and sampling effort, we removed all populations with a sampling protocol listed as `harvest` as well populations with the words `harvest` or `fur` in the cited reference title.

2. We removed all populations with uneven sampling intervals. I.e. we removed populations that didn't have a constant difference between the "decimal year begin" and “decimal year end” columns.

3. We removed all populations rated as $< 2$ in the GPDD quality assessment (on a scale of $1$ to $5$, with $1$ being the lowest quality data) [following @sibly2005; @ziebarth2010]

4. Populations with negative abundance values were assumed to be log values and transformed by taking $10$ to the power of the recorded abundance. In many cases this was noted for the population, but not in all cases. We inspected each of these \totalAssumedLog\ time series to make sure our assumption made sense (Fig. \ref{fig:log10-assumed}).

5. We filled in all missing time steps with `NA` values and imputed single missing values with the geometric mean of the previous and following values. We chose a geometric mean to be linear on the log scale that the Gompertz and Ricker-logistic models were fit on.

6. We filled in single recorded values of zero with the lowest non-zero value in the time series [following @brook2006a]. This assumes that single values of zero result from abundance being low enough that censusing missed present individuals. We turned multiple zero values in a row into `NA` values. This implies that multiple zero values were either censusing errors or caused by emigration. Regardless, our population models were fit on a multiplicative (log) scale and so could not account for zero abundance.

7. We removed all populations with four or more identical values in a row since these suggest either recording error or extrapolation between two observations.

8. We removed all populations without at least four unique values [following @brook2006a].

9. We then wrote an algorithm to find the longest unbroken window of abundance (no `NA`s) with at least 20 time steps in each population time series. If there were any populations with multiple windows of identical length, we took the most recent window. This is a longer window than used in some previous analyses [e.g. @brook2006a], but since our model attempts to capture the shape of the distribution tails, our model requires more data.

10. Finally, we removed GPDD Main ID `20531`, which we noticed was a duplicate of `10139` (a heron population).

We provide a supplemental figure of all the time series included in our analysis and indicate which values were interpolated (\percImputedPops\% of populations had at least one point interpolated and only \percImputedPoints\% of the total observations were interpolated) (Fig. \ref{fig:all-ts}). Table S1 shows the final taxonomic breakdown and the number of populations with interpolated values.

# Details on the heavy-tailed Gompertz probability model

For the Gompertz, the mostly weakly-informative priors I'm using are:

$$\begin{aligned}
b &\sim \mathrm{Uniform}(-1, 2)\\
\lambda &\sim \mathrm{Normal}(0, 10^2)\\
\nu &\sim \mathrm{Exponential}(0.01)\\
\sigma &\sim \mathrm{Half\mhyphen Cauchy} (2.5).\end{aligned}$$

See Fig. \ref{fig:priors}.

Priors justification: (TODO this is old)

-  $b$ is bounded just past stationary so we can detect if they are non-stationary will keeping the sampler from wandering off too far

-  $\lambda$ prior (variance = 100) is basically uninformative within the range of expected values for population growth... it allows even a X probability at a value of X

-  $\nu$ is based on @fernandez1998; they chose a more informative 0.1 value, we chose a less informative 0.01 and justify it based on performance in supplemental figure X (gives X probability of value less than 10) but constrains the sampling somewhat since above 40 or 50 the shape of the t distribution is almost identical and data are not usually informative at these data quantities... sampler would head off to infinity otherwise

-  $\phi$ given a standard deviation of 1 given prior analyses which suggest autocorrelation in these datasets is minimal; prior is weak enough to allow high or low values, but given little information will stay near our expectation; similar approach in [@thorson2014a]

-  $\sigma_\mathrm{proc}$ can be justified based on @gelman2006c and the expected range of this variable in nature from previous studies [e.g. @connors2014]

We fit our models with Stan 2.4.0 [@stan-manual2014], and R 3.1.1 [@r2014]. I'm starting with 4 chains and 2000 iterations with the first 1000 as warmup (i.e. 4000 total samples). If rhat is greater than 1.05 for any parameter or the minimum effective sample size is less than 200 for any parameter then I double both the total iterations and warmup and run again. This continued up to 8000 iterations (16000 total samples) …

# Simulation testing the model

Throughout all of this — show that the model, if anything, under-predicts heavy tails but with sufficient data is unbiased.

2 parts: how many samples from the true population t distribution do you need to detect low nu? And, given that you have a set of deviations in which nu is detectable (effective nu is within 0.5 CV of true nu), can the more complex Gompertz still capture this?

First part: We drew from t distributions with different nu values and mean of 0, scale of 1. We started with 1600 samples and then fitted again at the first 800, 400, 200, 100, 50, 25. Each time we recorded the nu posterior.

Second part: To generate a series of process deviations with an effective nu approximately equal to the true nu, we generated process deviation sets repeatedly and estimated the mean, scale, and nu values each time. We recorded when the estimated nu was within 0.5 CVs of the true nu and used this set of random seed values in our Gompertz simulation. We then fit AR1 Gompertz models to the simulated datasets with all parameters (except nu) set near the median values estimated in the GPDD.

# Alternative population models

## Ricker-logistic

We also fitted a Ricker-logistic model:

$$\begin{aligned}
x_t &= x_{t-1} + r_{\mathrm{max}}\left(1 - \frac{N_{t-1}}{K}\right) + b x_{t-1} + \epsilon_t\\
\epsilon_t &\sim \mathrm{Student\mhyphen t}_\nu(0, \sigma),\end{aligned}$$

where $r_\mathrm{max}$ represents the maximum population growth rate that is obtained when $N$ (abundance) $= 0$. The parameter $K$ represents the carrying capacity and, as before, $x_t$ represents the $\ln$ transformed abundance at time $t$. The Ricker-logistic model assumes a linear decrease in population growth rate ($x_t - x_{t-1}$) with increases in abundance ($N_t$). In contrast, the Gompertz model assumes a linear decrease in population growth rate with increases in $\ln$ abundance ($x_t$) (REF).

To fit the Ricker-logistic models, we chose a prior on $K$ uniform between zero and twice the maximum observed abundance (@clark2010 chose uniform between zero and maximum observed, which is less conservative). We set the prior on $r_\mathrm{max}$ as uniform between 0 and 20 as in @clark2010. We used the same priors on $\nu$ and $\sigma$ as in the Gompertz model.

## Autocorrelated residuals

We considered a version of the Gompertz model in which an autoregressive parameter was fit to the process noise residuals:

$$\begin{aligned}
x_t &= \lambda + b x_{t-1} + \epsilon_t\\
\epsilon_t &\sim \mathrm{Student\mhyphen t}_\nu(\phi \epsilon_{t-1}, \sigma).\end{aligned}$$

In addition to the parameters in the original Gompertz model, an additional parameter, $\phi$, is estimated that represents the relationship between of subsequent process noise residuals. Based on the results of previous analyses with the GPDD [e.g. @connors2014] and the chosen priors in previous analyses [e.g. @thorson2014a] and to greatly speed up chain convergence when running our model across all populations, we placed a weakly informative prior on $\phi$ that assumed the greatest probability density near zero with the reduced possibility of $\phi$ being near $-1$ or $1$. Specifically, we chose $\phi \sim \mathrm{Truncated\mhyphen Normal}(0, 1, \mathrm{min.} = -1, \mathrm{max.} = 1)$.

## Assumed observation error

We considered a version of the base Gompertz model that allowed for a specified level of observation error:

$$\begin{aligned}
U_t &= \lambda + b U_{t-1} + \epsilon_t\\
x_t &\sim \mathrm{Normal}(U_t, \sigma_\mathrm{obs})\\
\epsilon_t &\sim \mathrm{Student\mhyphen t}_\nu(0, \sigma_\mathrm{proc}),\end{aligned}$$

where $U$ represents the unobserved state vector, and $\sigma_\mathrm{obs}$ represents the standard deviation of observation error (on a log scale), which was set at $0.2$. TODO justify $0.2$.

\renewcommand{\thetable}{S\arabic{table}}
\setcounter{table}{0}

\include{stat-table.tex}

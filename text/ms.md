**Title: Evidence for black-swan events in animal populations**

**Authors**: Sean C. Anderson^1,2\*^, Trevor A. Branch^3^, Andrew B. Cooper^2^, Nicholas K. Dulvy^1^

^1^Earth to Ocean Research Group, Department of Biological Sciences, Simon Fraser University, Burnaby BC, V5A 1S6, Canada

^2^School of Resource and Environmental Management, Simon Fraser University, Burnaby, BC, V5A 1S6, Canada

^3^School of Aquatic and Fishery Sciences, University of Washington, Box 355020, Seattle, WA 98195, USA

^\*^Corresponding author: Sean C. Anderson; E-mail: sean\_anderson@sfu.ca

**Abstract**:

Black swans are statistically improbable events that nonetheless occur---often with profound consequences[@taleb2007; @sornette2009]. While extremes in the physical environment, such as monsoons and heat waves, are widely studied and increasing in magnitude and frequency[@meehl2004; @katz2005; @ipcc2012], it remains unclear the extent to which ecological populations buffer or suffer from such extremes. We derived a model for estimating the frequency of black-swan events and applied it to abundance time series of 609 animal populations. We found strong evidence for black-swan events, but they were rare, occurring in 3–5% of populations: most frequently for birds (4–8%), mammals (4–6%), and insects (2–3%). Black-swan events were predominantly (\percBSDown %) downward, implying that unexpected population crashes occur far more frequently than unexpected population increases. Black swan events were not explained by any life-history covariates, but tended to be driven by external perturbations such as climate or severe winters. Data and simulations revealed that black swan-events are more likely to be detected in longer and less noisy time series. Accurate forecasts of extinction risk require accounting for the possibility of black swans and the predominance of crashes over increases, or predicted abundance will underestimate the probability of crashes by a factor of XX–XX. Our results demonstrate the critical importance of allowing for downward heavy-tailed events in population modelling and developing robust conservation and management strategies. Rare bad events can and do happen, and ignoring them will underestimate extinction risk.

\clearpage

Major surprises happen more often than expected in financial, social, and natural systems[@taleb2007; @sornette2009; @may2008]. Massive unpredictable market swings are responsible for the majority of financial gains and losses[@taleb2007], fatalities from the largest wars dwarf those from other wars[@newman2005], and the most damaging earthquakes greatly exceed their expected frequency[@sornette2009]. In ecology, background rates of global extinction are punctuated by mass extinction[@harnik2012], viruses can mutate suddenly to infect new hosts, and billions of animals can die at once in mass mortality events[@fey2015]. Indeed, catastrophes may be the most important element affecting population persistence[@mangel1994] and their importance is likely to increase given projected increases in the frequency and magnitude of climate-related extremes[@ipcc2012]. Despite this compelling anecdotal evidence for ecological black swans, systematic evidence in population time series of higher taxa has been elusive[@keitt1998; @allen2001; @halley2002] and the vast majority of model fitting and risk forecasting assumes that assumes that uncertainty can be represented by a normal distribution[@brook2006a; @dennis2006; @knape2012]. If black swans occur, though, a normal distribution would under-estimate the probability of extreme events occurring.  

Here, we assess the frequency and magnitude of black-swan dynamics across time series of 609 populations from a wide array of taxonomic groups---mostly birds, mammals, insects, and fishes (Extended Data Table 1). We then identify characteristics of time series or life-history traits associated with the detection of black-swan events and verify known causes. To accomplish this, we develop a framework for identifying heavy-tailed (black-swan) process noise in population dynamics, i.e. whether the largest stochastic jumps in log abundance from one time step to the next are more extreme than typically seen with a normal distribution. Our framework allows for a range of population dynamic models, can incorporate observation uncertainty and skewness in process noise, and can be readily applied to abundance time series.

Specifically, we fit population dynamics models in which the process noise is drawn from a Student-t distribution. By estimating the degrees of freedom parameter, $\nu$, we can estimate the degree to which the process deviations have heavy tails and are thereof evidence of black-swan events (Fig. 1, Extended Data Fig. \ref{fig:priors}). For example, at $\nu = 3$, the probability of drawing a value more than five standard deviations below the mean is $0.02$, whereas the probability of drawing such a value from a normal distribution is tiny ($2.9\cdot10^{-7}$). As $\nu$ approaches infinity, the distribution approaches the normal distribution (Fig. 1).

We found strong, but rare, evidence for black-swan dynamics: most frequently for birds (7%) followed by mammals (5%), and insects (3%) (Fig. 2). Black-swan events were taxonomically widespread, occurring in 38% of taxonomic orders. Accounting for time series length and partially pooling inference across taxonomic class and order with a hierarchical model, we found stronger evidence for black swans in insect populations than these statistics suggest---four of eight orders with the highest median probability of heavy tails were insect orders (Fig. 3a).

The majority of our heavy-tailed estimates were robust to alternative population models, observation error, and choice of Bayesian priors. Our conclusions were not systematically altered when we included an autocorrelation structure in the residuals, modelled population growth rates with or without density dependence, or modelled the population dynamics as Ricker-logistic instead of Gompertz (Extended Data Fig. \ref{fig:alt}). Introducing moderate observation error (CV = 0.2) slightly decreased the estimated prevalence of black-swan events (Extended Data Fig. \ref{fig:alt}). The strength of the prior on $\nu$ had little influence on estimates of black-swan dynamics (Extended Data Fig. \ref{fig:alt-priors}). Our simulation testing shows that, if anything, our models underpredict the true magnitude and probability of heavy tailed events---especially given the length of the time series in our analysis (Extended Data Figs \ref{fig:sim-nu}, \ref{fig:sim-prob}).

Across populations, the probability of observing black-swan dynamics was positively related to time-series length and negatively related to magnitude of process noise but not clearly related to population growth rate, density dependence, or maximum lifespan (Fig. 3b, Extended Data Fig. \ref{fig:correlates}). Longer time-series length was the strongest covariate of observing black-swan dynamics. For instance, the expected probability density suggesting black-swan dynamics was about 1.6 times greater for a population with 60 time steps compared to one with only 30. However, the absolute change in probability density with increased time series length was small (0.20 vs. 0.13 in this example, Extended Data Fig. \ref{fig:correlates}).

We examined all time series with published explanations of why black-swan events occurred (Extended Data Table 2). The majority (\percBSDown %) were downward events and involved a combination of multiple factors. For example, a synchronization of environmental- and predation-mediated population cycles is thought to have caused a downward black-swan event for a water vole (*Arvicola terrestris*) population[@saucy1994]. Other black swans were the result of a sequence of extreme climate events on their own. For instance, severe winters in 1929, 1940--1942, and 1962--1963 were associated with black-swan downswings in grey heron (*Ardea cinerea*) abundance in the United Kingdom[@stafford1971] (Fig. 4c). Our analysis finds that the last event was a combination of two heavy-tailed events in a row and the population took three times longer to recover than predicted[@stafford1971]. Downward black swans were sometimes followed by upward black swans. For example, during a period of population crowding and nest shortages, a population of European shag cormorants (*Phalacrocorax aristotelis*) on the Farne Islands, United Kingdom, declined suddenly following a red tide event in 1968[@potts1980]. This freed quality nest sites for first-time breeders, productivity increased, and the population experienced a rapid upswing in abundance[@potts1980].

Given the prevalence of downwards events, we refit our heavy-tailed models estimating the skewness of process noise ($\gamma$; Fig. 4a) and used these models to make near-term risk forecasts. Aggregated across populations with strong evidence of heavy tails (median $\nu < 10$), 85% of the $\gamma$ probability density was below 1, indicating strong evidence of downwardly skewed process noise (Fig. 4b, Extended Data Fig. \ref{fig:skew-nu}). In contrast, populations with little evidence of heavy tails (median $\nu \geq 70$) had little evidence of skewed process noise (Fig. 4b; 90% of 95% credible intervals overlapped 1). Projecting these heavy-tailed populations forward five years revealed that assuming the standard normal process noise underestimated risk by X--X fold (ratio of interquartile range of 99% lower credible interval of abundance; Fig. 4c, Extended Data Fig. TODO?).

<!-- or something about many ecological populations are better modelled allowing for heavy-tailed downwardly skewed events -->
Our analysis provides strong evidence for the presence of downwardly skewed black-swan events in ecological process noise of higher taxa. Unmodelled intrinsic properties of populations or extrinsic forces acting on populations could generate black-swan dynamics. Given that a mixture of normal distributions with different variances can generate a Student-t distribution with heavy tails[@allen2001], we could observe black-swan dynamics if we miss an underlying mixture of processes. Those processes might be an aggregation of populations across space, population diversity[@anderson2015], or an intrinsic change in population variability through time[@carpenter2006]. Extrinsic ecological forces could also cause black-swan dynamics[@nunez2012]. These forces could be extreme themselves: extreme climate[@meehl2004; @katz2005; @ipcc2012], predation from (or competition with) other species experiencing black swans, or sharp changes in human pressures such hunting, fishing, or habitat destruction might cause black-swan events. Alternatively, the synchrony of multiple extrinsic forces could give rise to black-swan dynamics. This could occur with synergistic interactions[@kirby2009] or even if non-synergistic forces experience a rare alignment[@denny2009].

<!--Recognizing the prevalence of heavy-tailed dynamics suggests a number of policy directions.-->
Ecological resource management can draw from disciplines that focus on heavy tails. For example, earthquake preparedness and response is focussed on black-swan events. Similarly to ecological black swans, we can rarely predict the specific timing of large earthquakes. But, earthquake preparedness involves spatial planning based on forecast probabilities to focus early detection efforts and develop disaster response plans[@nrc2007]. The presence of black swans also suggests that we develop management policy that is robust to heavy tails. For instance, setting target population abundances appropriately far back from critical limits may buffer black-swan events[@caddy1996], and maintaining genetic, phenotypic, and behavioural diversity may allow some components of populations to persist when others are affected by disease or extreme environmental forces[@hilborn2003; @schindler2010; @anderson2015]. Finally, extreme and unexpected, surprising, or counterintuitive ecological dynamics offer a tremendous opportunity to learn about ecological systems, evaluate when models break down, and adjust future management policy[@doak2008; @pine-iii2009; @lindenmayer2010].

Rare catastrophes can have a profound influence on population persistence[@mangel1994]. In recent decades, ecology has moved toward focussing on aspects of variance in addition to mean responses[@loreau2010a; @thompson2013]. Our results suggest that an added focus on ecological extremes represents the next frontier, particularly in the face of increased climate extremes[@meehl2004; @ipcc2012; @thompson2013]. Financial analysts are concerned with the shape of financial return downward tails because these directly impact estimates of risk---the probability of a specific magnitude of undesired event occurring[@rachev2008]. Similarly, ecologists should focus more on estimating and predicting downward tails of population abundance, since these increase true extinction risk.

**References**

\renewcommand{\section}[2]{}%
\bibliographystyle{naturemag}
\bibliography{jshort,ms}

**Online Content:** Methods, along with any Extended Data display items are available in the online version of the paper; references unique to these sections appear only in the online paper.

**Acknowledgements**: We thank J.W. Moore, A.O. Mooers, L.R. Gerber, J.D. Yeakel, C. Minto and members of the Earth to Ocean Group for helpful discussions and comments. We are grateful to the contributors and maintainers of the Global Population Dynamics Database and to Compute Canada's WestGrid high-performance computing resources. Silhouette images were obtained from phylopic.org under Creative Commons licenses; sources are listed in the Supporting Information. Funding was provided by a Simon Fraser University Graduate Fellowship (S.C.A.), the Natural Sciences and Engineering Research Council of Canada (N.K.D., A.B.C.), the Canada Research Chairs Program (N.K.D.).

**Author Contributions**: S.C.A. and T.A.B. conceived the project; S.C.A., T.A.B., A.B.C., and N.K.D. designed the study; S.C.A. analyzed the data and wrote the paper with input from all authors. 

**Author information**: The authors declare no competing financial interests. Reprints and permissions information is available at www.nature.com/reprints. Correspondence and requests for materials should be addressed to S.C.A. (sean_anderson@sfu.ca).

\clearpage

**Figures**

\begin{center}
\includegraphics[width=\textwidth]{../analysis/t-nu-eg2.pdf}
\end{center}

**Figure 1: Illustration of population dynamic models that allow for heavy tails. a, b,** Probability density for the Student-t distribution with scale parameter of 1 and different values of $\nu$. Small values of $\nu$ create heavy tails. As $\nu$ approaches infinity the distribution approaches the normal distribution. **c–e,** Simulated population dynamics from a Gompertz model with process noise drawn from Student-t distributions with three values of $\nu$. Coloured dots in panels **c** and **d** represent jumps with less than a 1 in 1000 chance of occurring in a normal distribution. **f–h,** Estimates of $\nu$ from models fit to the times series in **c–e**. Shown are posterior samples (histograms), median and interquartile range of the posterior (IQR, dots and line segments), and the exponential prior on $\nu$ (dashed lines). Colour shading behind **f–h** illustrates the approximate region of heavy tails.

\clearpage

\begin{center}
\includegraphics[width=0.36\textwidth]{../analysis/nu-coefs-2.pdf}
\end{center}

**Figure 2: Estimates of population dynamic heavy-tailedness for populations of birds, mammals, insects, and fishes.** Small values of $\nu$ suggest heavy-tailed black-swan dynamics. Vertical points and line segments represent posterior medians and 50% / 90% credible intervals for individual populations. Inset plots show probability that $\nu < 10$ for populations arranged by taxonomic order and sorted by decreasing mean Pr($\nu < 10$). Taxonomic orders with three or fewer populations in **a** are omitted for space. Red to yellow points highlight populations with a high to moderately high probability of heavy-tailed dynamics.

\clearpage

\begin{center}
\includegraphics[width=0.35\textwidth]{../analysis/order-posteriors-covariates.pdf}
\end{center}

**Figure 3: Standardized probabilities and covariates of heavy-tailed dynamics. a,** Taxonomic-order-level posterior densities of Pr($\nu < 10$) after accounting for time-series length in a hrarchical model. Estimates are at the geometric mean of time series length across all the data (approximately 27 time steps). Colour shading refers to taxonomic class (yellow: fishes, green: insects, purple: birds, and red: mammals). Dotted vertical line in **a** indicates the median expected Pr($\nu < 10$) from the prior distribution. **b,** Posterior densities for potential covariates of Pr($\nu < 10$). In both panels, short vertical line segments within the density polygons indicate medians.

\clearpage

\begin{center}
\includegraphics[width=0.65\textwidth]{../analysis/skew-fig.pdf}
\end{center}

**Figure 4: Evidence for downward skewness in process noise for heavy-tailed populations. a,** Illustrations of Student-t distributions with 3 levels of skewness ($\gamma$) and heavy-tails ($\nu$). **b,** Posterior density of the skewness parameters aggregated across populations grouped into heavy-tailed ($\hat{\nu} < 10$), slightly heavy-tailed ($10 \leq \hat{\nu} \leq 70$), and normal-tailed populations ($\hat{\nu} > 70$). Approximate mid-values from **b** are illustrated in **a**. **c,** Example time series of heavy-tailed populations with different levels of skewness. Red dots highlight likely heavy-tailed events. Forecasts (grey regions) show median (solid lines) and lower 99% credible intervals (dotted lines) of abundance. Black and red lines indicate forecasts from Gompertz models with lognormal error, and process noise with estimated tail shape and skewness, respectively.

\clearpage

**Methods**

We selected abundance time series from the Global Population Dynamics Database[@gpdd2010] (GPDD), which contains nearly 5000 time series of abundance from $\sim 1000$ species and $\sim 100$ taxonomic orders. We filtered the data (Supporting Information) to remove populations from less reliable data sources, and those without sufficient data for our models, and then interpolated some missing values (sensu Brook \textit{et al.}[@brook2006a]). Our interpolation affected only $\sim$ \interpPointsPerc \% of the final data points (Extended Data Table 1) and none of the data points that were later considered black-swan events. Our final dataset contained \NPops\ populations across \NOrders\ taxonomic orders and seven taxonomic classes, with a median of \medianTimeSteps\ time steps (range of \minTimeSteps--\maxTimeSteps) (Extended Data Table 1).

We fit heavy-tailed Gompertz population dynamics models to data from the GPDD. The Gompertz model represents population growth as a linear function in log space. If $x_{t}$ represents the log abundance ($N$) at time *t* then 

\begin{align}
x_t &= \lambda + bx_{t - 1} + \epsilon_t\\
\epsilon_t &\sim \mathrm{Student\mhyphen t}(\nu, 0, \sigma)
\end{align}

The growth parameter $\lambda$ represents the expected growth rate if $N_{t} = 1$. The model is density independent if $b = 1$, maximally density dependent if $b = 0$, and inversely density dependent if $b < 0$. Usually, the process noise $\epsilon_{t}$ is modelled as normally distributed, but we allow it to be drawn from a Student-t distribution with scale parameter $\sigma$ and degrees of freedom $\nu$. We can then estimate the degree to which the process deviations have heavy tails and are thereof evidence of black-swan events (Fig. 1a, b).

For the Gompertz model, our weakly-informative priors (Extended Data Fig. \ref{fig:priors}) were: 

\begin{align}
b &\sim \mathrm{Uniform}(-1, 2)\\ \lambda &\sim \mathrm{Normal}(0, 10^2)\\
\sigma &\sim \mathrm{Half\mhyphen Cauchy} (0, 2.5)\\ \nu &\sim
\mathrm{Truncated\mhyphen Exponential}(0.01, \mathrm{min.} = 2).
\end{align}

Our prior on $b$ was uninformative between values of $-1$ and $2$. We would not expect values of $b$ with levels of density dependence as low as $-1$ (very strong inverse density dependence), nor would we generally expect values above $1$. We allowed values of $b$ above $1$ to allow for non-stationary time series of growth rates. The estimates of $b$ were well within these bounds. Our prior on $\lambda$ was very weakly informative within the range of expected values for population growth and is similar to the default priors suggested by Gelman *et al.*[@gelman2008d] for intercepts of regression models. Our Half-Cauchy prior on $\sigma$ follows Gelman *et al.*[@gelman2006c] and Gelman *et al.*[@gelman2008d] and the specific scale parameter of $2.5$ is based on our expected range of the value in nature from previous studies[@connors2014]. In our testing of a subsample of populations, our parameter estimates were not qualitatively changed by switching to an uninformative uniform prior on $\sigma$, but the weakly informative Half-Cauchy prior sped up chain convergence.

We chose a weakly informative exponential prior on $\nu$ with rate parameter of $0.01$ truncated at values above two---a slightly less informative prior than suggested by Fernandez and Steel[@fernandez1998]. We truncated the distribution at $2$, since at $\nu < 2$ the variance of the t distribution is undefined. This prior gives only a $7.7$\% probability that $\nu < 10$ but constrains the sampling sufficiently to avoid drifting towards infinity. In any case, for $\nu > 20$ the t distribution is almost indistinguishable from the normal distribution (Fig. 1). Based on the shape of the t distribution, we chose the probability that $\nu < 10$, Pr($\nu < 10$), to define the probability of heavy-tailed dynamics. When categorizing a population as heavy or normal tailed, we applied a 0.5 probability threshold. In the scenario where the data are uninformative about heavy tails (e.g. Fig. 1e, h), the posterior will approximately match the prior (prior median $= 71$, mean $= 102$) and a metric of Pr$(\nu < 10) > 0.5$) is unlikely to flag the population as heavy tailed.

We fit our models with Stan 2.4.0[@stan-manual2014], and R 3.1.1[@r2014]. We began with four chains and $2000$ iterations, discarding the first $1000$ as warm up (i.e. 4000 final samples). If $\hat{R}$ (the potential scale reduction factor---a measure of chain convergence) was greater than $1.05$ for any parameter or the minimum effective sample size, $n_\mathrm{eff}$, (a measure of the effective number of uncorrelated samples) for any parameter was less than $200$, we doubled both the total iterations and warm up period and sampled from the model again. These thresholds are in excess of the minimums recommended by Gelman *et al.*[@gelman2006a] of $\hat{R} < 1.1$ and effective sample size $> 100$ for reliable point estimates and confidence intervals. In the majority of cases our minimum thresholds were greatly exceeded. We continued this procedure up to $8000$ iterations ($16000$ total samples) by which all chains were deemed to have sufficiently converged. The No-U-Turn Hamiltonian Markov chain Monte Carlo sampler in Stan generally requires far fewer iterations to obtain equivalent effective sample sizes than Gibbs or Metropolis–Hastings algorithms[@stan-manual2014].

There are a number of caveats when considering the generality of our results. The GPDD data represent a taxonomically and geographically biased sample of populations---the longer time series we focus on are dominated by commercially and recreationally important species and a disproportionate number of populations are located in the United Kingdom. Although we would expect to find qualitatively similar evidence for black swans in other large taxonomic or geographic samples of populations, the common forces driving those black swans would likely differ. Additionally, some apparent black-swan events could be recording mistakes, or conversely, some extreme observations may have been discarded or altered if they were erroneously suspected of being recording mistakes. Indeed, we discarded three of the populations that our method initially identified as heavy tailed because they turned out to be data-entry errors (Supporting Information). Finally, the temporal scales of observation and population dynamics vary considerably across populations in the GPDD and these likely influence the detection of heavy tails. For example, if we make frequent observations relative to generation time (e.g. for many large-bodied mammals) we will average across generations and perhaps miss black-swan events. Conversely, if we census populations infrequently relative to generation time (e.g. many insects in the GPDD) the recorded data may average across extreme and less-extreme events and also dampen black-swan dynamics.

**Alternative priors** To test if the prior on $\nu$ influenced our estimate of black-swan dynamics, we refit our models with weaker and stronger priors. Our base model used a prior on $\nu$ of Truncated-Exponential(0.02, min.\ = 2). For a weaker prior we used Truncated-Exponential(0.005, min.\ = 2) and for a stronger prior we used Truncated-Exponential(0.02, min.\ = 2) (Extended Data Fig. \ref{fig:priors}). Note that the base and weaker priors are relatively flat within the region of $\nu < 20$, which is the region we are mostly concerned about when categorizing populations as heavy- or thin-tailed. We furthermore considered a Uniform($0, 1$) prior on $1/\nu$ as used in Gelman *et al.*[@gelman2014], but this prior puts considerably more weight on lower (heavy-tailed) values of $\nu$ than the truncated exponential priors considered here (not shown).

Our results show that these weaker and stronger priors would have little influence on our conclusions about heavy-tailed dynamics (Extended Data Fig. \ref{fig:alt-priors}). When the data are informative about tail behaviour (i.e.  when there is strong evidence of low $\nu$ values, upper-right of Extended Data Fig. \ref{fig:alt-priors}), the prior has little impact on the estimate of $\nu$. When the data are less informative about $\nu$ (i.e.\ when there are no or few tail events and time series are short or noisy), the prior can pull the estimate of $\nu$ towards larger or smaller values (Extended Data Fig. \ref{fig:alt-priors}). The vast majority of the populations with Pr$(\nu < 10)$ in the base prior were not altered qualitatively by this range of prior strength.

**Alternative population models** We fit alternative population models to test if four key phenomena systematically changed our conclusions. The range of percentages of black swans by taxonomic class cited in the abstract are based on lower and upper limits across our main Gompertz model and these four alternative models.

*Autocorrelated residuals* We considered a version of the Gompertz model in which an autoregressive parameter was fit to the process noise residuals:

\begin{align}
x_t &= \lambda + b x_{t-1} + \epsilon_t\\
\epsilon_t &\sim \mathrm{Student\mhyphen t}(\nu, \phi \epsilon_{t-1}, \sigma).
\end{align}

In addition to the parameters in the original Gompertz model, this model estimates an additional parameter $\phi$, which represents the correlation of subsequent residuals. Based on the results of previous analyses with the GPDD[@connors2014] and the chosen priors in previous analyses[@thorson2014a] and to speed up chain convergence when running our model across all populations, we placed a weakly informative prior on $\phi$ that assumed the greatest probability density near zero with the reduced possibility of $\phi$ being near $-1$ or $1$. Specifically, we chose $\phi \sim \mathrm{Truncated\mhyphen Normal}(0, 1, \mathrm{min.} = -1, \mathrm{max.} = 1)$. The MCMC chains did not converge for \modelsNoConvergeAROne\ populations according to our criteria ($\widehat{R} < 1.05, n_\mathrm{eff} > 200$) after 8000 iterations of four chains. This included only \modelsNoConvergeAROneHeavyBase\ population in which Pr($\nu < 10$) $> 0.5$ categorized it as heavy in the main Gompertz model. We did not include this population in Extended Data Fig. \ref{fig:alt}.

*Assumed density independence* We fit a simplified version of the Gompertz model in which the density dependence parameter $b$ was fixed at $1$ (density independent). This is equivalent to fitting a random walk model (with drift) to the $\log$ abundances or assuming the growth rates are drawn from a stationary distribution. The model was as follows:

\begin{align}
x_t &= \lambda + x_{t-1} + \epsilon_t\\
\epsilon &\sim \mathrm{Student\mhyphen t}(\nu, 0, \sigma).
\end{align}

We fit this model for three reasons: (1) it is computationally simpler and so provides a check that our more complicated full Gompertz model was obtaining reasonable estimates of $\nu$, (2) it provides a test of whether density dependence was systematically affecting our perception of heavy tails, (3) it matches how some previous authors have modelled heavy tails without accounting for density dependence[@segura2013].

*Assumed observation error* Observation error can bias parameter estimates[@knape2012] and is known to affect the ability to detect extreme events[@ward2007]. In our main analysis, we fit a model that ignored observation error. One way to account for observation error would be to fit a full state-space model that simultaneously estimates the magnitude of process noise and observation error. However, simultaneously estimating observation and process noise is a challenging problem (e.g.\ because the observation and process noise parameters tend to negatively covary in model fitting) and is known to sometimes result in identifiability issues with the Gompertz population model[@knape2008]. Furthermore, our model was applied to hundreds of time series, often of short length (as few as \minTimeSteps\ time steps) and our model estimates an additional parameter---the shape of the process deviation tails---potentially making identifiability and computational issues even greater. Therefore, we considered a version of the base Gompertz model that allowed for a fixed level of observation error:

\begin{align}
U_t &= \lambda + b U_{t-1} + \epsilon_t\\
x_t &\sim \mathrm{Normal}(U_t, \sigma_\mathrm{obs}^2)\\
\epsilon_t &\sim \mathrm{Student\mhyphen t}(\nu, 0, \sigma_\mathrm{proc}),
\end{align}

where $U$ represents the unobserved state vector, $\sigma_\mathrm{obs}$ represents the standard deviation of observation error (on a log scale), and $\sigma_\mathrm{proc}$ represent the process noise scale parameter. We set $\sigma_\mathrm{obs}$ to $0.2$, which represents the upper limit of values often used in simulation analyses[@valpine2002; @thorson2014b].

*Ricker-logistic* We also fit a Ricker-logistic model:

\begin{align}
x_t &= x_{t-1} + r_{\mathrm{max}}\left(1 - \frac{N_{t-1}}{K}\right) + \epsilon_t\\
\epsilon_t &\sim \mathrm{Student\mhyphen t}(\nu, 0, \sigma),
\end{align}

where $r_\mathrm{max}$ represents the theoretical maximum population growth rate that is obtained when $N_t = 0$. The parameter $K$ represents the carrying capacity and, as before, $x_t$ represents the $\log$ transformed abundance at time $t$. The Ricker-logistic model assumes a linear decrease in population growth rate with increases in abundance. In contrast, the Gompertz model assumes a linear decrease in population growth rate with increases in \textit{log} abundance.

To fit the Ricker-logistic models, we chose a prior on $K$ uniform between zero and twice the maximum observed abundance (similar and less informative than in Clark *et al.*[@clark2010]). We set the prior on $r_\mathrm{max}$ as uniform between 0 and 20 as in Clark *et al.*[@clark2010]. We used the same priors on $\nu$ and $\sigma$ as in the Gompertz model.

**Simulation testing** We performed two types of simulation testing. First, we tested how easily the Student-t distribution $\nu$ parameter could be recovered given different true values of $\nu$ and different sample sizes. Second, we tested the ability of the heavy-tailed Gompertz model to obtain unbiased parameter estimates of $\nu$ given that a set of process deviations was provided in which the effective $\nu$ value was close to the true $\nu$ value.

We separated our simulation into these two components to avoid confounding two issues. (1) With smaller sample sizes, there may not be a stochastic draw from the tails of a distribution. In that case, no model, no matter how perfect, will be able to detect the shape of the tails. (2) Complex models may return biased parameter estimates if there are conceptual, computational, or coding errors. Our first simulation tested the first issue and our second simulation tested the latter. In general, our simulations show that, if anything, our model under predicts the magnitude and probability of heavy tailed events---especially given the length of the time series in the GPDD.

First, we tested the ability to estimate $\nu$ given different true values of $\nu$ and different sample sizes. We took stochastic draws from t distributions with different $\nu$ values ($\nu = 3, 5, 10,$ and $10^6$ [$\approx$ normal]), with central tendency parameters of $0$, and scale parameters of $1$. We started with $1600$ stochastic draws and then fit the models again at the first $800, 400, 200, 100, 50,$ and $25$ draws. Each time we recorded the posterior samples of $\nu$.

We found that we could consistently and precisely recover median posterior estimates of $\nu$ near the true value of $\nu$ with large samples ($\ge 200$) (Extended Data Fig. \ref{fig:sim-nu} upper panels). At smaller samples we could usually qualitatively distinguish heavy from not-heavy tails, but the model tended to underestimate how heavy the tails were. At the same time, at smaller sample sizes, the model tended to overestimate how large the scale parameter was (Extended Data Fig. \ref{fig:sim-nu} lower panels).

Secondly, we tested the ability of the heavy-tailed Gompertz model to obtain unbiased parameter estimates when the process noise was chosen so that appropriate tail events were present. To generate these process deviations for the $\nu = 3$ and $\nu = 5$ scenarios, we repeatedly drew proposed candidate process deviations and estimated the central tendency, scale, and $\nu$ values each time. We recorded when $\hat{\nu}$ (median of the posterior) was within $0.2$ CVs (coefficient of variations) of the true $\nu$ value and used this set of random seed values in our Gompertz simulation (Supporting Information).

We then fit our Gompertz models to the simulated datasets with all parameters (except $\nu$) set near the median values estimated in the GPDD. We repeated this with $50$ and $100$ samples without observation error, $50$ samples with observation error ($\sigma_\mathrm{obs} = 0.2$), and $50$ samples with the same observation error and a Gompertz model that allowed for correctly specified observation error magnitude. Our results indicate that the Gompertz model can recapture the true value of $\nu$ when the process noise was chosen so that appropriate tail events were present (Extended Data Fig. \ref{fig:sim-prob} upper panels). The addition of observation error caused the model to tend to underestimate the degree of heavy-tailedness. Fitting a model with correctly specified observation error did not make substantial improvements to model bias (Extended Data Fig. \ref{fig:sim-prob}).

When converting the posterior distributions of $\nu$ into Pr($\nu < 10$), the models distinguished heavy and not-heavy tails reasonably well (Extended Data Fig. \ref{fig:sim-prob} lower panels). Without observation error, and using a probability of $0.5$ as a threshold, the model correctly classified all simulated systems with normally distributed process noise as not heavy tailed. The model would have miscategorized only one of $40$ simulations at $\nu = 5$ across simulated populations with $50$ or $100$ time steps (Extended Data Fig. \ref{fig:sim-prob}, scenarios 1 and 2 in lower row, second panel from left). The model would have correctly categorized all cases where the process noise was not heavy tailed (Extended Data Fig. \ref{fig:sim-prob} bottom-right panel) and all cases where $\nu = 3$ and there was not observation error. With $0.2$ standard deviations of observation error, the model still categorized \obsErrorNuFivePerc\% of cases as heavy tailed when $\nu = 5$ and all but one case when $\nu = 3$. Allowing for observation error made little improvement to the detection of heavy tails. Therefore, we chose to focus on the simpler model without observation error in the main text, particularly given that the true magnitude of observation error was unknown in the empirical data.

**Covariates of heavy-tailed dynamics** We fit a hierarchical beta regression model to the predicted probability of heavy tails, Pr($\nu < 10$), to investigate potential covariates of heavy-tailed dynamics. We obtained maximum lifespan and body-size data from Brook *et al.*[@brook2006a]. The beta distribution is useful when response data range on a continuous scale between zero and one[@ferrari2004]. The model was as follows:

\begin{align}
\mathrm{Pr}(\nu_i < 10) &\sim \mathrm{Beta}(A_i, B_i)\\
\mu_i &= \mathrm{logit}^{-1}(\alpha
  + \alpha^\mathrm{class}_{j[i]}
  + \alpha^\mathrm{order}_{k[i]}
  + \alpha^\mathrm{species}_{l[i]}
  + X_i \beta),
  \: \text{for } i = 1, \dots, 617\\
A_i &= \phi_\mathrm{disp} \mu_i\\
B_i &= \phi_\mathrm{disp} (1 - \mu_i)\\
\alpha^\mathrm{class}_j &\sim
  \mathrm{Normal}(0, \sigma^2_{\alpha \; \mathrm{class}}),
  \: \text{for } j = 1, \dots, 6\\
\alpha^\mathrm{order}_k &\sim
  \mathrm{Normal}(0, \sigma^2_{\alpha \; \mathrm{order}}),
  \: \text{for } k = 1, \dots, 38\\
\alpha^\mathrm{species}_l &\sim
  \mathrm{Normal}(0, \sigma^2_{\alpha \; \mathrm{species}}),
  \: \text{for } l = 1, \dots, 301,
\end{align}

where $A$ and $B$ represent the beta distribution shape parameters; $\mu_i$ represents the predicted value for population $i$, class $j$, order $k$, and species $l$; $\phi_\mathrm{disp}$ represents the dispersion parameter; and $X_i$ represents a vector of predictors (such as lifespan) for population $i$ with associated $\beta$ coefficients. The intercepts are allowed to vary from the overall intercept $\alpha$ by taxonomic class ($\alpha^\mathrm{class}_j$), taxonomic order ($\alpha^\mathrm{order}_k$), and species ($\alpha^\mathrm{species}_l$) with standard deviations $\sigma_{\alpha \; \mathrm{class}}$, $\sigma_{\alpha \; \mathrm{order}}$, and $\sigma_{\alpha \; \mathrm{species}}$. Where possible, we also allowed for error distributions around the predictors by incorporating the standard deviation of the posterior samples for the Gompertz parameters $\lambda$, $b$, and $\log \sigma$ around the mean point value as normal distributions (not shown in the above equation).

We log transformed $\sigma$, time-series length, and lifespan to match the way they are visually represented in Extended Data Fig. \ref{fig:correlates} and to make the relationship approximately linear on the logit-transformed response scale. All input variables were standardized by subtracting their mean and dividing by two standard deviations to make their coefficients comparable in magnitude[@gelman2008c]. We excluded body length as a covariate because it was highly correlated with lifespan, and lifespan exhibited more overlap across taxonomy than body length. Lifespan is also more directly related to time and potential mechanisms driving black-swan dynamics.

We incorporated weakly informative priors into our model: $\mathrm{Cauchy}(0, 10)$ on the global intercept $\alpha$, $\mathrm{Half\mhyphen Cauchy}(0, 2.5)$ on all standard deviation parameters, $\mathrm{Half\mhyphen Cauchy}(0, 10)$ on the dispersion parameter $\phi_\mathrm{disp}$, and $\mathrm{Cauchy}(0, 2.5)$ on all other parameters[@gelman2006c; @gelman2008d]. Compared to normal priors, the Cauchy priors concentrate more probability density around expected parameter values while allowing for a higher probability density far into the tails, thereby allowing the data to dominate the posterior more strongly if it disagrees with the prior. We fit our models with 5000 total iterations per chain, 2500 warm-up iterations, four chains, and discarding every second sample. We checked for chain convergence visually and with the same criteria as before ($\widehat{R} < 1.05$ and $n_\mathrm{eff} >200$ for all parameters).

To derive taxonomic-order-level estimates of the probability of heavy tails accounting for time-series length (Fig. 3a), we fit a separate hierarchical model with the same structure but with only $\log$ time-series length as a predictor---in this case, we did not want to control for intrinsic population characteristics such as density dependence. Since our predictors were centered by subtracting their mean value, we obtained order-level estimates of the probability of heavy tails at mean log time-series length by adding the posteriors for $\alpha$, $\alpha^\mathrm{class}_j$, and $\alpha^\mathrm{order}_k$.

**Skewed Student-t forecasts** To formally evaluate the apparent skewness of heavy-tailed process noise we fit Gompertz models with skewed Student-t distributed process noise[@fernandez1998] (Supporting Information). The distribution adds one parameter to the Student-t distribution, $\gamma$, which controls the skewness. The distribution is symmetrical if $\gamma = 1$, left skewed if $0 < \gamma < 1$, and right skewed if $1 < \gamma < \infty$. We placed a weakly informative prior of $\mathrm{Cauchy}(0, 2.5)$ on $\log \gamma$. We aggregated 1000 randomly selected posterior samples from the $\gamma$ parameter of each model at three levels of evidence for heavy-tails: *heavy tailed* median $\nu < 10$; *slightly heavy tailed* $10 \le \mathrm{median}\, \nu \le 70$, and *normal tails*, median $\nu > 70$.

To generate 5-year forecasts of abundance, we combined the posterior parameter samples with stochastically generated process noise. We compared these forecasts to those generated by a standard Gompertz model with normally distributed process noise fit to the same data. We calculated the ratio of the lower 99% quantile credible interval between the two model projections. To ensure the results had stabilized in the tails of the forecast posterior, we increased the number of posterior samples in Stan. We ran four chains of 20,000 iterations and discarded the first 10,000 as warmup for a total of 30,000 samples.

<!-- TODO: for which pops did we run the skew ratio bit? -->

**Code availability.** Data and code to reproduce our analysis is available at <https://github.com/seananderson/heavy-tails>.

\clearpage

**Extended Data Tables**

\newenvironment{helvetica}{\fontfamily{phv}\selectfont}{\par}


\singlespacing

\textbf{Extended Data Table 1: Summary statistics for the filtered Global Population Dynamics Database time series arranged by taxonomic class.} Columns are: number of populations, number of taxonomic orders, numbers of species, median time series length, total number of interpolated time steps, total number of substituted zeros, and total number of time steps.

\onehalfspacing

\begin{helvetica}
\smallskip
\begin{scriptsize}
\begin{tabular}{lrrrrrrrr}
\toprule
\input{../analysis/stat-table.tex}
\label{tab:stats}
\end{tabular}
\end{scriptsize}
\end{helvetica}

\clearpage

\LTcapwidth=\textwidth
<!--\bibpunct{}{}{;}{a}{}{;}-->

\renewcommand{\tablename}{\textbf{Extended Data Table}}
\setcounter{table}{1}

\renewcommand{\arraystretch}{0.1}% Tighter

\begin{helvetica}
\singlespacing
\begin{scriptsize}
\begin{longtable}{>{\RaggedRight}m{1.4cm}>{\RaggedRight}p{6.2cm}>{\RaggedRight}p{0.7cm}>{\RaggedRight}p{0.4cm}>{\RaggedRight}p{4.4cm}>{\RaggedRight}p{1.4cm}}
\caption{\textbf{Populations with a high probability of heavy-tailed dynamics in the base heavy-tailed Gompertz population dynamics model.} Shown are the log abundance time series, population descriptions, Global Population Dynamics Database Main IDs, citation for the data source or separate verification literature, a description of the cause of the black swan events (if known), the probability of heavy tails as calculated by our model, and median estimate of $\nu$ from our model with 90\% quantile credible intervals indicated in parentheses. Red dots on the time series indicate downward black-swan events and blue values indicate upward black-swan events that have a $10^{-4}$ probability or less of occurring if the population dynamics were explained by a Gompertz model with normally distributed process noise with a standard deviation equal to the scale parameter in the fitted t distribution.}\\

\toprule
\input{../analysis/cause-table.tex}
\label{tab:causes-supp}
\end{longtable}
\end{scriptsize}
\onehalfspacing
\end{helvetica}


\clearpage

**Extended Data Figures**

\renewcommand{\figurename}{\textbf{Extended Data Fig.}}
<!--\def\fnum@figure{\figurename\nobreakspace\textbf{\thefigure}}-->

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=0.6\textwidth]{../analysis/priors-gomp-base.pdf}

\caption{\textbf{Probability density of the Bayesian priors for the Gompertz models.} \textbf{a,} per capita growth rate at $\log$(abundance) = $0$: $\lambda \sim \mathrm{Normal}(0, 10^2)$. \textbf{b,} scale parameter of t-distribution process noise: $\sigma \sim \mathrm{Half\mhyphen Cauchy} (0, 2.5)$ \textbf{c,} Student-t distribution degrees of freedom parameter: $\nu \sim \mathrm{Truncated\mhyphen Exponential}(0.01, \mathrm{min.} = 2)$. \textbf{d,} AR1 correlation coefficient of residuals: $\phi \sim \mathrm{Truncated \mhyphen Normal}(0, 1, \mathrm{min.} = -1, \mathrm{max.} = 1)$. Not shown is $b$, the density dependence parameter: $b \sim \mathrm{Uniform}(-1, 2)$. Panel \textbf{c} also shows two alternative priors: a weaker prior $\nu \sim \mathrm{Truncated\mhyphen Exponential}(0.005, \mathrm{min.} = 2)$, and a stronger prior $\nu \sim \mathrm{Truncated\mhyphen Exponential}(0.02, \mathrm{min.} = 2)$. The inset panel shows the same data but with a log-transformed x axis. Note that the base and weaker priors are relatively flat within the region of $\nu < 20$ that we are concerned with.}

\label{fig:priors}
\end{center}
\end{figure}

\clearpage

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=0.8\textwidth]{../analysis/gomp-comparison.pdf}

\caption{\textbf{Estimates of $\nu$ from alternative models plotted against the base Gompertz model estimates of $\nu$.} Shown are medians of the posterior (dots) and 50\% credible intervals (segments). The diagonal line indicates a one-to-one relationship. Different colours indicate various taxonomic classes. The grey-shaded regions indicate regions of disagreement if $\nu = 10$ is taken as a threshold of heavy-tailed dynamics. The Gompertz observation error model assumes a fixed standard deviation of observation error of $0.2$ on a log scale.}

\label{fig:alt}
\end{center}
\end{figure}

\clearpage

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=0.8\textwidth]{../analysis/gomp-prior-comparison.pdf}

\caption{\textbf{Estimates of $\nu$ from Gompertz models with alternative priors on $\nu$.} Shown are medians of the posterior (dots) and 50\% credible intervals (segments). The diagonal line indicates a one-to-one relationship. Different colours indicate various taxonomic classes. The grey-shaded regions indicate regions of disagreement if $\nu = 10$ is taken as a threshold of heavy-tailed dynamics. The base, weaker, and stronger priors on $\nu$ are illustrated in Fig. \ref{fig:priors}. In general, the estimates are nearly identical in cases where the data are informative about low values of $\nu$. When the data are less informative about low values of $\nu$, the prior can slightly pull the estimates of $\nu$ towards higher or lower values.}

\label{fig:alt-priors}
\end{center}
\end{figure}

\clearpage


\begin{figure}[htbp]
\begin{center}
\includegraphics[width=0.8\textwidth]{../analysis/t-dist-sampling-sim-prior-exp0point01.pdf}
\includegraphics[width=0.8\textwidth]{../analysis/t-dist-sampling-sim-sigma-prior-exp0point01.pdf}

\caption{\textbf{Testing the ability to estimate $\nu$ (top panels) and the scale parameter of the process deviations (bottom panels) for a given number of samples (columns) drawn from a distribution with a given true $\nu$ value (rows).} The red lines indicate the true population value. When a small number of samples are drawn there may not be samples sufficiently far into the tails to recapture the true $\nu$ value; however, heavy tails are still distinguished from normal tails in most cases, even with only 25 or 50 samples.}

\label{fig:sim-nu}
\end{center}
\end{figure}

\clearpage

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=\textwidth]{../analysis/sim-gompertz-median-dist.pdf}
\includegraphics[width=\textwidth]{../analysis/sim-gompertz-p10.pdf}

\caption{\textbf{Simulation testing the Gompertz estimation model when the process deviation draws were chosen so that $\nu$ could be estimated close to the true value outside the full population model} (``effective $\nu$'' within a CV of 0.2 of specified $\nu$). Upper panels show the distribution of median $\widehat{\nu}$ across 20 simulation runs. Lower panels show the distribution of Pr($\nu < 10$) across 20 simulation runs. We ran the simulations across three population (``true'') $\nu$ values (3, 5, and $1\cdot 10^9$, i.e.\ approximately normal) and four scenarios: (1) 100 time steps and no observation error, (2) 50 time steps and no observation error, (3) 50 time steps and observation error drawn from $\mathrm{Normal} (0, 0.2^2)$ but ignored, and (4) 50 time steps with observation error in which the quantity of observation error was assumed known. Within each scenario the dots represent stochastic draws from the true population distributions combined with model fits. Underlayed boxplots show the median, interquartile range, and $1.5$ times the interquartile range.}

\label{fig:sim-prob}
\end{center}
\end{figure}

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=0.8\textwidth]{../analysis/correlates-p10.pdf}

\caption{\textbf{Potential covariates of heavy-tailed population dynamics (indicated by a high probability that $\nu < 10$).} Shown are \textbf{a--c} parameters from the Gompertz heavy-tailed population model ($\sigma$, $\lambda$, $b$), \textbf{d} number of time steps, \textbf{e} body length, and \textbf{f} lifespan. For the Gompertz parameters, $\sigma$ refers to the scale parameter of the Student-t process-noise distribution, $\lambda$ refers to the expected log abundance at the next time step at an abundance of one, $b$ refers to the density dependence parameter. Circles representing a few sharks, crustaceans, and gastropods are filled in white. Median and 90\% credible interval posterior predictions of a beta regression hierarchical model are shown in panels \textbf{a} and \textbf{d} where there was a high probability the slope coefficient was different from zero (Fig. 3b).}

\label{fig:correlates}
\end{center}
\end{figure}

\clearpage

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=6in]{../analysis/skew-t-illustration}
\includegraphics[width=4in]{../analysis/skewness-vs-nu}

\caption{\textbf{something} something TODO and/or have the time series projections of a series of heavy and semi-heavy tailed time series.}

\label{fig:skew-nu}
\end{center}
\end{figure}

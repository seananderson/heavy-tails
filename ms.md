*Authors are listed alphabetically (ABCD!) and open to reordering.*

^1^Earth to ocean research group, department of biological sciences, Simon Fraser University, Burnaby BC, V5A 1S6, Canada

^2^School of Aquatic and Fishery Sciences, University of Washington, Box 355020, Seattle, WA 98195, USA

^3^School of Resource and Environmental Management, Simon Fraser University, Burnaby, BC, V5A 1S6, Canada

^\*^Corresponding author: Sean C. Anderson; Earth to Ocean Research Group, Department of Biological Sciences, Simon Fraser University, Burnaby BC, V5A 1S6; Phone: 1-778-782-3989; E-mail: sean_anderson@sfu.ca

\clearpage

*Main message, to be deleted (25 words)*: Black swans are present & taxonomically widespread but rare in population dynamics. Extreme climate, predation, parasites, and their interactions are common causes; intrinsic drivers are unclear.

# Abstract

Black swans are statistically improbable events that nonetheless occur --- often with profound implications. While extremes in the physical environment, such as monsoons and heat waves, are widely studied and increasing in magnitude and frequency, it remains unclear the extent to which ecological populations buffer or suffer from extremes. Here, we develop a probability model to estimate the degree of heavy-tailedness (presence of black swans) in ecological process noise. We apply our model to \NPops\ time series from around the world across \NOrders\ taxonomic orders and seven classes. We find strong evidence of black swans, but they are rare, occurring in \overallMinPerc--\overallMaxPerc\% of populations; most frequently for birds (\AvesRangePerc\%) followed by mammals (\MammaliaRangePerc\%), insects (\InsectaRangePerc\%), and fishes (\OsteichthyesRangePerc\%). When they occur, they tend to be driven by climate and severe winters, cycles of parasites and predators, and interactions between these elements. We find little evidence of intrinsic life-history correlates. *close with modelling and policy implications*

# Introduction

1. Black swans are... One of the most striking black swans in ecology was the asteroid marking the mass extinctions at the K-T boundary. Today, it is the extremeness of climate that is expected to cause the greatest societal damage. While extremes in the physical environment are present and widely understood, it remains unclear the extent to which ecological systems buffer or suffer from black swans. (re-write that)

2. What is the current evidence for heavy tails in ecological time series?

3. There are two main possible reasons for why we may find little evidence of ecological black swans. (1) don’t exist or (2) a paucity of appropriate methods or length/quality of time series.

4. Here we develop... our research questions:
How frequent and strong are black swans in population dynamics across taxonomy? Are there characteristics of time series or intrinsic life-history characteristics associated with heavy-tailed dynamics? What are some verified causes of black swans in ecological timeseries?
Main conclusions...

Papers to work in:

[@inchausti2002; @halley2002; @inchausti2001]
[@jentsch2007]
[@ward2007]
[@garcia-carreras2011] [@sornette2009]
[@nunez2012]
[@thompson2013] [@beaugrand2012] [@pine-iii2009]
[@doak2008]
Extreme rain and monsoons over Inda increasing [@goswami2006]
[@smale2013]
[@easterling2000] [@scheffer2003] [@katz2005]
[@taleb2007]
[@vasseur2014]
[@vert-pre2013] [@lindenmayer2010]
[@valpine2002] [@gregory2010] @garcia-carreras2011 @brook2006 [@herrandoprez2014]
[@sibly2005; @ziebarth2010]

# Methods

To obtain estimates of the probability and magnitude of population dynamic black swans, we fit multiple commonly applied phenomenological population dynamics models to abundance time series from around the world. Compared to previous analyses of heavy tails, (1) our analysis is broad in taxonomic scope, (2) allows for population dynamics rather than  modelling abundance or growth rates, and (3) is based on a model that can estimate the probability of heavy tails directly, rather than competing multiple models. Compared to previous population dynamics analyses, the key addition of our models was estimating the shape of the tails of the process noise --- how extreme the largest jumps from time step to time step were beyond what a population model allows for. We estimate how heavy tailed process noise is across populations and then investigate how these estimates relate to taxonomy, time series properties, and life-history characteristics.

We obtained abundance time series from the Global Population Dynamics Database (GPDD) (REF). To derive a high-quality subset of populations that were suitable for our analyses, we implemented a number of filtering rules (Supplementary Material). Our filtering procedure removed populations from less reliable data sources, removed those without sufficient data, and interpolated some missing values (Supplementary Material). Our final dataset contained \NPops\ populations across \NOrders\ taxonomic orders and seven taxonomic classes (Table S1, Fig. \ref{fig:all-ts}).

Our main analysis focusses on the Gompertz population model. The Gompertz model is a commonly applied phenomenological population dynamics model (e.g. REF). In previous analyses of the GPDD, it has proven to be most frequently parsimonious model (REF). We extend the traditional Gompertz model by allowing the process noise to be drawn from a Student-t distribution with degrees of freedom parameter $\nu$. If $\nu$ is small ($\lesssim 10$) the t distribution has much heavier tails than a normal distribution. For example, at $\nu = 2$, the probability of drawing a value more than five standard deviations from the mean is 1.8%, whereas the probability of drawing such a value from a normal distribution is nearly zero ($2.9\cdot10^{-5}$%). As $\nu$ approaches infinity the distribution approaches the normal distribution (Fig. \ref{fig:didactic}). By estimating the value of $\nu$, constrained by a weakly informative exponential prior, we can quantity how heavy-tailed the process noise deviations are. Specifically, we chose an exponential prior with rate parameter of $0.01$ (Fig. \ref{fig:priors}), a slightly less informative prior than that suggested by @fernandez1998.

The heavy-tailed Gompertz model is:
\begin{align*}
x_t &= \ln N_t\\
x_t &= \lambda + b x_{t-1} + \epsilon_t\\
\epsilon_t &\sim \mathrm{Student\mhyphen t}(\nu, 0, \sigma),
\end{align*}
where $x_t$ is the $\ln$ abundance $N$ at time $t$. The model is density independent if $b = 1$, maximally density dependent if $b = 0$, and inversely density dependence if $b < 0$. The parameter $\lambda$ represents the expected population growth rate at $x_t = 0$. The process noise is modelled as a Student-t distribution with degrees of freedom parameter of $\nu$, mean of  $0$, and a scale parameter of $\sigma$. We fit the Gompertz models in a Bayesian framework using Stan [@stan-manual2014] assuring that chains had sufficiently converged and the sampler had obtained sufficient independent samples from the posterior (Supporting Material).

We fit three alternative population models to test if they systematically changed our conclusions (Supporting Material). (1) Autocorrelation has been suggested as a reason for increased observed variability of abundance time series through time, which could create apparent heavy tails [@inchausti2002]; therefore, we fit a model that modelled the serial correlation in the residuals. (2) Observation error could conceivably create false heavy tails or mask our ability to detect heavy tails. However, simultaneously estimating observation error and process noise in a state space framework is a challenging computational problem, particularly in our case where many time series were relatively short and we were attempting to estimate the shape of the process noise tails. Therefore, we fit a model in which we assumed a level of observation error ($0.2$ standard deviation on a log scale). (3) The Gompertz model assumes that per capita growth rate declines linearly with log abundance. We also fit an alternative model, the Ricker-logistic model, which assumes that per capita growth rate declines linearly with abundance itself (Supporting Material).

We investigated possible correlates of heavy-tailed population dynamics visually and through multilevel modelling. We plotted characteristics of the time series ($\sigma$, $\lambda$, $b$, and time-series length) along with two life-history characteristics (body length and maximum lifespan obtained from @brook2006a) against our estimates of $\nu$. We formally investigated these relationships by fitting multilevel models to two response variables: (1) the continuous response of the probability that $\nu$ was less than $10$ (i.e. approximately the probability of heavy tails), and (2) the simpler dichotomous response variable p($\nu < 10$) $> 0.5$ (heavy tailed) vs.\ p($\nu < 10$) $\le 0.5$ (not heavy tailed). To account for broad patterns of phylogenetic relatedness, we allowed for varying intercepts at the taxonomic class level and nested varying intercepts at the taxonomic order level.

Finally, we investigated a sample of populations that our method categorized as having a high probability of heavy tails. Where possible, we found the primary literature the dataset was derived from. As necessary we found related literature to verify the causes. These populations were chosen in a haphazard manner with the purpose of hypothesis generating about the causes of population dynamic black swans.

---

Supporting material to be noted:

1. Table showing the taxonomic breakdown, number of interpolated points, etc (Table S1).

2. Plot showing time series we assumed were recorded as log transformed (Fig. \ref{fig:log10-assumed}). (remove this?)

3. Plot showing all the time series we used with interpolated and zero imputed values (Fig. \ref{fig:all-ts}).

4. Plot showing the priors (Fig. \ref{fig:priors}).

# Results

Main figures and tables:

1. Fig. \ref{fig:didactic} illustrates the method: it shows the t-distribution tails, example simulated time series, and model fits to those time series.

2. Fig. \ref{fig:nu-coefs} shows the posterior distributions of the estimated $\nu$ values. These are split by taxonomic class and order.

3. Fig. \ref{fig:correlates} shows possible biological and time-series-property correlates of heavy-tailed behaviour. Modelling output in supplement.

4. Table. \ref{tab:sparks} selects populations that were categorized as heavy-tailed and digs into the causes. These are a non-random sample that I was able to verify in the literature.

Supporting material to be noted quickly:

1. Simulation testing: the ability to recover $\nu$ when randomly sampling from various distributions (Fig. \ref{fig:sim-nu}); heavy-tailed Gompertz model performance and confidence interval coverage given the process noise has deviations that are known to have effective $\nu$ equal to true $\nu$ (Fig. \ref{fig:sim-gompertz}); boxplots of the same output (Fig. \ref{fig:sim-gompertz-boxplots}); probability that $\nu < 10$ for the Gompertz simulation testing (Fig. \ref{fig:sim-prob}).

2. Time series of all heavy tailed populations (Fig. \ref{fig:heavy-ts}).

3. The effect of alternative population dynamics models on $\nu$ estimates (Fig. \ref{fig:alt}).

4. Coefficients from modelling covariates of the probability of heavy tails (Fig. \ref{fig:correlate-coefs}).

5. Example Stan code for heavy-tailed Gompertz model with AR1 residuals and observation error (a model including everything).

6. GPDD IDs used.

<!--Modelling result: an increase of 1 time step of data (given that you are ballpark around the mean — 30 time steps — to start with) results in an approximately 1% increase in the expected probability of observing heavy tails. (Using Gelman’s ‘divide by 4’ rule for interpreting logistic regression coefficients.)-->

# Discussion

1. Summary of our findings

2. How do our results mesh with previous related analyses?

3. Why might we expect to see heavy tails in ecological time-series? (mixture
   of normals, red-coloured noise, tipping points, interactions, extreme
   drivers)

4. Are the observed frequencies by taxonomic class real or an observational
   phenomena?

5. Ways forward

6. Policy implications

# Acknowledgements

Funding: SCA: SFU Graduate Fellowship, NKD: NSERC SFU ..., ABC: ..., TAB: ...

Earth to Ocean research group for helpful discussions

Global Population Dynamics Database

Silhouettes: `phylopic.org`: rabbit by Sarah Werning,
Aves photo by Jean-Raphaël Guillaumin (photography) and T. Michael Keesey (vectorization)
, hoverfly by Gareth Monger. All under Creative Commons Attribution 3.0 Unported license.

Compute Canada’s WestGrid high-performance computing resources

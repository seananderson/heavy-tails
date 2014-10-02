*Authors are listed alphabetically (ABCD!) and open to reordering.*

^1^Earth to ocean research group, department of biological sciences, Simon Fraser University, Burnaby BC, V5A 1S6, Canada

^2^School of Aquatic and Fishery Sciences, University of Washington, Box 355020, Seattle, WA 98195, USA

^3^School of Resource and Environmental Management, Simon Fraser University, Burnaby, BC, V5A 1S6, Canada

^\*^Corresponding author: Sean C. Anderson; Earth to Ocean Research Group, Department of Biological Sciences, Simon Fraser University, Burnaby BC, V5A 1S6; Phone: 1-778-782-3989; E-mail: sean_anderson@sfu.ca

\clearpage

*Main message, to be deleted (25 words)*: Black swans are present & taxonomically widespread but rare in population dynamics. Extreme climate, predation, parasites, and their interactions are common causes; intrinsic drivers are unclear.

# Abstract

Black swans are statistically improbable events that nonetheless occur --- often with profound implications. While extremes in the physical environment, such as monsoons and heat waves, are widely studied and expected to increase over the next century, it remains unclear the extent to which ecological populations buffer or suffer from extremes. Here, we develop a probability model to estimate the degree of heavy-tailedness (presence of black swans) in ecological process noise. We apply our model to \NPops\ time series from around the world across \NOrders\ taxonomic orders and seven classes. We find strong evidence of black swans, but they are rare, occurring in \overallMinPerc--\overallMaxPerc\% of populations; most frequently for birds (\AvesRangePerc\%) followed by mammals (\MammaliaRangePerc\%), insects (\InsectaRangePerc\%), and fishes (\OsteichthyesRangePerc\%). When they occur, they tend to be driven by climate and severe winters, cycles of parasites and predators, and interactions between these elements. We find little evidence of intrinsic life-history correlates. *close with modelling and policy implications*

# Introduction

1. Black swans are statistically improbable events whose occurrence has major ramifications. One of the most profound black swans in ecology was the asteroid marking the mass extinctions at the K-T boundary. Today, it is the extremeness of climate that is expected to cause the greatest societal damage. While extremes in the physical environment are present and widely understood, it remains unclear the extent to which ecological systems buffer or suffer from black swans.

2. What is the existing evidence for heavy tails in ecological time series?

3. There are two main possible reasons for why we may find little evidence of ecological black swans. (1) don’t exist, (2) a paucity of appropriate methods or length/quality of time series.

4. Here we develop …our research questions and main conclusions.

Papers to work in:

[@inchausti2002; @halley2002; @inchausti2001]

<!--[@jentsch2007]-->

[@ward2007]

[@garcia-carreras2011] [@sornette2009]

[@nunez2012]

[@thompson2013] [@beaugrand2012] [@pine-iii2009]

[@doak2008]

[@smale2013]

[@easterling2000] [@scheffer2003] [@katz2005]

[@taleb2007]

[@vasseur2014]

[@vert-pre2013] [@lindenmayer2010]

[@valpine2002] [@gregory2010] @garcia-carreras2011 @brook2006 [@herrandoprez2014]

[@sibly2005; @ziebarth2010]

# Methods

1. Short overview paragraph

2. The data: global population dynamics database; briefly introduce the GPDD and how we filtered it. Expand in the supplement. Mention the breakdown of taxonomy and interpolation (Table S1).

3. The heavy-tailed Gompertz population model. Mention the simulation testing but don’t go into details here.

4. Mention the alternative population models we fit.

5. Modelling the possible correlates (mixed effects models with binomial or beta distributions and logit links; nested random effects for taxonomic class and taxonomic order). Life-history data from @brook2006a.

The heavy-tailed Gompertz model is:

$$\ln N_t = x_t = \mathrm{Student\mhyphen t}( {\color{red} \nu }, \lambda + b x_{t-1}, \sigma),$$

where $x_t$ is the $\ln$ abundance $N$ at time $t$. The model is density independent if $b = 1$, maximally density dependent if $b = 0$, and inversely density dependence if $b < 0$. The parameter $\lambda$ represents the expected population growth rate at $x_t = 0$. The process noise is modelled as a Student-t distribution with a scale parameter of $\sigma_\mathrm{proc}$, and a degrees of freedom parameter of $\nu$ (highlighted in red). If $\nu$ is small ($\lesssim 10$) the distribution has much heavier tails than a normal distribution. For example, at $\nu = 2$, the probability of drawing a value less than -5 is 1.8%, whereas the probability of drawing such a value from a normal distribution is nearly zero ($2.9\cdot10^{-5}$%). As $\nu$ approaches infinity the distribution approaches the normal distribution (Fig. \ref{fig:didactic}). By estimating the value of $\nu$ we can quantity how heavy-tailed the process noise deviations are.

We fit alternative models that allowed for autocorrelation of the residuals, allowed for observation error, and allowed the functional form of the population dynamics to be represented as Ricker-logistic (Supporting Material).

Supporting material to be noted:

1. Table showing the taxonomic breakdown, number of interpolated points, etc (Table S1).

2. Plot showing time series we assumed were recorded as log transformed (Fig. \ref{fig:log10-assumed}).

3. Plot showing all the time series we used with interpolated and zero imputed values (Fig. \ref{fig:all-ts}).

4. Plot showing the priors (Fig. \ref{fig:priors}).

# Results

Main figures and tables:

1. Fig. \ref{fig:didactic} illustrates the method: it shows the t-distribution tails, example simulated time series, and model fits to those time series.

2. Fig. \ref{fig:nu-coefs} shows the posterior distributions of the estimated $\nu$ values. These are split by taxonomic class and order.

3. Fig. \ref{fig:correlates} shows possible biological and time-series-property correlates of heavy-tailed behaviour.

4. Table. \ref{tab:sparks} selects populations that were categorized as heavy-tailed and digs into the causes. These are a non-random sample that I was able to verify in the literature.

Supporting material to be noted quickly:

1. Simulation testing: the ability to recover $\nu$ when randomly sampling from various distributions (Fig. \ref{fig:sim-nu}); heavy-tailed Gompertz model performance and confidence interval coverage given the process noise has deviations that are known to have effective $\nu$ equal to true $\nu$ (Fig. \ref{fig:sim-gompertz}); boxplots of the same output (Fig. \ref{fig:sim-gompertz-boxplots}); probability that $\nu < 10$ for the Gompertz simulation testing (Fig. \ref{fig:sim-prob}).

2. Time series of all heavy tailed populations (Fig. \ref{fig:heavy-ts}).

3. The effect of alternative population dynamics models on $\nu$ estimates (Fig. \ref{fig:alt}).

4. Coefficients from modelling covariates of the probability of heavy tails (Fig. \ref{fig:correlate-coefs}).

5. Example Stan code for heavy-tailed Gompertz model with AR1 residuals and observation error (a model including everything).

6. GPDD IDs used.

<!--Modelling result: an increase of 1 time step of data (given that you are ballpark around the mean — 30 time steps — to start with) results in an approximately 1% increase in the expected probability of observing heavy tails. (Using Gelman’s ‘divide by 4’ rule for interpreting logistic regression coefficients.)-->

# Discussion

1. Summary of our findings

2. How do our results mesh with previous related analyses?

3. Why might we expect to see heavy tails in ecological time-series? (mixture of normals or extreme drivers)

4. Are the observed frequencies by taxonomic class real or an observational phenomena?

5. Ways forward

6. Policy implications

# Acknowledgements

Funding: SCA: SFU Graduate Fellowship, NKD: NSERC SFU ..., ABC: ..., TAB: ...

Earth to Ocean research group for helpful discussions

Global Population Dynamics Database

Silhouettes: `phylopic.org`: rabbit by Sarah Werning, grey heron by Ardea cinerea, hoverfly by Gareth Monger. All under Creative Commons Attribution 3.0 Unported license.

Compute Canada’s WestGrid high-performance computing resources

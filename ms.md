*Main message, to be deleted (25 words)*: Black swans are present & taxonomically widespread but rare in population dynamics. Extreme climate, predation, parasites, and their interactions are common causes; intrinsic drivers are unclear.

# Abstract

Black swans are statistically improbable events that nonetheless occur --- often with profound implications. While extremes in the physical environment, such as monsoons and heat waves, are widely studied and increasing in magnitude and frequency, it remains unclear the extent to which ecological populations buffer or suffer from extremes. Here, we develop a probability model to estimate the degree of heavy-tailedness (presence of black swans) in ecological process noise. We apply our model to \NPops\ time series from around the world across \NOrders\ taxonomic orders and seven classes. We find strong evidence of black swans, but they are rare, occurring in \overallMinPerc--\overallMaxPerc\% of populations; most frequently for birds (\AvesRangePerc\%) followed by mammals (\MammaliaRangePerc\%), insects (\InsectaRangePerc\%), and fishes (\OsteichthyesRangePerc\%). When they occur, they tend to be driven by climate and severe winters, cycles of parasites and predators, and interactions between these elements. 
They are more frequently detected for populations with longer time series and lower levels of process noise.
We find little evidence of intrinsic life-history correlates. *close with modelling and policy implications*

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
To obtain estimates of the probability and magnitude of population dynamic black swans, we fitted population dynamic models to hundreds of abundance time series from around the world.
The key addition of our models was estimating the shape of the process noise tails --- how extreme the largest jumps from time step to time step were beyond what a standard population model allows for.
We then compared how these estimates related to taxonomy, time series properties, and life-history characteristics, and for a subset of populations, we documented the established causes of ecological black swans.

## Time-series data
We obtained abundance time series from the Global Population Dynamics Database (GPDD) [@gpdd2010]. The GPDD contains nearly 5000 time series of abundance across $\sim$ 1000 species and 100 taxonomic orders. To derive a high-quality subset of populations that were suitable for our analyses, we filtered the data (Supporting Material). Our filtering procedure removed populations from less reliable data sources, removed those without sufficient data, and interpolated some missing values. Our final dataset contained \NPops\ populations across \NOrders\ taxonomic orders and seven taxonomic classes with a median of XX time steps (range of XX--XX) (Table S1, Fig. \ref{fig:all-ts}).

## Population models
Our main analysis focusses on the commonly applied phenomenological Gompertz population dynamics model. In previous analyses of the GPDD, the Gompertz has proved to be the most frequently parsimonious model [@brook2006]. We extended the traditional Gompertz model by allowing the process noise to be drawn from a Student-t distribution with degrees of freedom parameter $\nu$. If $\nu$ is small ($\lesssim 10$) the t distribution has much heavier tails than a normal distribution. For example, at $\nu = 2$, the probability of drawing a value more than five standard deviations below the mean is 1.8%, whereas the probability of drawing such a value from a normal distribution is nearly zero ($2.9\cdot10^{-5}$%). As $\nu$ approaches infinity the distribution approaches the normal distribution (Fig. \ref{fig:didactic}). By estimating the value of $\nu$ we can quantity how heavy-tailed the process noise deviations are. 

The heavy-tailed Gompertz model is:
$$\begin{aligned}
x_t &= \ln N_t\\
x_t &= \lambda + b x_{t-1} + \epsilon_t\\
\epsilon_t &\sim \mathrm{Student\mhyphen t}(\nu, 0, \sigma),
\end{aligned}$$
where $x_t$ is the $\ln$ abundance ($N$) at time $t$. The model is density independent if $b = 1$, maximally density dependent if $b = 0$, and inversely density dependence if $b < 0$. The parameter $\lambda$ represents the expected population growth rate at $x_t = 0$. The process noise is modelled as a Student-t distribution with degrees of freedom parameter $\nu$, a mean of $0$, and a scale parameter of $\sigma$. 

We chose weakly informative priors to incorporate our understanding of plausible population dynamics [@gelman2014, Supporting Material]. For $\nu$, we chose an exponential prior with rate parameter of $0.01$ truncated at values above two (Fig. \ref{fig:priors}), a slightly less informative prior than that suggested by @fernandez1998. This prior gives only a 7.7% probability that $\nu < 10$ but constrains the sampling sufficiently to avoid wandering off towards infinity --- above approximately $\nu = 20$ the t distribution is so similar to the normal distribution that time series of the length considered here are unlikely to be informative about the precise value of $\nu$ (Fig. \ref{fig:didactic}).

We estimated all models in a Bayesian framework using Stan [@stan-manual2014]. Stan samples from the posterior distribution with a version of Hamiltonian Markov chain Monte Carlo called the No-U-Turn Sampler [@hoffman2014]. We assured that chains had sufficiently converged and the sampler had obtained sufficient independent samples from the posterior ($\widehat{R} < 1.05$, $n_\mathrm{eff} > 200$; Supporting Material). We used simulated data to test how how easily we could detect $\nu$ given different sample sizes and to ensure we could recover unbiased parameter estimates from the Gompertz model (Supporting Material, Figs \ref{fig:sim-nu}, \ref{fig:sim-gompertz}, \ref{fig:sim-gompertz-boxplots}).

We fit four alternative population models to test if they systematically changed our conclusions. (1) Autocorrelation has been suggested as a reason for increased observed variability of abundance time series through time, which could create apparent heavy tails [@inchausti2002]; therefore, we fit a model that modelled the serial correlation in the residuals. (2) Previous work has focussed modelled growth rates without accounting for density dependence [@segura2013]; therefore, we fit a simpler model modelling growth rate as a random walk with drift. (3) Observation error could conceivably create false heavy tails or mask our ability to detect heavy tails [@ward2007]. However, simultaneously estimating observation error and process noise in a state space framework is a challenging computational problem even with longer time series and without estimating the shape of the process error tails (REF); therefore, we fit a model where we assumed a level of observation error ($0.2$ standard deviations on a log scale). (4) The Gompertz model assumes that per capita growth rate declines linearly with log abundance. We also fit an alternative model, the Ricker-logistic model, which assumes that per capita growth rate declines linearly with abundance itself (Supporting Material).

## Covariates of population dynamic black swans
We investigated possible covariates of heavy-tailed population dynamics visually and through multilevel modelling. We plotted characteristics of the time series ($\sigma$, $\lambda$, $b$, and time-series length) along with two life-history characteristics (body length and maximum lifespan obtained from @brook2006a) against our estimates of $\nu$. We formally investigated these relationships by fitting beta regression multilevel models to the probability that $\nu$ was less than $10$ (i.e. approximately the probability of heavy tails). For covariates that were derived from Gompertz model parameter estimates, we incorporated standard deviations around the means. To account for broad patterns of phylogenetic relatedness, we allowed for varying intercepts at the taxonomic class, order, and species level (Supporting Material).

Finally, we investigated a sample of populations that our method categorized as having a high probability of heavy tails (Pr$(\nu < 10) > 0.5$). Where possible, we found the documented causes of ecological black swans in the primary data source cited in the GPDD or in other literature describing the population. These populations were chosen haphazardly with the purpose of generating hypotheses about the causes of population dynamic black swans.

# Results
We found strong but rare evidence for black-swan dynamics. Defining black-swan dynamics as a greater than $0.5$ probability that $\nu < 10$, we observed black swans most frequently for birds (X%) followed by mammals (X%), insects (X%), and fishes (X%) (Figs \ref{fig:nu-coefs}, \ref{fig:heavy-ts}). Black swans were taxonomically widespread, occurring for at least one population in X% of the taxonomic orders recorded. Accounting for time series length and partially pooling inference across taxonomic class and order with a multilevel model, there was stronger evidence for black swans in insect populations than is apparent in Fig. \ref{fig:nu-coefs} --- four of 10 orders with the highest median probability of heavy tails were insect orders --- however, there was considerable uncertainty in these estimates (Fig. \ref{fig:order-estimates}).

The majority of our heavy-tailed estimates were robust to alternative population models and observation error (Fig. \ref{fig:alt}). Including an autocorrelation structure in the residuals, modelling population growth rates as a random walk with drift, or modelling the population dynamics as Ricker-logistic did not systematically alter our conclusions. Allowing for a fixed quantity of observation error decreased our estimates of black-swan dynamics, increasing the median estimate of $\nu$ from $<10$ to $>10$ in X of XX populations although the majority of $\nu$ estimates remained similar (Fig. \ref{fig:alt})

Across populations, the probability of observing black swan dynamics was positively related to time-series length and negatively related to magnitude of process noise ($\sigma$) but not clearly related to population growth rate ($\lambda$), density dependence ($b$), or maximum lifespan (Figs \ref{fig:correlates}, \ref{fig:correlate-coefs}). Longer time-series length was the strongest covariate of observing black swan dynamics (as was similarly observed in simulation testing, Fig. \ref{fig:sim-nu}). For example, we were about 15% TODO more likely to observe heavy tails in a population with 40 time steps than in one with 30 time steps (Fig. \ref{fig:perc-inc-p}). However, the absolute change in probability with increased time series length was small ($0.X$ vs. $0.X$ in the previous example, Fig. \ref{fig:correlates}).

We assembled the documented causes of black swan population dynamics in seven populations with high probabilities of heavy tails (Table 1). The majority of documented events were downward black swans and the majority involved a confluence of events. For example, a shortage of nest sites for the shag in the UK was thought to reduce population productivity and make the population vulnerable when a red-tide hit in 1968. The population experienced a black swan downswing, but the new availability of nest sites resulted in a rapid population upswing [@potts1980]. As another example, an interaction between population cycles caused by predation combined with cycles caused by the environment are thought to have caused a downward black swan for a water vole population [@saucy1994]. Other black swans were the result of a sequence of extreme climate events on their own. For example, severe winters in 1929, 1940--1942, and 1962--1963 were associated with black swan downswings in grey heron abundance in the UK [@stafford1971]. Our analysis finds that the last event was a combination of two black swan events in a row and population recovery took three times longer than predicted [@stafford1971].

# Discussion

1. Summary of our findings
  
	<!--there is evidence of heavy tails/black swans... in ecological time series... but they are rare-->
	<!--Taxonomically wide spread-->
	<!--we don't find compelling... intrinsic obvious biological reason-->
	<!--commonly extrinsic extreme events and interactions - climate and severe winters, cycles of natural enemies (parasites + predators)-->
	<!--Estimating the degrees of freedom parameter in a t distribution is a viable method of quantifying black swan events in ecological time series-->

2. How do our results mesh with previous related analyses?

We might expect to observe black swan dynamics in ecological time series because of unmodelled intrinsic properties of populations or extrinsic forces acting on populations. A t distribution can be formed through a mixture of normal distributions (in which the variances are inverse-gamma distributed) [@gelman2014]. Therefore, if we miss some underlying mixture of processes we might expect to observe heavy tails. That process might be an aggregation of populations across space, or population diversity, or some intrinsic change in population variability through time. Extrinsic forces could also cause black swan dynamics. These forces could be extreme themselves. For example, extreme climate, predation from or competition with other species experiencing black swans, or sharp changes in human pressure such hunting, fishing, or habitat destruction. Alternatively, the interaction of multiple "normal" forces could give rise to black swan ecological dynamics. This could occur if the interaction nature is synergistic or even if non-synergistic interactions experience a rare alignment [@denny2009].

Are the observed frequencies by taxonomic class real or an observational phenomena?

Bring in dynamics-observation scale mismatch issue. Also note that variance and heavy tails are not the same thing, in fact they are somewhat inversely related in the same population.

Taxonomically and geographically biased sample and different kinds of census data are used within taxonomic groupings.

Scale mismatch hypothesis: 
The details are most likely to be observed if the scale of observation matches the scale of population dynamics
If we observed frequently relative to generation time (e.g. many large bodied mammals) we will average across many generations and perhaps not observe the tails
If we observe infrequently relative to generation time (e.g. many insects) we may miss tails or aggregate across population events
The mostly salmon populations present an enigma then --- population scale (annual) is right, low observation error, but low evidence of heavy tails, but the combination of different cycles creates high lambda and, which may allow for more extreme events.


Ways forward. 
Can we forecast the probability of black swans in space and time? 
How can we detect them quickly after they happen?
Move from phenomenological to mechanistic (e.g.\ recruitment) models. 
Investigation into the mechanisms and covariates in a geographic and taxonomic subsets where higher quality data are available.
What is the impact of allowing for black swans in forecasting of ecological risk?

6. Policy implications

# Acknowledgements

We thank members of the Earth to Ocean research group for helpful discussions and comments.
We are grateful to the contributors and maintainers of the Global Population Dynamics Database and to Compute Canada's WestGrid high-performance computing resources.
Three silhouette images were obtained from `phylopic.org`: rabbit (Sarah Werning), hoverfly (Gareth Monger), and bird (Jean-Raphaël Guillaumin [photography] and T. Michael Keesey [vectorization]) under a Creative Commons Attribution 3.0 Unported license.
Funding was provided by a Simon Fraser University Graduate Fellowship (SCA), NSERC (NKD, ABC), the Canada Research Chairs Program (NKD), and TAB TODO?.

\bibliographystyle{apalike}
\bibliography{/Users/seananderson/Dropbox/tex/jshort,/Users/seananderson/Dropbox/tex/ref3}

\clearpage 

# Tables

\input{sparks.tex}


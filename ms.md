*Main message, to be deleted (25 words)*: Black swans are present & taxonomically widespread but rare in population dynamics. Extreme climate, predation, parasites, and their interactions are common causes; intrinsic drivers are unclear.

# Abstract

Black swans are statistically improbable events that nonetheless occur --- often with profound implications. While extremes in the physical environment, such as monsoons and heat waves, are widely studied and increasing in magnitude and frequency, it remains unclear the extent to which ecological populations buffer or suffer from extremes. Here, we develop a probability model to estimate the degree of heavy-tailedness (presence of black swans) in ecological process noise. We apply our model to \NPops\ time series from around the world across \NOrders\ taxonomic orders and seven classes. We find strong evidence of black swans, but they are rare, occurring in \overallMinPerc--\overallMaxPerc\% of populations; most frequently for birds (\AvesRangePerc\%) followed by mammals (\MammaliaRangePerc\%), insects (\InsectaRangePerc\%), and fishes (\OsteichthyesRangePerc\%). When they occur, they tend to be driven by climate and severe winters, cycles of parasites and predators, and interactions between these elements.
They are more frequently detected for populations with longer time series and lower levels of process noise.
We find little evidence of intrinsic life-history correlates. *close with modelling and policy implications*

# Introduction

Black swans are unexpected extreme events with dramatic consequences. One of the most striking black swans in ecology was the asteroid marking the mass extinctions of the K-T boundary (REF). Today, it is the extremeness of climate --- particularly in concert with shifts in mean temperature --- that is expected to cause the greatest societal damage (REF). However, while extremes in the physical environment are present and widely accepted, it remains unclear the extent to which ecological systems buffer or suffer from black swans.

<!-- (re-write that) work in @jentsch2007 ("A new generation of climate-change experiments: events, not trends")-->

Existing evidence for black swans or heavy tails in population dynamics is thus far limited. Considerable work has focussed on heavy tails in physical environmental processes in, for example, rain fall (REF), wave height (REF), earth quake magnitude, climate temperature (REF). time to extinction, plankton time series, GPDD work found little evidence... subsequent work suggested, but ... Ward applied a data-intense method and found evidence for ...

There are two key reasons why we may find little evidence of ecological black swans. First, they might not exist. Indeed, the vast majority of population dynamics model fitting and risk forecasting assumes that population dynamics are not heavy tailed (e.g. REF, REF, REF). Alternatively, population dynamics might be heavy tailed, but a paucity of statistical methods and shortness of available time series has limited our ability to detect them. For example, a commonly applied probability distribution to detect the shape of time series tails (the generalized extreme value distribution) requires lumping data into the most extreme event per unit time (e.g. heaviest rainfall per year) and therefore requires a considerable quantity of data (REF). Previous work into the shape of population dynamic tails, has focused on model comparison to identify parsimonious explanations from a set of particular hypotheses, but does not quantify the magnitude of heavy tails on a continuous scale or put the conclusions in a probabilistic framework (REF). Furthermore, there is evidence that information criteria are poor at distinguishing subtle differences in population dynamics for the relatively short time series that are most readily available (REF).

- Or bring back in Ward - data intensive, can't be applied to most individual time series
- plankton --- hideous distribution to fit, even harder in population dynamics framework -and what about seasonal cycles, density dependence?
- GPDD - model competition and finding most parsimonious explanation from set of hypotheses - but only *very* heavy from not heavy
- GEV - super data intensive
- tipping point and regime shift methods - say something... not appropriate for most of the data available

<!--A large body of literature focusses on tipping points or regime shifts.-->

Here we develop a probabilist framework for identifying ecological black swans in population dynamic process noise --- the stochastic jumps from time step to time step. Our framework quantifies both the magnitude and probability of heavy-tailed dynamics, allows for a range of population dynamic models, can incorporate observation uncertainty, and is applicable to the lengths of abundance time series most commonly available. We apply our model to hundreds of populations from around the world across four taxonomic classes to address three questions: (1) how frequent and strong are black swans in population dynamics across taxonomy, (2) what are some verified causes of black swans in ecological time series, and (3) are there characteristics of time series or intrinsic life-history characteristics associated with heavy-tailed dynamics. We find a high probability of ecological black swans across a wide range of taxa but in a limited number of populations. When they occurred, extreme climate, predation, parasites, and their interactions were common documented causes. We found little evidence that intrinsic life history characteristics were related to the probability of black swans.

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

To obtain estimates of the probability and magnitude of population dynamic black swans, we fitted population dynamic models to hundreds of abundance time series from around the world. The key addition of our models was estimating the shape of the process noise tails --- how extreme the largest jumps from time step to time step were beyond what a standard population model allows for. We then compared how these estimates related to taxonomy, time series properties, and life-history characteristics, and for a subset of populations, we documented the established causes of ecological black swans.

## Time-series data

We obtained abundance time series from the Global Population Dynamics Database (GPDD) [@gpdd2010]. The GPDD contains nearly 5000 time series of abundance across $\sim$ 1000 species and 100 taxonomic orders. To derive a high-quality subset of populations that were suitable for our analyses, we filtered the data (Supporting Material). Our filtering procedure removed populations from less reliable data sources, removed those without sufficient data, and interpolated some missing values. Our final dataset contained \NPops\ populations across \NOrders\ taxonomic orders and seven taxonomic classes with a median of XX time steps (range of XX--XX) (Table S1, Fig. \ref{fig:all-ts}).

## Population models

Our main analysis focusses on the commonly applied phenomenological Gompertz population dynamics model. In previous analyses of the GPDD, the Gompertz has proved to be the most frequently parsimonious model [@brook2006]. We extended the traditional Gompertz model by allowing the process noise to be drawn from a Student-t distribution with degrees of freedom parameter $\nu$. If $\nu$ is small ($\lesssim 10$) the t distribution has much heavier tails than a normal distribution. For example, at $\nu = 2$, the probability of drawing a value more than five standard deviations below the mean is 1.8%, whereas the probability of drawing such a value from a normal distribution is nearly zero ($2.9\cdot10^{-5}$%). As $\nu$ approaches infinity the distribution approaches the normal distribution (Fig. \ref{fig:didactic}). By estimating the value of $\nu$ we can quantity how heavy-tailed the process noise deviations are.

<!--TODO cite [@knape2012] Are patterns of density dependence in the Global Population Dynamics Database driven by uncertainty about population abundance - used Gompertz - similar style paper-->

<!--see [@knape2013] introduction: "A fundamental estimation problem in state-space modeling is that the two sources of variation are often only weakly distinguished by the data (Millar and Meyer 2000, Dennis et al. 2006, Knape 2008). "-->

<!--another gompertz example: @dennis2014-->

The heavy-tailed Gompertz model is:
$$\begin{aligned}
x_t &= \ln N_t\\
x_t &= \lambda + b x_{t-1} + \epsilon_t\\
\epsilon_t &\sim \mathrm{Student\mhyphen t}(\nu, 0, \sigma),
\end{aligned}$$
where $x_t$ is the $\ln$ abundance ($N$) at time $t$. The model is density independent if $b = 1$, maximally density dependent if $b = 0$, and inversely density dependence if $b < 0$. The parameter $\lambda$ represents the expected population growth rate at $x_t = 0$. The process noise is modelled as a Student-t distribution with degrees of freedom parameter $\nu$, a mean of $0$, and a scale parameter of $\sigma$.

We chose weakly informative priors to incorporate our understanding of
plausible population dynamics [@gelman2014, Supporting Material]. For $\nu$, we chose an exponential prior with rate parameter of $0.01$ truncated at values above two (Fig. \ref{fig:priors}), a slightly less informative prior than that suggested by @fernandez1998. This prior gives only a 7.7% probability that $\nu < 10$ but constrains the sampling sufficiently to avoid wandering off towards infinity --- above approximately $\nu = 20$ the t distribution is so similar to the normal distribution that time series of the length considered here are unlikely to be informative about the precise value of $\nu$ (Fig. \ref{fig:didactic}).

We estimated all models in a Bayesian framework using Stan [@stan-manual2014]. Stan samples from the posterior distribution with a version of Hamiltonian Markov chain Monte Carlo called the No-U-Turn Sampler [@hoffman2014]. We assured that chains had sufficiently converged and the sampler had obtained sufficient independent samples from the posterior ($\widehat{R} < 1.05$, $n_\mathrm{eff} > 200$; Supporting Material). We used simulated data to test how how easily we could detect $\nu$ given different sample sizes and to ensure we could recover unbiased parameter estimates from the Gompertz model (Supporting Material, Figs \ref{fig:sim-nu}, \ref{fig:sim-gompertz}, \ref{fig:sim-gompertz-boxplots}).

We fit four alternative population models to test if they systematically changed our conclusions. (1) Autocorrelation has been suggested as a reason for increased observed variability of abundance time series through time, which could create apparent heavy tails [@inchausti2002]; therefore, we fit a model that modelled the serial correlation in the residuals. (2) Previous work has focussed modelled growth rates without accounting for density dependence [@segura2013]; therefore, we fit a simpler model in which we assumed density independence. (3) Observation error could bias parameter estimates [@knape2012] or mask our ability to detect heavy tails [@ward2007]. However, simultaneously estimating observation error and process noise in a state space framework is a challenging computational problem [e.g. dennis2006; @knape2008] even with longer time series and without estimating the shape of the process error tails; therefore, we fit a model where we assumed a level of observation error ($0.2$ standard deviations on a log scale). (4) The Gompertz model assumes that population growth rate declines linearly with log abundance. We also fit an alternative model, the Ricker-logistic model, which assumes that population growth rate declines linearly with abundance itself (Supporting Material).

## Covariates of population dynamic black swans

We investigated possible covariates of heavy-tailed population dynamics visually and through multilevel modelling. We plotted characteristics of the time series ($\sigma$, $\lambda$, $b$, and time-series length) along with two life-history characteristics (body length and maximum lifespan obtained from @brook2006a) against our estimates of $\nu$. We formally investigated these relationships by fitting beta regression multilevel models to the probability that $\nu$ was less than $10$ (i.e. approximately the probability of heavy tails). For covariates that were derived from Gompertz model parameter estimates, we incorporated standard deviations around the means. To account for broad patterns of phylogenetic relatedness, we allowed for varying intercepts at the taxonomic class, order, and species level (Supporting Material).

Finally, we investigated a sample of populations that our method categorized as having a high probability of heavy tails (Pr$(\nu < 10) > 0.5$). Where possible, we found the documented causes of ecological black swans in the primary data source cited in the GPDD or in other literature describing the population. These populations were chosen haphazardly with the purpose of generating hypotheses about the causes of population dynamic black swans.

# Results

We found strong but rare evidence for black-swan dynamics. Defining black-swan dynamics as a greater than $0.5$ probability that $\nu < 10$, we observed black swans most frequently for birds (X%) followed by mammals (X%), insects (X%), and fishes (X%) (Figs \ref{fig:nu-coefs}, \ref{fig:heavy-ts}). Black swans were taxonomically widespread, occurring for at least one population in X% of the taxonomic orders recorded. Accounting for time series length and partially pooling inference across taxonomic class and order with a multilevel model, there was stronger evidence for black swans in insect populations than is apparent in Fig. \ref{fig:nu-coefs} --- four of 10 orders with the highest median probability of heavy tails were insect orders --- however, there was considerable uncertainty in these estimates (Fig. \ref{fig:order-estimates}).

The majority of our heavy-tailed estimates were robust to alternative population models and observation error (Fig. \ref{fig:alt}). Including an autocorrelation structure in the residuals, modelling population growth rates as a random walk with drift, or modelling the population dynamics as Ricker-logistic did not systematically alter our conclusions. Allowing for a fixed quantity of observation error decreased our estimates of black-swan dynamics, increasing the median estimate of $\nu$ from $<10$ to $>10$ in X of XX populations although the majority of $\nu$ estimates remained qualitatively similar (Fig. \ref{fig:alt})

Across populations, the probability of observing black swan dynamics was positively related to time-series length and negatively related to magnitude of process noise ($\sigma$) but not clearly related to population growth rate ($\lambda$), density dependence ($b$), or maximum lifespan (Figs \ref{fig:correlates}, \ref{fig:correlate-coefs}). Longer time-series length was the strongest covariate of observing black swan dynamics (as was similarly observed in simulation testing, Fig. \ref{fig:sim-nu}). For example, we were about 15% TODO more likely to observe heavy tails in a population with 40 time steps than in one with 30 time steps (Fig. \ref{fig:perc-inc-p}). However, the absolute change in probability with increased time series length was small ($0.X$ vs. $0.X$ in the previous example, Fig. \ref{fig:correlates}).

We assembled the documented causes of black swan population dynamics in seven populations with high probabilities of heavy tails (Table 1). The majority of documented events were downward black swans and the majority involved a confluence of events. For example, a shortage of nest sites for the shag in the UK was thought to reduce population productivity and make the population vulnerable when a red-tide hit in 1968. The population experienced a black swan downswing, but the new availability of nest sites resulted in a rapid population upswing [@potts1980]. As another example, an interaction between population cycles caused by predation combined with cycles caused by the environment are thought to have caused a downward black swan for a water vole population [@saucy1994]. Other black swans were the result of a sequence of extreme climate events on their own. For example, severe winters in 1929, 1940--1942, and 1962--1963 were associated with black swan downswings in grey heron abundance in the UK [@stafford1971]. Our analysis finds that the last event was a combination of two black swan events in a row and population recovery took three times longer than predicted [@stafford1971].

# Discussion

We found strong evidence that the multiplicative process error jumps from time step to time step were greater than we'd expect under a normal distribution and a Gompertz or Ricker-logistic model in XX--XX% of populations assessed. We detected black swans more frequently for longer time series and for populations with less typical magnitude of process error, as we might expect. However, we failed to find strong evidence that black swans were associated with density dependence, population growth rate, or lifespan. In cases that we verified, black swans were often a result of the interaction between elements of extreme climate, predators and parasite cycles, and strong changes in human pressures. Our empirical results, sensitivity tests, and simulation tests suggest that estimating the tail shape of process error is a viable method of detecting heavy-tailed population dynamics. The presence of black swans highlights the importance of developing management strategies that detect quickly, respond to, and are robust to extremes in population dynamics --- particularly as the frequency and magnitude of climatic extremes increase over the next century.

2. How do our results mesh with previous related analyses?
- Not in disagreement with @halley2002 but different focus (@inchausti2002 and @inchausti2001 on reddening?)
- But tail events are necessarily rare and our time series generally very short
- @keitt1998 find heavy (power law) tails in breeding bird populations; but @allen2001 point out this is because of a mixture of normals --- feeds into our use of a t-distribution
- @segura2013 find heavy tails, but only modelling growth rates and used Levy stable distributions... a good example of a fruitful way of doing this kind of analysis --- longterm high frequency time series... could be mixture of seasons?
- plenty of evidence in environment: @katz2002, @katz2005, @denny2009
- in line with @doak2008 --- surprises inevitable?

We might expect to observe black swan dynamics in ecological time series because of unmodelled intrinsic properties of populations or extrinsic forces acting on populations. A t distribution can be formed through a mixture of normal distributions [in which the variances are inverse-gamma distributed, @gelman2014]. Therefore, if we miss some underlying mixture of processes we might expect to observe heavy tails. That process might be an aggregation of populations across space, or population diversity, or some intrinsic change in population variability through time. Extrinsic forces could also cause black swan dynamics. These forces could be extreme themselves. For example, extreme climate, predation from or competition with other species experiencing black swans, or sharp changes in human pressure such hunting, fishing, or habitat destruction might cause black swans. Alternatively, the interaction of multiple "normal" forces could give rise to black swan ecological dynamics. This could occur if the interaction nature is synergistic (REF) or even if non-synergistic interactions experience a rare alignment [@denny2009].

There a number of important caveats when considering the generality of our results. First, the GPDD data represent a taxonomically and geographically biased sample of populations. The longer time series we focus on are dominated by commercially and recreationally important species and a disproportionate number of populations are located in the United Kingdom. Although we expect that we would find strong but rare evidence for black swans for other taxa around the world, the common forces driving those black swans (e.g. severe winters in Table 1) likely differ. Second, in a large database with disparate datasources, we cannot ensure the quality of every datapoint. Recording mistakes undoubtedly comprise some of the black swans we detected. Conversely, recording mistakes could mask heavy tails if initial observations are discarded or altered as assumed recording errors. Third, the temporal scale of observation and population dynamics vary considerably across populations in the GPDD and this likely influences the detection of heavy tails. For example, if we make frequent observations relative to generation time (e.g. for many large bodied mammals) we will average across  generations and perhaps not observe the tails. Furthermore, if we record population censuses infrequently relative to generation time (e.g. many insects in the GPDD) the recorded data may average across extreme and less extreme events and dampen observed black swan dynamics.

Ways forward.
Can we forecast the probability of black swans in space and time?
How can we detect them quickly after they happen?
Move from phenomenological to mechanistic (e.g.\ recruitment) models.
Investigation into the mechanisms and covariates in a geographic and taxonomic subsets where higher quality data are available.
What is the impact of allowing for black swans in forecasting of ecological risk?

6. Policy implications

# Acknowledgements

We thank members of the Earth to Ocean research group for helpful discussions and comments. We are grateful to the contributors and maintainers of the Global Population Dynamics Database and to Compute Canada's WestGrid high-performance computing resources. Three silhouette images were obtained from `phylopic.org`: rabbit (Sarah Werning), hoverfly (Gareth Monger), and bird (Jean-Raphaël Guillaumin [photography] and T. Michael Keesey [vectorization]) under a Creative Commons Attribution 3.0 Unported license. Funding was provided by a Simon Fraser University Graduate Fellowship (SCA), NSERC (NKD, ABC), the Canada Research Chairs Program (NKD), and TAB TODO?.

\bibliographystyle{ecologyletters}
\bibliography{/Users/seananderson/Dropbox/tex/jshort,/Users/seananderson/Dropbox/tex/ref3}

\clearpage

# Tables

\input{sparks.tex}


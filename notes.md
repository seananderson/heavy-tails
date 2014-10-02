

1 - there is evidence of heavy tails/black swans... in ecological time series... but they are rare

- taxonomically widespread
- we don't find compelling... intrinsic obvious biological reason
- commonly extrinsic extreme events and interactions - climate and severe winters, cycles of natural enemies (parasites + predators),
- as climate gets more variable, extremes more frequent and greater magnitude (and skewed), might expect more of these
- and as our time series get longer and censusing gets more accurate... more ability to detect these
- calls for policy that is be robust to the possibility of black swans
- lead with ecological application

grouse are hunted - shags have erupted

2 - method
- probabilistic
- you can get away with 'these kinds of mistakes if you allow for heavy-tailed process errors'
-  first with a population dynamics model
- para on why novel

hook - to what degree to extreme environmental events..
population dynamics can often be connected to
will be increase in...
and then, however it remains unclear the degree to which animal populations also exhibit (too neutral, something like suffer) from extreme events...

ecology letters

- argue for why we do what we're doing

some we were able to track down

- digg into Barbary macaque

decreasing median - order plot

4 rows one column

add little years on spark lines at beginning and end at bottom?

switch to nu = 20 (or 10?) in first figure
l
combine fig 2 and 3 into new fig 2

new fig3 is cross plots
b, lambda, sigma_proc, N,
log10length, log10Lifespan or min age at maturity,

add 'black swan' simulation in which single value chosen at deviation = 8 or something like that

new supp fig that shows cross plot of accounting or not for observation error

show boxes of misclassification in cross plot figure of various model estimates

line up black swans, say with birds in UK - NAO

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

citation with good justification of gompertz and other references: \citep{herrandoprez2014}

The Gompertz was by far the most favoured model in \citet{brook2006} with the GPDD.

I could see an argument being made that we're just seeing them because of autocorrelation in the residuals. Autocorrelation plays a big role in \citet{inchausti2002} justifying why they think they're seeing increasing CV of time series in the GPDD with time and an increasing CV could create heavy tails. They just model the spectral properties and CV of the abundance time series --- not a population dynamics model. I tried estimating a coefficient for the first-order autocorrelation of the residuals. There is some autocorrelation it picks up sometimes, but interestingly this doesn't change things that much with respect to heavy tail results. Nothing systematic at least (Another important reference on spectral analysis and the GPDD: \citet{garcia-carreras2011}.)

little evidence for benefit of including annual climate or climate extremes with these data:
\citep{gregory2010}
and little evidence for allee effect with them

the natural logorithm version of the logistic is often convenient for parameter estimation
\citep{valpine2002} (and they cite a couple people too)

cite Limited evidence for the demographic Allee effect from numerous species across taxa
Stephen D. Gregory1,4, Corey J. A. Bradshaw2,3, Barry W. Brook2, and Franck Courchamp1
Read More: http://www.esajournals.org.proxy.lib.sfu.ca/doi/full/10.1890/09-1128.1

###

\begin{enumerate}

\item illustration of how the exponential prior slightly constrains the sampling of $\nu$\ldots illustrating what the prior probability of heavy tails is given uninformative data (about 5\% probability $\nu < 10$.

\item plots of the priors overlayed with estimates across populations for all parameters

\item cross plot of $\nu$ estimates from the models: no density dependence, add density dependence, add AR1, add assumed observation error

\item time series plots of all heavy-tailed populations

\item simulation testing of detecting nu given reduced sampling from true distribution; show repeated estimates of nu with different sample sizes, and show a few example samples with different observation window lengths and highlight the deviations that are beyond the 0.001 and 0.999 probabilities for a normal distribution.

\item simulation testing of the Gompertz model with process error deviations fixed to have an effective $\nu$ estimate at the true level --- this asks how well is the Gompertz model able to partition the various parameters in a precise and unbiased way given that we know the heavy tails are there (also include a scenario with a massive black swan to show how it behaves there)

\end{enumerate}

%Also, we may not be crazy seeing the heavier tails for heavier longer lived animals. There's some precedence for that in a paper that found more red noise (longer term fluctuations dominating the variance and therefore making the CV grow with time) in large-bodied animals. They attributed it to a mismatch between dynamics and the observation scale. They used a different life-history dataset too. But, there are other systematic differences with the body size-nu relationship that make this sketchy for us. E.g.\ taxonomic orders and interval frequency of data collection vary systematically and in clumps with body size (not shown in the below figure). And basically the only stuff with very \textit{low} probability of heavy tails are the the small ones with generation times less than a year and small collection intervals. It's basically a step function. We may not be able to tease that all apart... at least not in this paper.

black swans:

Pederson, N., Dyer, J.M., McEwan, R.W., Hessl, A.E., Mock, C.J., Orwig, D.A. et al. (in press). The legacy of episodic climatic events in shaping temperate, broadleaf forests. Ecol. Monogr., 10.1890/13– 1025.1.
Cohn, J.P. (2000). Saving the Salton Sea: researchers work to understand its problems and provide possible solutions. Bioscience, 50, 295–301.

Gelman 2006: "Large but finite values of A represent prior distributions which we call “weakly informative” because, even in the tail, they have a gentle slope (unlike, for example, a half-normal distribution) and can let the data dominate if the likelihood is strong in that region. "

TODO - see this paper: http://onlinelibrary.wiley.com/doi/10.1046/j.1365-2656.2003.00738.x/full Population-level mechanisms for reddened spectra in ecological time series H. Resit Akçakaya1,*, John M. Halley1,2 andPablo Inchausti1,3



- why we think we see them in some taxa but not others (observation scale dynamics mismatch; observation errors...)

Ways forward:

- hierarchical modelling of nu and other parameters

- better datasets for specific taxa (e.g.\ fish recruitment)

- joint prior on scale and nu parameters

- formally investigating what heavy tails are associated with

- do heavy tails have conservation importance?

- can we forecast their probability and develop correlates?

- how can we detect them relatively quickly after they happen?

- how often are heavy tails real? observation error, recording errors...

- what does this mean for how we model time series (should we consider using the t as a default even if we fix nu?)

- what does this mean for conservation management (crashes might be more frequent and greater in magnitude than our typical models tell us) --- importance of establishing polices that are robust to this


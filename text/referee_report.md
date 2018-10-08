# Referee Report 
## Reviewer's Comments:
I have completed my review of this manuscript and I apologize to the authors for the long delay. My recommendation is that this manuscript be revised prior to publication and that the revised version also be reviewed. I would be happy to review the revised manuscript. In the remainder of this report, I give the reasons for this recommendation.

## Major comments.

I have a small number of comments that I consider to be significant and that could, in principle, affect the quantitative results of the manuscript. It is my recommendation that the authors be required to address these issues in any revised manuscript.

1. The author's describe the covariance matrix from the measurements. However, as the predictions are made using a finite simulation, there should be an analogous contribution to the covariance from the predictions themselves. As the simulated and observed volumes are not very dissimilar for many of the samples the authors study, this error on the theoretical predictions should not necessarily be negligible.

2. Figure 8 shows some features that indicate MCMC chains that are not very well converged (at least in the sense of making the particular statements that the authors aim to make). The bimodality of the univariate posteriors evident in some panels (for example, for Acen, the most important parameter studied in this analysis) indicates possible poor convergence as it is extremely unlikely that this bimodality is physical. Likewise, the analogous criticism holds for the "isthmuses" of probability evident in several of the bivariate posterior plots. This comment may seem knit-picky, but is important for several reasons. A standard convergence test, such as Gelman-Rubin, more or less checks that the mean of the posterior is converged. However, these posterior plots suggest that the low-probability tails of the posterior are not well converged. This is important because it is the shape of the posterior in these low-probability tails that determines things such as the errors on the inferred parameters and the compatibility of the inferred parameters with various hypotheses. Based on these figures, I would not be surprised if these quantitative statements could change with a larger MCMC sample. In the context of the present work, convergence is likely to be slow because each region of parameter space must be sampled many times due to the fact that the likelihood function itself is stochastic.

3. Related to point 2 above, the authors should specify how they determine the maximum likelihood (or minimum chi^2). The current manuscript does not describe this calculation, but it may have important consequences, especially for model comparison. The calculation of the maximum likelihood in this context is non-trivial for two reasons. First, MCMC is not particularly good at determining the maximum likelihood. Do the authors simply choose the maximum likelihood realized in their chains? That may not be a particularly good estimate of maximum likelihood, especially if the likelihood is strongly peaked or has features on a scale smaller than the 68% regions. It would be better to get close to the maximum likelihood using the results of the MCMC chain and then to find the maximum likelihood using a function maximization procedure. 

Second, the maximum likelihood is difficult to determine in this type of model because the likelihood itself is stochastic. The maximum likelihood should be estimated from an average of many likelihood calculations at a single point in the parameter space. 

These possibilities should be checked by the authors and the method for estimating the maximum likelihood should be specified.

4. I am rather concerned that the resolution of the Bolshoi-P and SMDP simulations is insufficient to analyze the samples with Mr > -19.5. This issue cannot remain unaddressed in the manuscript. My reasoning is as follows. For these samples, log Mmin is on the order of 11.5 to 11.7, while sigma_log_M is on the order of 0.6 or larger (much larger in some cases). In the models studied by the authors, a significant number of halos are populated that are ~ 2*sigma_logM from the value of log Mmin. That means that halos with log-masses of log(Mmin)-2*(sigma_log_M) must be well resolved. For the -18.5 sample (just as an example) this corresponds to a mass of ~10^10 solar masses! For these simulations, that is about 100 particles. This may be too few to determine clustering reliably in dense regions as poorly resolved halos will have their structures disturbed by neighboring halos to a much larger degree than well resolved halos, enhancing the "assembly bias effect." More specifically, determining halo concentration from halos with few particles is quite uncertain, and it is not clear how this additional uncertainty may alter the results quoted in the present manuscript. There are many papers on this subject, recent examples being Poveda-Ruiz et al. 2016 and Klypin et al. 2016, but see the references therein as well. The rule-of-thumb that is often used is that the halos should have a few thousand particles in order to have concentration determined reliably. Therefore, it is probably necessary to consider only halos with masses above ~10^11 solar masses in the present work. Alternatively, the authors may attempt to construct an argument that demonstrates that resolution does not alter their results and that this rule-of-thumb criterion may not be applicable in this case. 

## Minor comments.

> 1.The two-population model of assembly bias within the context of the halo model is only a special case of the decorated HOD models described in Hearin et al. 2016. This is not clear from the text of the manuscript, particularly around the 9th paragraph of section 1.

The 9th paragraph now explicitly states that the two-population model used in this investigation is a special case of the decorated HOD model from Hearin et al. (2016). 

> 2.In Section 2.1, the larger volumes are necessary because of lower number densities. The phrasing in this section (that some samples occupy "larger volumes") is strange and potentially misleading.

The phrasing of the section has been clarified. 

> 3.In Section 2.2, the authors state that they choose the concentrations of the halos from the Dutton & Maccio median relation. This is acceptable. However, this has the potential to mitigate one source of assembly bias that may (or may not) be present in nature. In particular, because individual halos cluster as a function of concentration, satellite populations may also be distributed within their halos in a manner that is correlated spatially. By using the Dutton & Maccio relation, the authors eliminate this possibility. This is likely to be a small effect, but this detail should at least be pointed out by the authors. 

The halo concentrations now come directly from the halo catalog.

> 4.Most of the interpretation of the clustering results, especially toward the beginning of section 5.1, follows directly from the modeling paper of Hearin et al. 2016. This work should be cited as the source for this interpretation.

> 5.In section 5.1, the authors attribute abundance matching by Vmax to Hearin and Watson 2013. These authors did use abundance matching, but did not introduce Abundance matching via Vmax. Abundance matching via Vmax goes back at least to a paper by Conroy et al. 2006. 

Conroy et al. (2006) has been included in the citation for abundance matching via V_max.

> 6.In comparing results from different simulations in section 5.3, the authors draw particular attention to the -20.5 and -19 samples. However, drawing attention to these particular samples seems hardly warranted given the uncertainties in each inference of Asat (e.g., in Fig. 6). Rather, these seem to be specific examples of the general trend that Asat inferred from Bolshoi-P tends to be marginally higher than Asat inferred from SMDP.

> 7.There are a number of awkward phrasings in the text. As an example, take the first sentence of Section 5.1. "The constraints on the assembly bias parameters fall into two main categories. First, the satellite assembly bias parameter Asat and the second, the central assembly bias parameter Acen." This is awkward because category really isn't the correct word here. There are simply two parameters. Second, the second "sentence" is actually a sentence fragment as it has no verb. This is one example of awkwardness that is sprinkled throughout the paper. The manuscript could stand to be proofread more thoroughly, both by the authors and by a native English speaker if neither of the authors speak English as a native language.

This comment is bit offensive... 


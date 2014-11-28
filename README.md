# bPCA - Bayesian Principal Components Analysis

*Jan Smyčka*, *Petr Keil*

### Description
This is an R package that enables to perform Bayesian Principal Components Analysis, which we call bPCA. 

### Installation
For easy installation directly from this repository **just type this in R**:

==	library(devtools)

	install_github("bPCA", username="petrkeil")
	
	library(bPCA)==

### The idea
The idea is to fit a **MultiVariate Normal (MVN) distribution** to set of continuous variables. The fitting is done using MCMC sampler in JAGS. The means and covariances (parameters of the MVN) are monitored during the MCMC sampling and stored as MCMC chains, and subsequently subjected to various summary procedures.

The potential advantages over classical PCA are:
- **Prior information**  about associations between variables can be provided.
- **Stability of the PCA** can be assessed, especially when only small sample sizes are available.  
- The **posterior distributions** for the PCA scores and loadings can be extracted for further use, e. g. anywhere where propagation of uncertainty is of interest.
  

### History
The project was initiated by **Jan Smyčka** who approached **Petr Keil**  with the idea of performing principal components analysis (PCA) in Bayesian setting. After some discussions Jan put together the first JAGS code and essentially the core function of the project. Petr then build the rest of the package, and together they wrote the documentation. 

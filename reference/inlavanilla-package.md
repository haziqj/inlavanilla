# inlavanilla: INLA-Style Skew-Normal Marginals for Any Log-Posterior

Approximates a posterior distribution from its log density alone, using
the approximation pipeline of the 'INLAvaan' package (Jamil and Rue,
2026a)
[doi:10.48550/arXiv.2603.25690](https://doi.org/10.48550/arXiv.2603.25690)
, an adaptation of the Integrated Nested Laplace Approximation (Rue et
al., 2009)
[doi:10.1111/j.1467-9868.2008.00700.x](https://doi.org/10.1111/j.1467-9868.2008.00700.x)
, with no latent field and nothing model-specific. A Laplace
approximation with a variational Bayes mean correction gives the log
marginal likelihood, a skew-normal is fitted to each parameter's
marginal by scanning the log-posterior, and a Gaussian copula couples
the marginals for posterior draws. See also Jamil and Rue (2026b)
[doi:10.48550/arXiv.2604.00671](https://doi.org/10.48550/arXiv.2604.00671)
.

## See also

Useful links:

- <http://haziqj.ml/inlavanilla/>

- <https://github.com/haziqj/inlavanilla>

- Report bugs at <https://github.com/haziqj/inlavanilla/issues>

## Author

**Maintainer**: Haziq Jamil <haziq.jamil@gmail.com>
([ORCID](https://orcid.org/0000-0003-3298-1010)) \[copyright holder\]

Authors:

- Haziq Jamil <haziq.jamil@gmail.com>
  ([ORCID](https://orcid.org/0000-0003-3298-1010)) \[copyright holder\]

Other contributors:

- Håvard Rue ([ORCID](https://orcid.org/0000-0002-0222-1881))
  (Statistical and computational methodology) \[contributor\]

- Thomas Luu (Skew-normal quantile approximation (Apache License 2.0))
  \[copyright holder\]

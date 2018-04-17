# LNA++: a Fast C++ Implementation of the Linear Noise Approximation with first- and second-order sensitivities

The linear noise approximation (LNA) provides an approximate description of the statistical moments of stochastic chemical reaction networks, a commonly used modeling paradigm describing the probability distribution of systems of biochemical species in the intracellular environment. Unlike exact formulations, the LNA remains computationally feasible even for CRNs with many reaction events. The tractability of the LNA for approximation of CRNs makes it a common choice for inference of unknown chemical reaction parameters. However, this task is impeded by a lack of suitable inference tools for arbitrary CRN models. In particular, no available tool provides temporal cross-correlations, parameter sensitivities and efficient numerical integration.

LNA++ provides the means for fast derivation and simulation of the LNA including the computation of means, covariances, and temporal cross-covariances of species in the CRN. For efficient parameter estimation and uncertainty analysis, LNA++ implements 1st and 2nd order sensitivity equations. Interfaces are provided for easy integration with Matlab and Python.

More detailed information is available in [Tutorial.docx](Tutorial.docx) and [UserGuide.docx](UserGuide.docx).

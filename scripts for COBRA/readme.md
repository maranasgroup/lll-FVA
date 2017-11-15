readme.md

This is a program associated with the publication 'Accelerating flux balance calculations in genome-scale metabolic models by localizing the applicaitons of loopless constraints' by Chan SHJ, Wang L, Dash S and Maranas CD (2017).
The functions in this folder provide a general implementation of the presented method in the COBRA toolbox. COBRA toolbox must be installed for the functions to be properly used.

`fluxVariabilityLLC.m` was adapted from the latest `fluxVariability.m` from the COBRA toolbox (as of Nov 10, 2017) with Fast-SNP (Saa and Nielsen, 2016) and the localized loopless FVA implemented. Please see the documentation therein.

`addLoopLawConstraints.m` was also adapted from the latest function with the same name from the COBRA toolbox with the new methods implemented.

`fastSNPcobra.m` is the Fast-SNP preprocessing implemented in COBRA toolbox

`findMinNull.m` is the new method to quickly find a minimal null-space for ll-FVA presented in the paper

`getRxnLink.m` is to find the reactions connected to each other by any EFMs in loops by calculating EFMs. It calls `CalculateFluxModes.m` which is a funciton of the EFMtool (Terzer and Stelling, 2008)

`testLooplessFVAmethods.m` tests `fluxVariabilityLLC.m` with a model shipped with the COBRA toolbox.


References:
Saa,P.A. and Nielsen,L.K. (2016) Fast-SNP: A fast matrix pre-processing algorithm for efficient loopless flux optimization of metabolic models. Bioinformatics, 32, 3807–3814.
Terzer,M. and Stelling,J. (2008) Large-scale computation of elementary flux modes with bit pattern trees. Bioinformatics, 24, 2229–35.

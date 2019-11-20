# PoissonBalancedNets
Implementation of Balanced Spiking Networks with Poisson Neurons, from Rullán &amp; Pillow 2019. [Preprint here.](https://www.biorxiv.org/content/10.1101/836601v1l)

=========================

MATLAB code to implement the balanced spiking networks (BSN) and generate figures from Rullán & Pillow, 2019.

### Sample scripts
*  `figure_X.m`: Scripts to generate figures 2, 4, 5, 6, 7 and 8 from Rullán & Pillow, 2019. (e.g., figure_2.m for fig. 2)

### BSN model functions
*  `local_framework.m`: implementation of the local Poisson BSN 
*  `population_framework.m`: implementation of the population Poisson BSN 
*  `bsn.m`: implementation of the original BSN from [Boerlin et al.. 2013](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003258)

### Other 
*  `plotRaster.m`: plots raster of spikes
*  `colorcet.m`: generates perceptually uniform colormaps (for fig. 5B and 5C)
*  `test_approx.m`: calculates theoretical minimum error from exponential approximation (for fig. 7B)


Downloading the repository
------------

- **From command line:**

     ```git clone git@github.com:pillowlab/PoissonBalancedNets.git```

- **In browser:**   click to
  [Download ZIP](https://github.com/pillowlab/PoissonBalancedNets/archive/master.zip)
  and then unzip archive.

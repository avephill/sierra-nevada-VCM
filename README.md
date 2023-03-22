# Vegetation-climate mismatch in Sierra Nevada conifer forests

Hill, A. P. et al. [Low-elevation conifers in California's Sierra Nevada are out of equilibrium with climate.](https://doi.org/10.1093/pnasnexus/pgad004) *PNAS Nexus* (2023)


## Analysis
Code that was used for primary analysis are in the `analysis` directory. Because of the size of the source data, it's not available here.

Wieslander survey data are available thanks to the efforts of the [VTM Digitization Project](http://vtm.berkeley.edu/#/data/vegetation). Contemporary vegetation data are available from US Forest Service [EVeg maps](https://data.fs.usda.gov/geodata/edw/datasets.php?xmlKeyword=eveg). Historical and contemporary climate data at 800 m resolution are available from the [PRISM Climate Group](prism.oregonstate.edu). Future climate data are available from [CMIP6](esgf-node.llnl.gov/projects/cmip6/)


## Data Products
The VCM map shown in Figure 2 of the manuscript is available in the `data` folder as [SierraNevadaVCM_2015-2020.gpkg](https://github.com/avephill/sierra-nevada-VCM/blob/main/data/SierraNevadaVCM_2015-2020.gpkg). This was originally a raster at 800m resolution but was vectorized here for convenience. 

If you're interested in projections to future time periods or different time windows (e.g. 2000-2020) please email me and I'll be happy to help.


## License
All script and data are available with CC-BY-SA license.

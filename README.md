# thesis-pia
Whitefish and vendace SDMs
 
In this master’s thesis we apply state-of-the-art statistical species distribution models to predict the likely effects of climate change to whitefish and vendace larvae in
the Gulf of Bothnia, located in the northern part of the Baltic Sea.

The repository contains the main code implemented in the studies. The files are related to:

* nc-file: SmartSea data preparation, with comments

* whitefish: preliminary data cleaning and various models implemented through the studies, with comments

* SmartSeaData prep: function used to recover final raster (can't be run directly, since the source data are too big to be uploaded on github )

* dataprep: data preparation, to save the variables used in the models

* functions: all functions implemented and used in the other files

* models: fit of the 8 models reported in the thesis

* modelselect: kfoldcv and training-test validation of the models

* whitefish-final: two final models with outputs

* Effect of bottom coverage: extra material, not included in the thesis

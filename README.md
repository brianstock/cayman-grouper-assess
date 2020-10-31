## Pulse recruitment and recovery of Cayman Islands Nassau Grouper (**Epinephelus striatus**) spawning aggregations revealed by **in situ** length-frequency data

This repository holds the analysis for: Stock, BC, SA Heppell, L Waterhouse, IC Dove, CV Pattengill-Semmens, CM McCoy, PG Bush, G Ebanks-Petrie, and BX Semmens. Pulse recruitment and recovery of Cayman Islands Nassau Grouper (**Epinephelus striatus**) spawning aggregations revealed by **in situ** length-frequency data. **ICES Journal of Marine Science**. In press.
 
To recreate the analysis, you can run the `.R` files in the `code` folder in order. Intermediate model results are saved in `results`. Figures will be produced in the `plots` folder.

In addition to the `.R` code files, there are three `.cpp` files with TMB models:

1. `growth_REisland.cpp`: hierarchical growth model inspired by [Helser and Lai (2004)](https://www.sciencedirect.com/science/article/abs/pii/S0304380004001577). We used this model to estimate island-specific growth parameters as random effects, and then use these estimates in the LBSPR and LIME-fixed-k assessment models.

2. `LIME.cpp`: we made some minor modifications to [LIME](https://github.com/merrillrudd/LIME) appropriate for Little Cayman Nassau Grouper (see text). LIME assumes that growth parameters are fixed/known, and in the text we refer to this as the "LIME-fixed-k" model.

3. `LIME_integrated_REisland.cpp`: we merged the above two models to conduct an integrated assessment where growth parameters were estimated internally, the "LIME-integrated" model. This propagates the uncertainty in growth parameters into estimates of recruitment, depletion, etc. 

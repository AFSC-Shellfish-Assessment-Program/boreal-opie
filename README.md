# Authors
M.A. Litzow, E.J. Fedewa, M.J. Malick, E.R. Ryznar

# License
Software code created by U.S. Government employees is not subject to copyright in the United States (17 U.S.C. §105). The United States/Department of Commerce reserve all rights to seek and obtain copyright protection in countries other than the United States for Software authored in its entirety by the Department of Commerce. To this end, the Department of Commerce hereby grants to Recipient a royalty-free, nonexclusive license to use, copy, and create derivative works of the Software outside of the United States.

# About
This repository includes all data and code required to reproduce the results presented in: Litzow, M.A, Fedewa, E.J., Malick, M.J., Connors, B.M., Eisner, L., Kimmel, D.G., Kristiansen, T., Nielsen, J.M., and Ryznar, E.R. 2024. Human-induced borealization leads to the collapse of Bering Sea snow crab. Nature Climate Change, doi: 10.1038/s41558-024-02093-0

If you have any questions or find any errors, please contact the lead author: mike.litzow@noaa.gov.

# Scripts
The following brief descriptions serve as a guide to the different scripts used in the analysis.

*BCS Prevalence.r* -- Processes visual bitter crab syndrome (BCS) detection data from the bottom trawl survey for use in the borealization index.

*borealization_attribution_resampling.R* -- Calculates preindustrial and historical probabilities for estimation the Fraction of Attributable Risk (FAR) and Risk Ratio (RR).

*borealization_effects.R* -- Evaluates the statistical relationship between southeast Bering Sea borealization and snow crab abundance.

*borealization_effects_additional_model_evaluation.R* -- Conducts model evaluation presented in Extended Data Figs. 2 & 7.

*borealization_pdfs.R* -- Calculates probability densities for borealization index at different levels of North Pacific warming.

*bottom_temp_imputation.R* -- Imputes missing 2020 bottom temperature values for model evaluation presented in Extended Data Fig. 2.

*bycatch_biomass_plot.R* -- Creates Extended Data Fig. 6.

*chla_processing.R* -- Processes phytoplankton size ratio data for inclusion in borealization index.

*DFA.R* -- Creates borealization index using a Dynamic Factor Analysis model.

*ERSST query.R* -- Creates Extended Data Fig. 1.

*estimate_bottom_temp.R* -- Uses multiple imputation to fill in missing station values for bottom temperature and corrects bottom temperature estimates for annual differences in seasonal timing of sampling.

*Figs.R* -- Creates Figs. 1 & 2.

*Groundfish_biomass.R* -- Processes groundfish data for inclusion in borealization index.

*ice_query.R* -- Processes ERA5 ice cover data for inclusion in borealization index.

*Imm_Area.R* -- Calculates core range for immature snow crab.

*imputed_male_female_abundance_Fig1.R* -- Uses multiple imputation to estimate missing station values for snow crab abundance for plotting in Fig. 1 and use in "borealization_effects.R".

*sst_borealization.R* -- Creates Bayesian operating model of SST-borealization index relationship for use in attribution.

*stan_utils.R* -- Runs diagnostics for Bayesian regression models.

*zooplankton_processing.R* -- Processes *Calanus* and *Pseudocalanus* data for inclusion in borealization index.
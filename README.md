# boreal-opie

Supporting adaptation decisions for eastern Bering Sea snow crab fishery stakeholders through attribution of extreme climate events and assessment of historical and forward-looking assessment of the risk of borealizing climate change.

Initial thoughts re. analysis steps:

- Demonstrate borealization within the observed record. Use Bayes DFA to assess shared trend in immature snow crab abundance and other time series associated with borealization, including arctic:subarctic biomass ratios from bottom trawl survey, cod overlap with juvenile core habitat, BCS visual incidence, SC juvenile abundance change from previous year (or just juvenile SC abundance), Arctic cod abundance, cold pool extent, sea ice extent, phytoplankton/zooplankton size structure

- If we find strong shared variability among these time series, use SST (annual or winter / large spatial scale, full Bering) as a covariate in a Bayes regression model with the borealization trend as a response variable (and we should be able to use the uncertainty around that shared trend in regression model)

- Force the same model with FAR values for SST

- Use CMIP6 projections to compare historical and forward-looking perspectives on the risk of borealization
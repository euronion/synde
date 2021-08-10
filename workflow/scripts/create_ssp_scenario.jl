# Activate julia environment
using Pkg
Pkg.activate("./workflow/envs")
using GlobalEnergyGIS

scenario = snakemake.wildcards["ssp_scenario"]
ssp_year = parse(Int64, snakemake.wildcards["ssp_year"])

# Prepare the scenario data
create_scenario_datasets(scenario, ssp_year)

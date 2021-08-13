# SPDX-FileCopyrightText: 2021 Johannes Hampp <johannes.hampp@zeu.jlug.de>
#
# SPDX-License-Identifier: MIT

# Activate julia environment
using Pkg
Pkg.activate("./workflow/envs")
using GlobalEnergyGIS

scenario = snakemake.wildcards["ssp_scenario"]
ssp_year = parse(Int64, snakemake.wildcards["ssp_year"])

# Prepare the scenario data
create_scenario_datasets(scenario, ssp_year)

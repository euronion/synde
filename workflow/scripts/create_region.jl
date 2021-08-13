# SPDX-FileCopyrightText: 2021 Johannes Hampp <johannes.hampp@zeu.jlug.de>
#
# SPDX-License-Identifier: MIT

# Activate julia environment
using Pkg
Pkg.activate("./workflow/envs")
using GlobalEnergyGIS
import YAML

region_name = snakemake.wildcards["region"]
region_definitions = YAML.load_file(snakemake.input["region_definitions"])

region = region_definitions["regions"][region_name]

# Ensure members are alphabetically sorted by code
region = sort(collect(region))

# GEGIS expects a matrix with name + GADM(...) entry for each
# member of a region. Convert to appropriate shape.
members = [[k GADM(v)] for (k,v) in region]
members = vcat(members...)

# Create and save the GEGIS region
saveregions(region_name, members)
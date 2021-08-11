# Activate julia environment
using Pkg
Pkg.activate("./workflow/envs")
using GlobalEnergyGIS
import YAML

region_name = snakemake.wildcards["region"]
region_definitions = YAML.load_file(snakemake.input["region_definitions"])
gadm_members = region_definitions["regions"][region_name]

# Ensure members are alphabetically sorted
gadm_members = sort(gadm_members)

# GEGIS expects a matrix with name + GADM(...) entry for each
# member of a region. Convert to appropriate shape.
members = [[name GADM(name)] for name in gadm_members]
members = vcat(members...)

saveregions(region_name, members)
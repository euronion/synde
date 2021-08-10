# Activate julia environment
using Pkg
Pkg.activate("./workflow/envs")
using GlobalEnergyGIS

region = snakemake.wildcards["region"]
ssp = snakemake.wildcards["ssp"]
rcp = snakemake.wildcards["rcp"]
ssp_year = parse(Int64, snakemake.wildcards["ssp_year"])
era_year = parse(Int64, snakemake.wildcards["era_year"])

ssp_rcp_scenario = replace("$ssp-$rcp", "." => "") # GEGIS expects RCP without "."

# Use GEGIS ML approach to generated synthetic demand
predictdemand(gisregion=region,
        sspscenario=ssp_rcp_scenario,
        sspyear=ssp_year,
        era_year=era_year)

# Make sure output directory exists
mkpath(dirname(snakemake.output[1]))

# Move the output file from the default output location to a better location
mv("./resources/gegis/output/SyntheticDemand_" * region 
        * "_" * ssp_rcp_scenario 
        * "-" * string(ssp_year) 
        * "_" * string(era_year) 
        * ".jld",
   snakemake.output[1])
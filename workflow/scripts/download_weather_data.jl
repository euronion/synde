# SPDX-FileCopyrightText: 2021 Johannes Hampp <johannes.hampp@zeu.jlug.de>
#
# SPDX-License-Identifier: MIT

# Activate julia environment
using Pkg
Pkg.activate("./workflow/envs")
using GlobalEnergyGIS

era_year = parse(Int64, snakemake.wildcards["era_year"])

println("Downloading ERA5 data.")
download_and_convert_era5(era_year)

# SPDX-FileCopyrightText: 2021 Johannes Hampp <johannes.hampp@zeu.jlug.de>
#
# SPDX-License-Identifier: MIT

# Activate julia environment
using Pkg
Pkg.activate("./workflow/envs")
using GlobalEnergyGIS

println("Downloading aux. datasets.")
download_datasets()

println("Rasterising datasets.")
rasterize_datasets(cleanup=:all)

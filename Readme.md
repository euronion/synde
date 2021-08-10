

## Installation

1. Install `julia` (version 1.6.2 or above was used for this workflow),
    see [instructions and download here](https://julialang.org/downloads/platform/).

> Important:  
> Make sure `julia` can be executed from the shell, e.g. add `julia` to your `PATH` environment
> as described in the documentation.

2. From within the `SynDe` project directory open a shell and startup `julia` with the project environment enabled:

    ```
    SynDe$ julia --project=workflow/envs/
    ```
    Open the package manager to check if the environment is activated
    and build the packages (may take some time)
    ```
    julia> ]
    
    (envs) pkg> status
        Status `~/share/GitHub/SynDe/workflow/envs/Project.toml`
    [31bfc850] GlobalEnergyGIS v0.1.0 `https://github.com/niclasmattsson/GlobalEnergyGIS#master`
    
    (envs) pkg> instantiate
    ```
3. Configure `GlobalEnergyGIS` with user account data and work folder location
    ```
    julia> using GlobalEnergyGIS

    julia> saveconfig(string(pwd(),"/","resources/gegis"), <your CopernicusID>, "<your Copernicus API key>", agree_terms=true)
    ```
    For instruction what to enter for `<your ...>` see the [GEGIS documentation](https://github.com/niclasmattsson/GlobalEnergyGIS#2-create-config-files-and-agree-to-dataset-terms).




## Snakemake workflow

Important notes:

* In the used version (2021-08-09) during the `download_datasets()` step,
    the protected land database (WDPA) did not download correctly, possibly due to a wrong URL.
    Changing the URL inside the `GlobalEnergyGIS` code to
    ```
    "WDPA"         ("WDPA (protected areas):", "WDPA.zip",
    "https://d1gam3xoknrgr2.cloudfront.net/current/WDOECM_Aug2021_Public_shp.zip")
    ```
    solved the issue.
* In the used version (2021-08-09) during the `download_datasets()` step,
    the extraction step of the WRI database did not work properly.
    The files were extracted but failed to be moved to the correct location,
    this was then done manually and the faulty step removed from the `GlobalEnergyGIS`,
    then the step was repeated.


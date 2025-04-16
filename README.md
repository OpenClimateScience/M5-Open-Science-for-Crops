M5: Open Science for Crop Conditions
====================================

> How can NASA datasets be used to map crop conditions?

The fifth module of our [open climate-science curriculum](https://openclimatescience.github.io/curriculum) focuses on how to begin a reproducible computational science project, using crop conditions as a thematic example. **At the end of this module, you should be able to:**

- Access and utilize satellite-based datasets on plant productivity, condition, and evapotranspiration.
- Understand how plants link between the carbon and water cycles through photosynthesis and evapotranspiration.


Getting Started
---------------

[See our installation guide here.](https://github.com/OpenClimateScience/M1-Open-Climate-Data/blob/main/HOW_TO_INSTALL.md)

You can run the notebooks in this repository using [Github Codespaces](https://docs.github.com/en/codespaces/overview) or as [a VSCode Dev Container](https://code.visualstudio.com/docs/devcontainers/containers). Once your container is running, launch Jupyter Notebook by:

```sh
# Create your own password when prompted
jupyter server password

# Then, launch Jupyter Notebook; enter your password when prompted
jupyter notebook
```

**The Python libraries required for the exercises can be installed using the `pip` package manager:**

```sh
pip install xarray netcdf4 dask
```


Acknowledgements
----------------

This curriculum was enabled by a grant from NASA's Transition to Open Science (TOPS) Training program (80NSSC23K0864), part of [NASA's TOPS Program](https://nasa.github.io/Transform-to-Open-Science/)


# ## Ideas from QuantumToolbox.jl
# ### Language-specific (Julia) quirks (for users new to the language)
# ### Technical info, e.g. AutoDiff support for e.g. QuantumSystems and rollout/fidelity functions, parallelization/GPU support (particularly for bundles), etc.
# ### Package “settings” section (presuming top-level design choices are made to support global settings for solves, tolerances, backends, etc.; at present “PiccoloOptions” (and formerly “IpoptOptions”) are responsible for this sort of thing)
# ### Separate docs sections for each extension package (e.g. their docs on their Makie.jl extension do not cover plotting in detail (covered in a separate page on Bloch sphere visualizations) but rather detail the dependencies associated with the extension and the functions that are enabled by that extension being loaded) (notably theirs is incomplete; CUDA and Makie extensions are covered but GPUArrays and ChainRulesCore are not, since those extensions do not directly impact the user-facing API)

# ## Ideas from ControlSystems.jl

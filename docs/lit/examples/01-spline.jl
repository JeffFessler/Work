#=
# [Spline grid fitting](@id 01-spline)

This page illustrates finite-series modeling with splines
with the intention of comparing to INR.

[`URL`](https://github.com/JeffFessler/Work).

This page was generated from a single Julia file:
[01-spline.jl](@__REPO_ROOT_URL__/01-spline.jl).
=#

#md # In any such Julia documentation,
#md # you can access the source code
#md # using the "Edit on GitHub" link in the top right.

#md # The corresponding notebook can be viewed in
#md # [nbviewer](https://nbviewer.org/) here:
#md # [`01-spline.ipynb`](@__NBVIEWER_ROOT_URL__/01-spline.ipynb),
#md # and opened in [binder](https://mybinder.org/) here:
#md # [`01-spline.ipynb`](@__BINDER_ROOT_URL__/01-spline.ipynb).


# ### Setup

# Packages needed here.

using Template
using MIRTjim: jim, prompt
using InteractiveUtils: versioninfo
using Pltos


# The following line is helpful when running this file as a script;
# this way it will prompt user to hit a key after each figure is displayed.

isinteractive() ? jim(:prompt, true) : prompt(:draw);


# ### Setup

#=
=#


# ### Reproducibility

# This page was generated with the following version of Julia:

io = IOBuffer(); versioninfo(io); split(String(take!(io)), '\n')


# And with the following package versions

import Pkg; Pkg.status()

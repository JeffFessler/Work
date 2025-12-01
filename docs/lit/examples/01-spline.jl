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

using BSplineKit
using LinearAlgebra: I, norm
using MIRTjim: jim, prompt
using InteractiveUtils: versioninfo
using Plots
using Random: seed!
default(markerstrokecolor=:auto, label="")
seed!(0)


# The following line is helpful when running this file as a script;
# this way it will prompt user to hit a key after each figure is displayed.

isinteractive() ? jim(:prompt, true) : prompt(:draw);


# ### Data

N = 100
tmax = 3
t = tmax * sort(rand(N)) # random samples on [0,3]
f(t) = 1 < t < 2 # rect function
y = f.(t); # noiseless data for now

# plot data and true signal
yaxis = ((-0.2, 1.4), -0.2:0.2:1.4)
p0 = scatter(t, y; yaxis, color=:black, label="data")
plot!(p0, tf, f.(tf), label="true")


# B-spline settings
degree = 3
M = 40
knots = range(0, tmax, M) |> collect

# Make B-spline basis objects
basis = BSplineBasis(BSplineOrder(degree + 1), knots) # ; periodic=true

# Plot finely sampled B-spline basis functions
tf = range(-0.3, tmax+0.3, 1001)
Bf = hcat([b.(tf) for b in basis]...)
pb1 = plot(tf, Bf; title="Basis functions, degree=$degree")

# Sample B-spline basis objects at measurement points
B = hcat([b.(t) for b in basis]...)
pb2 = deepcopy(pb1)
scatter!(pb2, t, B; title="Basis function samples")

# LS fit of B-splines to data
xh = B \ y
yf = +([b.(tf) * xh[i] for (i,b) in enumerate(basis)]...)

p1 = deepcopy(p0)
plot!(p1, tf, yf; label="bspline")


#=
Evidently the unregularized B-spline fit here
is over-fitting, even though ``N < M``,
so next we try simple Tikhonov regularization.
=#
β = 2e-1
A = [B; sqrt(β)*I]

yz = [y; zeros(length(basis))]
xr = A \ yz
yfr = +([b.(tf) * xr[i] for (i,b) in enumerate(basis)]...)

p2 = deepcopy(p0)
plot!(p2, tf, yfr, label="bspline β=$β")


#=
But what regularization parameter β to choose?

Cross-validation is one approach to choose.

todo
=#



#todo include("../../../inc/reproduce.jl")

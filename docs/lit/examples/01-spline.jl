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

# plot data and true signal (finely sampled)
yaxis = ((-0.2, 1.4), -0.2:0.2:1.4)
tf = range(-0.3, tmax+0.3, 1001) # fine grid
ytf = f.(tf) # true function on a fine grid
p0 = scatter(t, y; yaxis, color=:black, label="data")
plot!(p0, tf, ytf, label="true")


# B-spline settings
degree = 3
M = 40
knots = range(0, tmax, M) |> collect

# Make B-spline basis objects
basis = BSplineBasis(BSplineOrder(degree + 1), knots) # ; periodic=true

# Plot finely sampled B-spline basis functions
Bf = hcat([b.(tf) for b in basis]...)
pb1 = plot(tf, Bf; title="Basis functions, degree=$degree")

# Sample B-spline basis objects at measurement points
B = hcat([b.(t) for b in basis]...)
pb2 = deepcopy(pb1)
scatter!(pb2, t, B; title="Basis function samples")

# LS fit of B-splines to data
xh = B \ y
yhf = +([b.(tf) * xh[i] for (i,b) in enumerate(basis)]...);

# Utility functions for NRMSE
do_nrmse(yr::AbstractVector) = norm(yr - ytf) / norm(ytf)
round3(x) = round(x, digits=3);

p1 = deepcopy(p0)
plot!(p1, tf, yhf; label="bspline, NRMSE=$(round3(do_nrmse(yhf)))")


#=
Evidently the unregularized B-spline fit here
is over-fitting, even though ``N < M``,
so next we try simple Tikhonov regularization.
=#
function do_fit(β::Real)
    A = [B; sqrt(β)*I]
    yz = [y; zeros(length(basis))]
    xr = A \ yz
    yr = +([b.(tf) * xr[i] for (i,b) in enumerate(basis)]...)
end
nrmse_fit(β::Real) = do_nrmse(do_fit(β));

β = 2e-1
yfr = do_fit(β);

p2 = deepcopy(p0)
plot!(p2, tf, yfr,
 label = "bspline β=$(round3(β)), NRMSE=$(round3(do_nrmse(yfr)))")


#=
But what regularization parameter β to choose?

Cross-validation is one approach to choose.

Here we use the simpler "oracle" approach
of finding the β value that leads
to the best fit to the true function
(which would be unknown in practice).

=#

#=
It turns out here that Tikhonov regularization
only slightly reduces NRMSE,
even in this very favorable case
where we use the ground-truth "oracle" to select β.

Apparently Tikhonov regularization is ineffective
for (uniformly spaced) B-spline fitting to non-uniformly data.
=#
βlist = 10 .^ (-5:0.2:1)
nrmse = nrmse_fit.(βlist)
βbest = argmin(nrmse_fit, βlist)
pn = plot(log10.(βlist), 100*nrmse; title="NRMSE (%)", xlabel="log10(β)",
 yaxis = ("NRMSE (%)", (10, 30)),
)
plot!(log10(βbest) * [1, 1], [10, 40], color=:red)

yfb = do_fit(βbest)
p3 = deepcopy(p0)
nrmse_best = nrmse_fit(βbest)
plot!(p3, tf, yfb,
 label="bspline β=$(round3(βbest)), NRMSE=$(round3(nrmse_best))")


include("../../inc/reproduce.jl")

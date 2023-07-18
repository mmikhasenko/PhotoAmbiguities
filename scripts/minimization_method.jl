cd(joinpath(@__DIR__, ".."))
import Pkg
Pkg.activate(".")
Pkg.instantiate()
# 
using PhotoAmbiguities
using Markdown
using LinearAlgebra
using ThreadsX

using Plots
import Plots.PlotMeasures:mm
theme(:wong2, frame=:box, grid=false, minorticks=true,
    guidefontvalign=:top, guidefonthalign=:right,
    xlim=(:auto, :auto), ylim=(:auto, :auto),
    lw=1.2, lab="", colorbar=false,
    bottom_margin=3mm, left_margin=3mm)

md"""
## Default model
"""

paperinput = "S0" => (0.499, 0),
    "D−1" => (0.201, 15.4),
    "D0" => (0.567, 174),
    "D1" => (0.624, -81.6)
# 
const wave2l = Dict('S'=>0, 'P'=>1, 'D'=>2, 'F'=>3)
# 
const waveset0 = [let
    Wave(v*cis(phi/180*π);
        ϵ = +1,
        l = wave2l[n[1]],
        m = eval(Meta.parse(n[2:end])))
end for (n,(v,phi)) in paperinput]

const mg = Model((waveset=waveset0, Pγ = 0.85))

# 
phsp = generate(10000);
w0 = intensity.(Ref(mg), eachrow(phsp))
data = phsp[w0 .> maximum(w0) .* rand(length(w0)), :]

md"""
## Local minima
"""


const data0 = data[1:1000,:]

@time mimimizarion_attempts = let Natt=10
    _mg = mg; #Model((waveset=mg.waveset, Pγ=0))
    ThreadsX.map(1:Natt) do _
        init = standardize(2rand(ComplexF64, 4).-(1+1im))
        go2min(_mg, data0, init)
    end
end



using Optim
"""
See also MIGrad in Chapter 4: Minuit Commands, https://root.cern.ch/download/minuit.pdf
""" 
function miunitstop(state)
    mt = state.metadata
    edm = mt["g(x)"]' * mt["~inv(H)"] * mt["g(x)"] / 2
    edm < 1e-3 * 0.1 * 1.0
end


init_test = standardize(2rand(ComplexF64, 4).-(1+1im))


@time let m=mg
    objective(x) = NNL(update(m, x |> unfold), data0)
    init = init_test |> fold
    opt_result = optimize(objective, init, BFGS())
    opt_result, opt_result.minimizer |> unfold 
end

@time let m=mg
    objective(x) = NNL(update(m, x |> unfold), data0)
    init = init_test |> fold
    opt_result = optimize(objective, init, BFGS(), Optim.Options(extended_trace=true, callback=miunitstop))
    opt_result, opt_result.minimizer |> unfold 
end

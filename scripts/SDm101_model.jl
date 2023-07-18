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

[a[1].minimum for a in mimimizarion_attempts]
[a[2] for a in mimimizarion_attempts]
[a[2] for a in mimimizarion_attempts] .|> standardize


selected_minima0 = [a[2] for a in mimimizarion_attempts] .|> standardize |> cluster |> drop_conjugate
sort!(selected_minima0, by=m->NNL(update(mg, m), data0))


md"""
# Experiment
"""

exp0 = let Ndata = 1000
    Experiment(mg, data, Ndata;
        initv = selected_minima0)
end

let 
    NNL0 = NNL(exp0, 0)
	plot(t->NNL(exp0, t)-NNL0, -0.5, nminima(exp0)-0.5)
	vline!(0:nminima(exp0)-1)
end
savefig(joinpath("plots","localminimas_SDm101.pdf"))


exps = let Nexp = 10, Ndata=200
    ThreadsX.map(1:Nexp) do _
        Experiment(mg, data, Ndata;
            initv = selected_minima0)
    end
end

let 
	plot()
    K = nminima(exp0)
	map(exps) do exp
		vt = range(-0.3, K-0.3, 200)
		NLL0 = NNL(exp,0)
		plot!(vt, NNL.(Ref(exp), (vt)) .- NLL0)
	end
    plot!(ylim=(-10,:auto))
    # 
    for i in 0:K-1
        lens!([-0.1,0.1] .+ i, [-5,20], inset=(1, bbox(0.05+0.95i/K,0.3,0.9/(K+1),0.2)))
        hline!(sp=i+2, [0.0], xaxis=nothing, yaxis=nothing, l=(:dash, :gray))
    end
	vline!(0:K-1, lc=2, ls=:dash)
end
savefig(joinpath("plots","bootstrap_SDm101.pdf"))

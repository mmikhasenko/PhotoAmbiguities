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
const wave2l = Dict('S'=>0, 'P'=>1, 'D'=>2)
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

function normalize(v)
    v .* cis(-arg(v[1])) ./ norm(v) 
end

const data0 = data[1:1000,:]

@time mimimizarion_attempts = let Natt=10
    _mg = mg; #Model((waveset=mg.waveset, Pγ=0))
    ThreadsX.map(1:Natt) do _
        init = normalize(2rand(ComplexF64, 4).-(1+1im))
        go2min(_mg, data0, init)
    end
end

[a[1].minimum for a in mimimizarion_attempts]
[a[2] for a in mimimizarion_attempts]
[a[2] for a in mimimizarion_attempts] .|> normalize

distance(a,b) = norm(a .- b)
function cluster(minima, ϵ=1e-5)
    selected = [minima[1]]
    for m in minima
        all(distance.(selected, Ref(m)) .> ϵ) && push!(selected, m)
    end
    return selected
end

function drop_conjugate(minima)
    return filter(minima) do m
        imag(m[2]) > 0
    end
end


selected_minima0 = [a[2] for a in mimimizarion_attempts] .|> normalize |> cluster |> drop_conjugate
sort!(selected_minima0, by=m->NNL(update(mg, m), data0))


md"""
# Experiment
"""

exp0 = let Ndata = 1000
    Experiment(mg, data, Ndata;
        initg = selected_minima[1],
        inita = selected_minima[2])
end

let 
    NNL0 = NNL(exp0, 0)
	plot(t->NNL(exp0, t)-NNL0, -0.5, 1.5)
	vline!([0,1])
    lens!([-0.1,0.1],      [-1,45], inset=(1, bbox(0.1,0.3,0.3,0.2)))
    lens!([-0.1,0.1] .+ 1, [-1,45], inset=(1, bbox(0.6,0.3,0.3,0.2)))
    hline!(sp=2, [0.0], xaxis=nothing, l=(:dash, :gray))
    hline!(sp=3, [0.0], xaxis=nothing, l=(:dash, :gray))
end
savefig(joinpath("plots","bootstrap_SDm101.pdf"))


exps = let Nexp = 10, Ndata=200
    ThreadsX.map(1:Nexp) do _
        Experiment(mg, data, Ndata;
            initg = selected_minima0[1],
            inita = selected_minima0[2])
    end
end

let 
	plot()
	map(exps) do exp
		vt = range(-0.5, 1.5, 200)
		NLL0 = NNL(exp,0)
		plot!(vt, NNL.(Ref(exp), (vt)) .- NLL0)
	end
	vline!([0,1])
end

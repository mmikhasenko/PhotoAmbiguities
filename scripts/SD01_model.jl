using PhotoAmbiguities
using Markdown
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

const waveset0 = [
	Wave(0.229; ϵ=+1, l=0, m=0),
	Wave(-0.217+0.31im; ϵ=+1, l=2, m=0),
	Wave(0.77+0.448im; ϵ=+1, l=2, m=1)
]
const mg = Model((waveset=waveset0,Pγ = 0.5))

md"""
Alternative model is defined using approximate minimum, calculated with Pγ=0.
"""
const ma = update(mg,
    [0.63,
    0.043+0.056im,
    0.28-0.713im])

md"""
## Generate data according to the general model
"""

phsp = generate(10000);
w0 = intensity.(Ref(mg), eachrow(phsp))
data = phsp[w0 .> maximum(w0) .* rand(length(w0)), :]

md"""
compare should be similar
"""
intensity(mg, data[1,:]), intensity(ma, data[1,:])


md"""
# Experiment
"""

exp0 = let Ndata = 1000
    Experiment(mg, data, Ndata;
        initg = getproperty.(mg.waveset, :value),
        inita = getproperty.(ma.waveset, :value))
end

let 
    NNL0 = NNL(exp0, 0)
	plot(t->NNL(exp0, t)-NNL0, -0.5, 1.5)
	vline!([0,1])
    lens!([-0.1,0.1],      [-1,20], inset=(1, bbox(0.1,0.3,0.3,0.2)))
    lens!([-0.1,0.1] .+ 1, [-1,20], inset=(1, bbox(0.6,0.3,0.3,0.2)))
    hline!(sp=2, [0.0], xaxis=nothing, l=(:dash, :gray))
    hline!(sp=3, [0.0], xaxis=nothing, l=(:dash, :gray))
end
savefig(joinpath("plots","localminimun.pdf"))
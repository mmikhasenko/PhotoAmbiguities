### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ c3d151e0-255f-11ee-0813-15d3d3fbf424
begin
	cd(joinpath(@__DIR__, ".."))
	import Pkg
	Pkg.activate(".")
	Pkg.instantiate()
	# 
	using PhotoAmbiguities
	using Markdown
	using LinearAlgebra
	using ThreadsX
	using Random
	# 
	using Plots
	import Plots.PlotMeasures:mm
end

# ╔═╡ ae60cdb8-e221-4921-a67e-bd725cbf2754
theme(:wong2, frame=:box, grid=false, minorticks=true,
	guidefontvalign=:top, guidefonthalign=:right,
	xlim=(:auto, :auto), ylim=(:auto, :auto),
	lw=1.2, lab="", colorbar=false,
	bottom_margin=3mm, left_margin=3mm)

# ╔═╡ 4beaa46a-78de-4f83-96a0-74d59ca91271
Random.seed!(4321)

# ╔═╡ b3ebefb0-dfb6-405e-8802-895e868005eb
md"""
## Default model
"""

# ╔═╡ 482e8bdb-3720-4fb6-9950-bf7700590183
const waveset0 = let
	# 
	paperinput = "S0" => (0.499, 0),
		"P0" => (-0.3, 80.1),
		"P1" => (-0.5, 190.1),
	    # "D−1" => (0.201, 15.4),
	    "D0" => (0.567, 174),
	    "D1" => (0.624, -81.6)
	# 
	wave2l = Dict('S'=>0, 'P'=>1, 'D'=>2, 'F'=>3)
	# 
	w0 = [Wave(v*cis(phi/180*π);
        ϵ = +1,
        l = wave2l[n[1]],
        m = eval(Meta.parse(n[2:end])))
		for (n,(v,phi)) in paperinput]
end

# ╔═╡ df66be2b-7608-4900-9f03-0dadd0df4c28
const mg = Model((; waveset=waveset0, Pγ = 0.185)) |> standardize

# ╔═╡ f174c8f7-ee0e-44f7-8b8f-098a5f523443
data = let
	phsp = generate(10000);
	w0 = intensity.(Ref(mg), eachrow(phsp))
	phsp[w0 .> maximum(w0) .* rand(length(w0)), :]
end

# ╔═╡ 56300d87-8323-483f-b3fc-4706e7728970
md"""
## Local minima
"""

# ╔═╡ e1a07689-fd52-4680-94ef-20c772b7a873
const data0 = data[1:1000,:];

# ╔═╡ a111d6f5-316a-4714-bf22-09af899eea93
@time mimimizarion_attempts = let Natt=30
	_mg = mg # Model((; waveset=mg.waveset, Pγ=0.0))
    ThreadsX.map(1:Natt) do _
        init = standardize(2rand(ComplexF64, length(_mg.waveset)).-(1+1im))
        go2min_precompute(_mg, data0, init)
    end
end;

# ╔═╡ c490aadf-e573-4507-bfc3-9591d460dc24
getindex.(mimimizarion_attempts, 2) .|> standardize

# ╔═╡ 2583396d-10f0-43ac-aee9-fae94a6d4a27
add_conjugate(minima) = vcat(minima, conj.(minima))

# ╔═╡ e8a80e40-7f31-49ab-98ac-6ea45f9628b6
selected_minima0 = getindex.(mimimizarion_attempts, 2) .|> 
	standardize |> add_conjugate |> cluster |> drop_conjugate

# ╔═╡ b1af62b7-0813-44be-aee7-2e23bf4a364a
sort!(selected_minima0, by=m->NNL(update(mg, m), data0))

# ╔═╡ c37b0ecb-9217-4957-9475-9c49eafccd07
map(selected_minima0) do m
	NNL(update(mg, m), data0)
end

# ╔═╡ a7aa4bd9-f87b-4e95-a7d6-40670f3389f0
md"""
## Pseudoexperiments
"""

# ╔═╡ e28de1ec-4039-4c13-b888-cdd4c224646b
exp0 = let Ndata = 1000
    Experiment(mg, data, Ndata;
        initv = selected_minima0)
end;

# ╔═╡ 4347ece5-a592-43ed-a115-8e0e0f255074
let 
    NNL0 = NNL(exp0, 0)
	plot(t->NNL(exp0, t)-NNL0, -0.5, nminima(exp0)-0.5)
	vline!(0:nminima(exp0)-1)
	plot!(xlab="across minima (t)", ylab="NNL")
end

# ╔═╡ d2d93852-909a-4f74-adf3-4a90a59efcdb
exps = let Nexp = 3, Ndata=200
    ThreadsX.map(1:Nexp) do _
        Experiment(mg, data, Ndata;
            initv = selected_minima0)
    end
end;

# ╔═╡ 986a472e-5140-457d-ab67-6e5b7275ea6d
let 
	plot()
    K = nminima(exps[1])
	map(exps) do exp
		vt = range(-0.3, K-0.3, 200)
		NLL0 = NNL(exp,0)
		plot!(vt, NNL.(Ref(exp), (vt)) .- NLL0)
	end
    plot!(ylim=(-10,:auto), ylab="NNL", xlab="across minima (t)")
    # 
    for i in 0:K-1
        lens!([-0.1,0.1] .+ i, [-5,20],
			inset=(1, bbox(0.05+0.95i/K,0.3,0.9/(K+1),0.2)))
        hline!(sp=i+2, [0.0], xaxis=nothing, yaxis=nothing, l=(:dash, :gray))
    end
	vline!(0:K-1, lc=2, ls=:dash)
end

# ╔═╡ Cell order:
# ╠═c3d151e0-255f-11ee-0813-15d3d3fbf424
# ╠═ae60cdb8-e221-4921-a67e-bd725cbf2754
# ╠═4beaa46a-78de-4f83-96a0-74d59ca91271
# ╟─b3ebefb0-dfb6-405e-8802-895e868005eb
# ╠═482e8bdb-3720-4fb6-9950-bf7700590183
# ╠═df66be2b-7608-4900-9f03-0dadd0df4c28
# ╠═f174c8f7-ee0e-44f7-8b8f-098a5f523443
# ╟─56300d87-8323-483f-b3fc-4706e7728970
# ╠═e1a07689-fd52-4680-94ef-20c772b7a873
# ╠═a111d6f5-316a-4714-bf22-09af899eea93
# ╠═c490aadf-e573-4507-bfc3-9591d460dc24
# ╠═2583396d-10f0-43ac-aee9-fae94a6d4a27
# ╠═e8a80e40-7f31-49ab-98ac-6ea45f9628b6
# ╠═b1af62b7-0813-44be-aee7-2e23bf4a364a
# ╠═c37b0ecb-9217-4957-9475-9c49eafccd07
# ╟─a7aa4bd9-f87b-4e95-a7d6-40670f3389f0
# ╠═e28de1ec-4039-4c13-b888-cdd4c224646b
# ╠═4347ece5-a592-43ed-a115-8e0e0f255074
# ╠═d2d93852-909a-4f74-adf3-4a90a59efcdb
# ╠═986a472e-5140-457d-ab67-6e5b7275ea6d

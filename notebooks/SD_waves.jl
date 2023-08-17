### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ c3d151e0-255f-11ee-0813-15d3d3fbf424
# ╠═╡ show_logs = false
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
	using DataFrames
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

# ╔═╡ b3ebefb0-dfb6-405e-8802-895e868005eb
md"""
## Default model
"""

# ╔═╡ 482e8bdb-3720-4fb6-9950-bf7700590183
const waveset0 = let
	# 
	paperinput = "S0" => (0.499, 0),
	    "D−1" => (0.201, 15.4),
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
end;

# ╔═╡ df66be2b-7608-4900-9f03-0dadd0df4c28
const mg = Model((; waveset=waveset0, Pγ = 0.85)) |> standardize

# ╔═╡ f174c8f7-ee0e-44f7-8b8f-098a5f523443
data = let
	phsp = generate(10000);
	w0 = intensity.(Ref(mg), eachrow(phsp))
	phsp[w0 .> maximum(w0) .* rand(length(w0)), :]
end;

# ╔═╡ 56300d87-8323-483f-b3fc-4706e7728970
md"""
## Local minima
"""

# ╔═╡ e1a07689-fd52-4680-94ef-20c772b7a873
const data0 = data[1:1000,:];

# ╔═╡ a111d6f5-316a-4714-bf22-09af899eea93
@time mimimizarion_attempts = let Natt=30
    _mg = mg; #Model((waveset=mg.waveset, Pγ=0))
    ThreadsX.map(1:Natt) do _
        init = standardize(2rand(ComplexF64, 4).-(1+1im))
        go2min_precompute(_mg, data0, init)
    end
end;

# ╔═╡ 93e42e31-d99f-4653-8d0d-00e1f2422643
function fraction2string(a,b; digits=2)
	f = round(100*(1-length(b)/length(a)); digits)
	"$(f)%"
end

# ╔═╡ e8a80e40-7f31-49ab-98ac-6ea45f9628b6
begin
	standartized0 = getindex.(mimimizarion_attempts, 2) .|> standardize;
	clustered0 = standartized0 |> cluster |> drop_conjugate;
	# 
	@info "Clustering removes $(fraction2string(standartized0,clustered0))"
	# 
	selected0 = clustered0 |> drop_conjugate;
	@info "Conjugation filtering removes $(fraction2string(clustered0,selected0))"
end

# ╔═╡ b1af62b7-0813-44be-aee7-2e23bf4a364a
sort!(selected0, by=m->NNL(update(mg, m), data0))

# ╔═╡ 3f31065a-41ac-4b88-8715-30403cebe9f6
map(selected0) do m
	NNL(update(mg, m), data0)
end

# ╔═╡ a7aa4bd9-f87b-4e95-a7d6-40670f3389f0
md"""
## Pseudoexperiments
"""

# ╔═╡ e28de1ec-4039-4c13-b888-cdd4c224646b
exp0 = let Ndata = 1000
    Experiment(mg, data, Ndata;
        initv = selected0)
end;

# ╔═╡ 9e9e17ca-943b-478f-bf59-805c31fc7029
exp0.pv[1], exp0.pv[2]

# ╔═╡ 4347ece5-a592-43ed-a115-8e0e0f255074
let 
    NNL0 = NNL(exp0, 0)
	plot(t->NNL(exp0, t)-NNL0, -0.5, nminima(exp0)-0.5)
	vline!(0:nminima(exp0)-1)
	plot!(xlab="across minima (t)", ylab="NNL")
end

# ╔═╡ d2d93852-909a-4f74-adf3-4a90a59efcdb
exps = let Nexp = 300, Ndata=100
    ThreadsX.map(1:Nexp) do _
        Experiment(mg, data, Ndata;
            initv = selected0)
    end
end;

# ╔═╡ 986a472e-5140-457d-ab67-6e5b7275ea6d
let 
	plot()
    K = nminima(exps[1])
	map(exps[rand(1:length(exps),10)]) do exp
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

# ╔═╡ fa3c69b7-d4f6-4847-bd3f-7e29374dc507
md"""
## How different are the minima?
The scalar product of the amplitude values
shows how different the wave values in the minimas are.
 - The product of the solution to itselt is 1, $(i|i) = 1$ due to normalization.
 - The production of minimization result to the starting point is close to one if the starting guess is given resonably, $(i|i_0) \approx 1$ (**orange** and **blue** points)
 - The product of solutions in two minimas should not be close to 1 if solutions are different, $(i|j) \neq 1$ (**black** points)
"""

# ╔═╡ ea1dd9b5-24d6-4bca-aea2-06fdaa92ba7e
let
	plot()
	#
	scatter!(
		map(exps) do e
			dot(e.pv[1], e.pv[2])
		end, lab="(1|2)")
	# 
	scatter!(
			map(exps) do e
			dot(e.pv[1], exp0.pv[1])
		end, lab="(1|o1)")
	# 
	scatter!(
			map(exps) do e
			dot(e.pv[2], exp0.pv[2])
		end, lab="(2|o2)")
	# 
	scatter!(
		map(exps) do e
			dot(e.pv[1], exp0.pv[2])
		end, lab="(1|o2)")
	# 
	plot!(xlab="Re(i|j)", ylab="Im(i|j)")
end

# ╔═╡ abff57af-bb32-4a49-bbaa-f2226cdeaabc
md"""
## Can "local" gets deeper that global?
"""

# ╔═╡ df4f4f64-075b-48d6-a951-55f91ade85c1
begin
	exps′ = filter(exps) do exp
		norm(exp.pv[1]-exp.pv[2]) > 0.1
	end
	close_fraction = 1-length(exps′)/length(exps)
	@info "The proximity cut removes $(round(100*close_fraction, digits=1))%"
	# 
	df_NLL = map(exps′) do e
		dmin = norm(e.pv[1]-e.pv[2])
		dprod = dot(e.pv[1],e.pv[2])
		dNLL = NNL(e,1)-NNL(e,0)
		(; dmin, dprod, dNLL)
	end |> DataFrame
	plot(layout=grid(1,2), size=(800,350))
	stephist!(sp=1, df_NLL.dmin, bins=20,
		title="distance between minima", ylab="#counts")
	f = round(count(df_NLL.dNLL .< 0) / length(df_NLL.dNLL); digits=2)
	stephist!(sp=2, df_NLL.dNLL, bins=25, lab="f_{<0} = $f",
		title="difference of NLL", xlab="local-global", ylab="#counts")
	# 
		
	f_swapped = round(count(df_NLL.dNLL .< 0) / length(df_NLL.dNLL), digits=3)
	@info "The chance swapping minima is $(fraction2string(df_NLL.dNLL, df_NLL.dNLL[df_NLL.dNLL .>= 0]))%"
	# 
	plot!()
end

# ╔═╡ Cell order:
# ╠═c3d151e0-255f-11ee-0813-15d3d3fbf424
# ╠═ae60cdb8-e221-4921-a67e-bd725cbf2754
# ╟─b3ebefb0-dfb6-405e-8802-895e868005eb
# ╠═482e8bdb-3720-4fb6-9950-bf7700590183
# ╠═df66be2b-7608-4900-9f03-0dadd0df4c28
# ╠═f174c8f7-ee0e-44f7-8b8f-098a5f523443
# ╟─56300d87-8323-483f-b3fc-4706e7728970
# ╠═e1a07689-fd52-4680-94ef-20c772b7a873
# ╠═a111d6f5-316a-4714-bf22-09af899eea93
# ╠═93e42e31-d99f-4653-8d0d-00e1f2422643
# ╠═e8a80e40-7f31-49ab-98ac-6ea45f9628b6
# ╠═b1af62b7-0813-44be-aee7-2e23bf4a364a
# ╠═3f31065a-41ac-4b88-8715-30403cebe9f6
# ╟─a7aa4bd9-f87b-4e95-a7d6-40670f3389f0
# ╠═9e9e17ca-943b-478f-bf59-805c31fc7029
# ╠═e28de1ec-4039-4c13-b888-cdd4c224646b
# ╠═4347ece5-a592-43ed-a115-8e0e0f255074
# ╠═d2d93852-909a-4f74-adf3-4a90a59efcdb
# ╠═986a472e-5140-457d-ab67-6e5b7275ea6d
# ╟─fa3c69b7-d4f6-4847-bd3f-7e29374dc507
# ╠═ea1dd9b5-24d6-4bca-aea2-06fdaa92ba7e
# ╟─abff57af-bb32-4a49-bbaa-f2226cdeaabc
# ╠═df4f4f64-075b-48d6-a951-55f91ade85c1

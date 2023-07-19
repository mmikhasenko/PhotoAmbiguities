using PhotoAmbiguities
using Test
using Markdown

md"""
## To speed up the calculations,
the angular functions are precomputed.
It is done with `ForwardDiff.hessian`.
The test demonstrate that the output using the precomputed matrices is the same as the original function.
"""
let
    waveset0 = [
        Wave(0.229; ϵ=+1, l=0, m=0),
        Wave(-0.217+0.31im; ϵ=+1, l=2, m=0),
        Wave(0.77+0.448im; ϵ=+1, l=2, m=1)
    ]
    mg = Model((waveset=waveset0, Pγ = 0.85))
    # 
    phsp = generate(10000);
    w0 = intensity.(Ref(mg), eachrow(phsp))
    data = phsp[w0 .> maximum(w0) .* rand(length(w0)), :]
    # 
    Hv = pullbilinear(mg, data)
    @test NNL(mg, data) == NNL(getproperty.(mg.waveset, :value) |> fold, Hv)
end
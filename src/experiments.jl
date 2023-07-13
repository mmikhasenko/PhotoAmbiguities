function NNL(m::Model, data)
	@unpack waveset = m
	I(τ) = intensity(m, τ)
	Nd = size(data,1)
	return -sum(log ∘ I, eachrow(data)) + Nd * sum(w->abs2(w.value), waveset)
end

function unfold(v::AbstractVector{<:Real})
    N = div(length(v)+1,2)
    return v[1:N] + 1im * [0, v[N+1:end]...]
end
function fold(v::AbstractVector{<:Number})
    [real(v)..., imag(v[2:end])...]
end

function go2min(m::Model, data, startvalues)
	objective(x) = NNL(update(m, x |> unfold), data)
	init = startvalues |> fold
	opt_result = optimize(objective, init, BFGS())
	return opt_result, opt_result.minimizer |> unfold 
end



struct Experiment{M,T}
    model::M
    data::T
    pg::Vector{<:Number}
    pa::Vector{<:Number}
end

function Experiment(m::Model, data, N::Int; initg::Vector{<:Number}, inita::Vector{<:Number})
    pickeditems = rand(1:size(data,1), N)
    _data = data[pickeditems,:]
    # 
    pg = go2min(m, _data, initg)[2]
    pa = go2min(m, _data, inita)[2]
    # 
    Experiment(m, _data, pg, pa)
end

track(p0,pa,t) = (p0+pa)/2 + cis(π*t) * (p0-pa)/2

function NNL(exp::Experiment, t)
    @unpack pg, pa = exp
    f(t) = track(pg,pa,t)
    return NNL(update(exp.model, f(t)), exp.data)
end

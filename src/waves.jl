
struct Wave{T}
	ϵ::Int
	l::Int
	m::Int
	value::T
end
Wave(value; l::Int, m::Int, ϵ::Int) = Wave(ϵ, l, m, value)
Wave(value, w::Wave) = Wave(w.ϵ, w.l, w.m, value)
update(waves::AbstractVector{Wave}, values::AbstractVector{<:Number}) = Wave.(values,waves)

const Model = NamedTuple{(:waveset,:Pγ), Tuple{Vector{<:Wave}, Float64}}
update(m::Model, values::AbstractVector{<:Number}) = Model((;waveset=Wave.(values,m.waveset), m.Pγ))

function standardize(v::AbstractArray)
    v .* cis(-arg(v[1])) ./ norm(v) 
end
standardize(m::Model) = update(m, getproperty.(m.waveset, :value) |> standardize)

function intensity(model::Model, τ)
	@unpack cosθ,ϕ,Φ = τ
	@unpack waveset, Pγ = model
	# 
	vmϵ_all = map(waveset) do w
		f = sqrt(2w.l+1)*d(w.l, w.m, cosθ) * w.value
		(f,w.m,w.ϵ)
	end
	# 
	i0 = sum(Iterators.product(vmϵ_all,vmϵ_all)) do ((f,m,ϵ), (f′,m′,ϵ′))
		f*conj(f′)*cos((m-m′)*ϕ) * (ϵ == ϵ′) |> real
	end
	i1 = -sum(Iterators.product(vmϵ_all,vmϵ_all)) do ((f,m,ϵ), (f′,m′,ϵ′))
		ϵ*f*conj(f′)*cos((m+m′)*ϕ) * (ϵ == ϵ′) |> real
	end
	i2 = -sum(Iterators.product(vmϵ_all,vmϵ_all)) do ((f,m,ϵ), (f′,m′,ϵ′))
		ϵ*f*conj(f′)*sin((m+m′)*ϕ) * (ϵ == ϵ′) |> real
	end
	return i0 - Pγ*(i1*cos(Φ) + i2*sin(Φ))
end


function ε(l::Int,m::Int,u)
	N = factorial(l)*sqrt(factorial(l-m)*factorial(l+m))
	X = sum(max(0,-m):min(l,l-m)) do k
		D = factorial(l-m-k)*factorial(l-k)*factorial(m+k)*factorial(k)
		(-1)^k*u^(2k+m-l) / D
	end
	return X*N
end

function d(l::Int,m::Int,cosθ)
	u = tan(acos(cosθ)/2)
	return (u/(1+u^2))^l*(-1)^m*ε(l,m,u)
end

Y(l::Int,m::Int,cosθ) = sqrt((2l+1)/(4π))*d(l,m,cosθ)
Y(l::Int,m::Int,cosθ,ϕ) = sqrt((2l+1)/(4π))*d(l,m,cosθ)*cis(m*ϕ)

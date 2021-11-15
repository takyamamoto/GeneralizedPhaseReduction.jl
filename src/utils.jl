# identity matrix
eye(T::Type, m, n) = Matrix{T}(LinearAlgebra.I, m, n)
eye(n) = eye(Float64, n, n)

# vectorization fucntion
vec(X) = vcat(X...)

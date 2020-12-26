# Define MyType
struct MyType  # Don't inherit from anything
    value
end
value(x::MyType) = x.value

# Define its methods
Base.one(::Type{MyType}) = MyType(1)
Base.zero(::Type{MyType}) = MyType(0)

Base.:+(x::MyType, y::MyType) = MyType(value(x)+value(y))
Base.:*(x::MyType, y::MyType) = MyType(value(x)*value(y))
Base.:-(x::MyType, y::MyType) = MyType(value(x)-value(y))
Base.:/(x::MyType, y::MyType) = MyType(value(x)/value(y))

Base.adjoint(x::MyType) = MyType(adjoint(value(x)))

function Base.isapprox(x::MyType, y::MyType; atol, rtol)
    return isapprox(value(x), value(y), atol=atol, rtol=rtol)
end
function Base.rtoldefault(::Type{MyType}, ::Type{MyType}, z)
    # I don't really want to be restricted to Float64 tolerance.
    return Base.rtoldefault(Float64, Float64, z)
end

# Shouldn't have to define these
Base.one(::MyType) = one(MyType)
Base.zero(::MyType) = zero(MyType)

Base.:+(x::MyType, y) = MyType(value(x)+y)
Base.:*(x::MyType, y) = MyType(value(x)*y)

Base.abs(x::MyType) = MyType(abs(value(x)))
Base.:<(x::MyType, y::MyType) = value(x) < value(y)
Base.inv(x::MyType) = MyType(inv(value(x)))

Base.promote_rule(::Type{Any}, ::Type{MyType}) = MyType
Base.convert(::Type{MyType}, x::MyType) = x
Base.convert(::Type{MyType}, x) = MyType(x)

Base.iterate(x::MyType) = (value(x), nothing)
Base.iterate(::MyType, ::Any) = nothing
Base.length(x::MyType) = 1

# Begin check for ((X-I)')^-1 ≈ Y
X = [
    0.4 1;
    -0.2 0.8;
]
X = map(x -> MyType(x), X)

Y = [
    -0.625 0.625;
    -3.125 -1.875;
]
Y = map(x -> MyType(x), Y)

using LinearAlgebra
println((X-I)')
println(inv((X-I)'))
println(inv((X-I)') ≈ Y)


dataframe = [
    1 2;
    1 3;
    2 3;
    2 4;
    2 5;
    3 2;
    4 2;
    5 3;
]

max_trial = maximum(dataframe[:, 1])
max_c1 = Array{Any}(undef, max_trial)
for c in 1:max_trial
    c1 = findall(dataframe[:, 1] .== c)
    max_c1[c] = maximum(dataframe[c1, 2])
end

max_c1 = [maximum(dataframe[c1, 2]) for c in unique()]

# Name of trials
naming = unique(dataframe[!,:9])
dataframe = dataframe[!,:1:8]
# Number of trials
max_trial = maximum(dataframe[!,:1])
# Number of variables on each trial
max_c1 = []
for c in 1:max_trial
    c1 = findall(dataframe[:,1] .== c)
    max_c1[c] = maximum(dataframe[c1,2])
end

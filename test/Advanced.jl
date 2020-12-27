using DiscreteMarkovChains
using Test
using LinearAlgebra

@testset "Custom eltype" begin
    struct MyType
        value
    end
    value(x::MyType) = x.value

    Base.one(::Type{MyType}) = MyType(1)
    Base.one(::MyType) = one(MyType)
    Base.zero(::Type{MyType}) = MyType(0)
    Base.zero(::MyType) = zero(MyType)

    Base.:+(x::MyType, y::MyType) = MyType(value(x)+value(y))
    Base.:*(x::MyType, y::MyType) = MyType(value(x)*value(y))
    Base.:-(x::MyType, y::MyType) = MyType(value(x)-value(y))
    Base.:/(x::MyType, y::MyType) = MyType(value(x)/value(y))
    Base.:+(x::MyType, y) = MyType(value(x)+y)
    Base.:*(x::MyType, y) = MyType(value(x)*y)
    Base.:-(x::MyType, y) = MyType(value(x)-y)
    Base.:/(x::MyType, y) = MyType(value(x)/y)

    Base.adjoint(x::MyType) = MyType(adjoint(value(x)))
    Base.abs(x::MyType) = MyType(abs(value(x)))
    Base.:<(x::MyType, y::MyType) = value(x) < value(y)
    Base.:>(x::MyType, y::MyType) = value(x) > value(y)
    Base.inv(x::MyType) = MyType(inv(value(x)))

    function Base.isapprox(x::MyType, y::MyType; atol, rtol)
        return isapprox(value(x), value(y), atol=atol, rtol=rtol)
    end
    function Base.isapprox(x::MyType, y::MyType)
        return isapprox(value(x), value(y))
    end
    function Base.rtoldefault(::Type{MyType}, ::Type{MyType}, z)
        return Base.rtoldefault(Float16, Float16, z)
    end
    Base.promote_rule(::Type{Any}, ::Type{MyType}) = MyType
    Base.convert(::Type{MyType}, x::MyType) = x
    Base.convert(::Type{MyType}, x) = MyType(x)
    Base.iterate(x::MyType) = x

    T = first.([
        ([0.5 0; 0 0.2],) ([0.5 0; 0 0.8],);
        ([0.7 0; 0 0.2],) ([0.3 0; 0 0.8],);
    ])

    T = [
        MyType(0.5) MyType(0.5);
        MyType(0.2) MyType(0.8);
    ]

    X = DiscreteMarkovChain(T)

    @test communication_classes(X) == ([[1, 2]], Any[true])
    @test is_ergodic(X)
    s = stationary_distribution(X)
    for (i, v) in enumerate(s)
        @test v ≈ [MyType(2/7), MyType(5/7)][i]
    end
    #@test stationary_distribution(X) ≈ [MyType(2/7), MyType(5/7)]

    T = first.([
        ([1 0; 0 1],) ([0 0; 0 0],);
        ([0 0; 0 0],) ([1 0; 0 1],);
    ])

    X = DiscreteMarkovChain(T)

    @test communication_classes(X) == ([[1], [2]], Any[true, true])
    @test !is_ergodic(X)
    # @test stationary_distribution(X) == T

end

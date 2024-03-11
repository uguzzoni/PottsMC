using SafeTestsets

@safetestset "model" begin include("model.jl") end
@safetestset "anneal" begin include("annealing.jl") end
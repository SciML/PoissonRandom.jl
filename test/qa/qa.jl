using SafeTestsets

@safetestset "Aqua" begin
    include("qa_aqua.jl")
end

@safetestset "JET static analysis" begin
    include("qa_jet.jl")
end

@safetestset "ExplicitImports" begin
    include("qa_explicitimports.jl")
end

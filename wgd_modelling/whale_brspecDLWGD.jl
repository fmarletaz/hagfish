using Whale, Turing, Random, NewickTree, Distributions, DataFrames, CSV
using Plots, StatsPlots
Random.seed!(562);

default(grid=false, size=(500, 800), titlefontsize=9, title_loc=:left, guidefont=8)

# Species tree
tree = readnw("((Sackow:5.65,Strpur:5.65):0.25,(Bralan:5.82,((Parata:4.58,Petmar:4.58):0.59,(Calmil:4.57,(Lepocu:4.17,(Homsap:3.19,Galgal:3.19):0.98):0.40):0.60):0.65):0.08);")

# Insert hypothethized WGD nodes
insertnode!(getlca(tree, "Petmar", "Petmar"), dist=2.29, name="wgd_lamp")
insertnode!(getlca(tree, "Parata", "Parata"), dist=2.29, name="wgd_hag")
insertnode!(getlca(tree, "Petmar", "Parata"), dist=0.295, name="wgd_cyclo")
insertnode!(getlca(tree, "Galgal", "Calmil"), dist=0.30, name="wgd_2r")
insertnode!(getlca(tree, "Parata", "Calmil"), dist=0.325, name="wgd_1r")


# Definition of the rates model
n = length(postwalk(tree))
rates = DLWGD(λ=zeros(n), μ=zeros(n), q=[0.2, 0.2, 0.2, 0.2, 0.2], η=0.9)

# Definition of the Whale model
model = WhaleModel(rates, tree, 0.1)

# Get the data
data  = read_ale("data_all/", model)

# Define the Bayesian probabilistic model using Turing.jl
@model branchrates(M, n, X, ::Type{T}=Float64) where T = begin
    η ~ Beta(3, 1)
    λ̄ ~ Normal(log(0.15), 2)
    μ̄ ~ Normal(log(0.15), 2)
    τ ~ Exponential(0.1)
    λ ~ MvNormal(fill(λ̄, n-1), τ)
    μ ~ MvNormal(fill(μ̄, n-1), τ)
    q1 ~ Beta()
    q2 ~ Beta()
    q3 ~ Beta()
    q4 ~ Beta()
    q5 ~ Beta()
    X ~ model((λ=λ, μ=μ, η=η, q=[q1, q2, q3, q4, q5]))
end

# Get an MCMC sample using NUTS
chain = sample(branchrates(model, n, data), NUTS(), 1000)

# Compute bayes factors
summarize(chain[[:q1, :q2, :q3, :q4, :q5]], Whale.bayesfactor)

# Plot posterior dist for retention params
plotd = plot(chain[[:q1, :q2, :q3, :q4, :q5]], size=(700, 900));
savefig(plotd, "brspec_DLWGD_params_posterior_8931trees.svg")

# Save posterior distributions as csv
posterior = DataFrame(chain)
CSV.write("brspec_DLWGD_params_posterior_8931trees.csv", posterior)

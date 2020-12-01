var documenterSearchIndex = {"docs":
[{"location":"references/#References","page":"References","title":"References","text":"","category":"section"},{"location":"references/","page":"References","title":"References","text":"This package is based on:","category":"page"},{"location":"references/","page":"References","title":"References","text":"McElreath: Statistical Rethinking 2nd edition","category":"page"},{"location":"references/","page":"References","title":"References","text":"There is no shortage of additional good books on Bayesian statistics. A few of my favorites are:","category":"page"},{"location":"references/","page":"References","title":"References","text":"Bolstad: Introduction to Bayesian statistics\nBolstad: Understanding Computational Bayesian Statistics\nGelman, Hill: Data Analysis using regression and multileve,/hierachical models\nKruschke: Doing Bayesian Data Analysis\nLee, Wagenmakers: Bayesian Cognitive Modeling\nGelman, Carlin, and others: Bayesian Data Analysis\nCausal Inference in Statistics - A Primer\nBetancourt: A Conceptual Introduction to Hamiltonian Monte Carlo","category":"page"},{"location":"acknowledgements/#Acknowledgements","page":"Acknowledgements","title":"Acknowledgements","text":"","category":"section"},{"location":"acknowledgements/","page":"Acknowledgements","title":"Acknowledgements","text":"Of course, without this excellent textbook by Richard McElreath, this package would not have been possible. The author has also been supportive of this work and gave permission to use the datasets.","category":"page"},{"location":"acknowledgements/","page":"Acknowledgements","title":"Acknowledgements","text":"Richard Torkar has taken the lead in developing the Turing versions of the models in chapter 8 and subsequent chapters. ","category":"page"},{"location":"acknowledgements/","page":"Acknowledgements","title":"Acknowledgements","text":"Tamas Papp has been very helpful during the development of the DynamicHMC versions of the models.","category":"page"},{"location":"acknowledgements/","page":"Acknowledgements","title":"Acknowledgements","text":"The TuringLang team and #turing contributors on Slack have been extremely helpful! The Turing examples by Cameron Pfiffer are followed closely in several example scripts.","category":"page"},{"location":"acknowledgements/","page":"Acknowledgements","title":"Acknowledgements","text":"The increasing use of Particles to represent quap approximations is possible thanks to the package MonteCarloMeasurements.jl. Soss.jl and related write-ups introduced me to that option.","category":"page"},{"location":"acknowledgements/","page":"Acknowledgements","title":"Acknowledgements","text":"Developing rethinking must have been an on-going process over several years, StatisticalRethinking.jl and associated packages will likely follow a similar path.","category":"page"},{"location":"","page":"Functions","title":"Functions","text":"CurrentModule = StatisticalRethinking","category":"page"},{"location":"#sr_path","page":"Functions","title":"sr_path","text":"","category":"section"},{"location":"","page":"Functions","title":"Functions","text":"sr_path(parts...)","category":"page"},{"location":"#StatisticalRethinking.sr_path-Tuple","page":"Functions","title":"StatisticalRethinking.sr_path","text":"sr_path\n\nRelative path using the StatisticalRethinking src/ directory.\n\nExample to get access to the data subdirectory\n\nsr_path(\"..\", \"data\")\n\nNote that in the projects, e.g. StatisticalRethinkingStan.jl and StatisticalRethinkingTuring.jl, the DrWatson approach is a better choics, i.e: sr_datadir(filename)\n\n\n\n\n\n","category":"method"},{"location":"#sr_datadir","page":"Functions","title":"sr_datadir","text":"","category":"section"},{"location":"","page":"Functions","title":"Functions","text":"sr_datadir(parts...)","category":"page"},{"location":"#StatisticalRethinking.sr_datadir-Tuple","page":"Functions","title":"StatisticalRethinking.sr_datadir","text":"sr_datadir\n\nRelative path using the StatisticalRethinking src/ directory.\n\nExample to access Howell1.csv in StatisticalRethinking:\n\ndf = CSV.read(sr_datadir(\"Howell1.csv\"), DataFrame)\n\n\n\n\n\n","category":"method"},{"location":"#link","page":"Functions","title":"link","text":"","category":"section"},{"location":"","page":"Functions","title":"Functions","text":"link(dfa::DataFrame, vars, xrange)","category":"page"},{"location":"#StatisticalRethinking.link-Tuple{DataFrame,Any,Any}","page":"Functions","title":"StatisticalRethinking.link","text":"link\n\nCompute the link function for standardized variables.\n\nlink(dfa, vars, xrange)\n\n\nExtended help\n\nRequired arguments\n\n* `df::DataFrame`                      : Chain samples converted to a DataFrame\n* `vars::Vector{Symbol}`               : Variables in DataFrame (2 variables)\n* `xrange::range`                      : Range over which link values are computed\n\nOptional arguments\n\n* `xbar::Float64`                      : Mean value of observed predictor\n* `ybar::Float64`                      : Mean value of observed outcome (requires xbar argument)\n\nReturn values\n\n* `result`                             : Vector of link values\n\n\n\n\n\n","category":"method"},{"location":"#rescale","page":"Functions","title":"rescale","text":"","category":"section"},{"location":"","page":"Functions","title":"Functions","text":"rescale(x::Vector{Float64}, xbar::Float64, xstd::Float64)","category":"page"},{"location":"#StatisticalRethinking.rescale-Tuple{Array{Float64,1},Float64,Float64}","page":"Functions","title":"StatisticalRethinking.rescale","text":"rescale\n\nRescale a vector to \"un-standardize\", the opposite of scale!().\n\nrescale(x, xbar, xstd)\n\n\nExtended help\n\nRequired arguments\n\n* `x::Vector{Float64}`                 : Vector to be rescaled\n* `xbar`                               : Mean value for rescaling\n* `xstd`                               : Std for rescaling\n\nReturn values\n\n* `result::AbstractVector`             : Rescaled vector\n\n\n\n\n\n","category":"method"},{"location":"#sample","page":"Functions","title":"sample","text":"","category":"section"},{"location":"","page":"Functions","title":"Functions","text":"sample(df::DataFrame, n; replace=true, ordered=false)","category":"page"},{"location":"#StatsBase.sample-Tuple{DataFrame,Any}","page":"Functions","title":"StatsBase.sample","text":"sample\n\nSample rows from a DataFrame\n\nMethod\n\nsample(df, n; replace, ordered) \n\nRequired arguments\n\n* `df::DataFrame`                      : DataFrame\n* `n::Int`                             : Number of samples\n\nOptional argument\n\n* `rng::AbstractRNG`                   : Random number generator\n* `replace::Bool=true`                 : Sample with replace \n* `ordered::Bool=false`                : Sort sample \n\nReturn values\n\n* `result`                             : Array of samples\n\n\n\n\n\n","category":"method"},{"location":"#hpdi","page":"Functions","title":"hpdi","text":"","category":"section"},{"location":"","page":"Functions","title":"Functions","text":"hpdi(x::Vector{T}; alpha::Real=0.05) where {T<:Real}","category":"page"},{"location":"#StatisticalRethinking.hpdi-Union{Tuple{Array{T,1}}, Tuple{T}} where T<:Real","page":"Functions","title":"StatisticalRethinking.hpdi","text":"hpdi\n\nCompute high density region.\n\nhpdi(x; alpha)\n\n\nDerived from hpd in MCMCChains.jl.\n\nBy default alpha=0.11 for a 2-sided tail area of p < 0.055% and p > 0.945%.\n\n\n\n\n\n","category":"method"},{"location":"#pairsplot","page":"Functions","title":"pairsplot","text":"","category":"section"},{"location":"","page":"Functions","title":"Functions","text":"pairsplot(df::DataFrame, vars::Vector{Symbol}, fig::AbstractString)","category":"page"},{"location":"#StatisticalRethinking.pairsplot-Tuple{DataFrame,Array{Symbol,1},AbstractString}","page":"Functions","title":"StatisticalRethinking.pairsplot","text":"pairsplot\n\nA simple version of StatsPlots' cornerplot.\n\npairsplot(df, vars, fig)\n\n\nRequired arguments\n\n* `df::DataFrame`                      : DataFrame containing the variables (as columns)\n* `vars::Vector{Symbol}`               : Vector of variables to include in plot\n* `fig::AbstractString`                : File to store the produced plot\n\n\n\n\n\n","category":"method"},{"location":"#plotbounds","page":"Functions","title":"plotbounds","text":"","category":"section"},{"location":"","page":"Functions","title":"Functions","text":"plotbounds(\n  df::DataFrame, \n  xvar::Symbol,\n  yvar::Symbol, \n  dfs::DataFrame, \n  linkvars::Vector{Symbol};\n  fig::AbstractString=\"\",\n  bounds::Vector{Symbol}=[:range, :hpdi],\n  title::AbstractString=\"\",\n  xlab::AbstractString=String(xvar),\n  ylab::AbstractString=String(yvar),\n  alpha::Float64=0.11,\n  colors::Vector{Symbol}=[:lightgrey, :grey],\n  stepsize::Float64=0.01\n)","category":"page"},{"location":"#StatisticalRethinking.plotbounds-Tuple{DataFrame,Symbol,Symbol,DataFrame,Array{Symbol,1}}","page":"Functions","title":"StatisticalRethinking.plotbounds","text":"plotbounds\n\nPlot regression line and intervals based on Stan samples of coeffficients.\n\nplotbounds(df, xvar, yvar, dfs, linkvars; fig, bounds, title, xlab, ylab, alpha, colors, stepsize, rescale_axis)\n\n\nRequired arguments\n\n* `df::DataFrame`                      : DataFrame with observed variables and scaled variables\n* `xvar::Symbol`                       : X variable in df\n* `yvar::Symbol`                       : Y variable in df\n* `dfs::DataFrame`                     : DataFrame with Stan samples\n* `linkvars::Vector{Symbol}`           : Initial 2 Symbols are regression coefficients,\n                                         3rd - only for :predict bounds indicates σ.\n\nOptional arguments\n\n* `fig::AbstractString=\"\"`             : File to store the plot. If \"\", a plot is returned\n* `bounds::Vector{Symbol}`             : Bounds to display, see below\n* `title::AbstractString=\"\"`           : Title for plot\n* `title::AbstractString=String(xvar)` : X axis variable label\n* `title::AbstractString=String(yvar)` : Y axis variable label\n* `alpha::Float64=0.11`                : Interval value, defaults to [0.045, 0.945]\n* `colors::Vector{Symbol}`             : Colors for regions, defaults to [:lightgrey, :grey]\n* `stepsize::Float64=0.01`             : Stepsize for boundary accuracy \n* `rescale_axis=true`                  : Display using un-standardized scale         \n\nThis method is primarily intended to display 2 regions, typically predicted values and  the quantile or hpdi region around the mean line. This is specified using the bounds keyword argument, e.g. bounds = [:predicted, :hpdi].\n\nFor the prediction interval, a 3rd parameter needs to be present in the Stan sample DataFrame (dfs) containing the σ value. This symbol needs to be added to linkvars, e.g.\n\nlinkvars = [:a, :bM, :sigma]\n\nFor other options, :quantile and :hpdi, two parameters suffice (typically the itercept and slope parameters).\n\n\n\n\n\n","category":"method"},{"location":"#simulate","page":"Functions","title":"simulate","text":"","category":"section"},{"location":"","page":"Functions","title":"Functions","text":"simulate(df, coefs, var_seq)\nsimulate(df, coefs, var_seq, coefs_ext)","category":"page"},{"location":"#StatisticalRethinking.simulate-Tuple{Any,Any,Any}","page":"Functions","title":"StatisticalRethinking.simulate","text":"simulate\n\nUsed for counterfactual simulations.\n\nsimulate(df, coefs, var_seq)\n\n\nRequired arguments\n\n* `df`                                 : DataFrame with coefficient samples\n* `coefs`                              : Vector of coefficients\n* `var_seq`                            : Input values for simulated effect\n\nReturn values\n\n* `m_sim::NamedTuple`                  : Array with predictions\n\n\n\n\n\n","category":"method"},{"location":"#StatisticalRethinking.simulate-NTuple{4,Any}","page":"Functions","title":"StatisticalRethinking.simulate","text":"simulate\n\nCounterfactual predictions after manipulating a variable.\n\nsimulate(df, coefs, var_seq, coefs_ext)\n\n\nRequired arguments\n\n* `df`                                 : DataFrame with coefficient samples\n* `coefs`                              : Vector of coefficients\n* `var_seq`                            : Input values for simulated effect\n* `ext_coefs`                          : Vector of simulated variable coefficients\n\nReturn values\n\n* `(m_sim, d_sim)`                     : Arrays with predictions\n\n\n\n\n\n","category":"method"},{"location":"srgithub/#Github-organization","page":"StatisticalRethinkingJulia","title":"Github organization","text":"","category":"section"},{"location":"srgithub/","page":"StatisticalRethinkingJulia","title":"StatisticalRethinkingJulia","text":"StatisticalRethinking.jl is part of the broader StatisticalRethinkingJulia Github organization.","category":"page"},{"location":"srgithub/","page":"StatisticalRethinkingJulia","title":"StatisticalRethinkingJulia","text":"Implementations of the models using Stan, DynamicHMC and Turing can currently be found in StanModels, DynamicHMCModels and TuringModels.","category":"page"}]
}

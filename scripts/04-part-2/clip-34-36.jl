# Load Julia packages (libraries) needed  for the snippets in chapter 0

using StatisticalRethinking, StanSample, LinearAlgebra

# ### Snippet 4.26

ProjDir = @__DIR__
df = CSV.read(rel_path("..", "data", "Howell1.csv"), delim=';')
df2 = filter(row -> row[:age] >= 18, df);
first(df2, 5)

heightsmodel = "
// Inferring the mean and std
data {
  int N;
  real<lower=0> h[N];
}
parameters {
  real<lower=0> sigma;
  real<lower=0,upper=250> mu;
}
model {
  // Priors for mu and sigma
  mu ~ normal(178, 20);
  sigma ~ uniform( 0 , 50 );

  // Observed heights
  h ~ normal(mu, sigma);
}
";

# ### Snippet 4.31

sm = SampleModel("heights", heightsmodel);

heightsdata = Dict("N" => length(df2[:, :height]), "h" => df2[:, :height]);

(sample_file, log_file) = stan_sample(sm, data=heightsdata);

if sample_file !== nothing
	println()
	chn = read_samples(sm)
	dfa = DataFrame(chn)

	q = quap(dfa)

	# Check equivalence of Stan samples and Particles.
	mu_range = 152.0:0.01:157.0
	plot(mu_range, ecdf(dfa[:, :mu])(mu_range),
		xlabel="ecdf", ylabel="mu", lab="Stan samples")
	plot!(mu_range, ecdf(systematic_sample(4000, q[:mu]))(mu_range),
		lab="Particles samples")
	savefig("$ProjDir/Fig-34-36.png")
end

# End of `clip-34-36.jl`
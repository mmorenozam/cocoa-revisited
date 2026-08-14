// Hierarchical (partial-pooling) version of baseline_model.stan.
//
// baseline_model.stan fits the 24 kinetic parameters of the cocoa bean
// fermentation ODE separately to each of the 28 trials. This model instead
// treats the 28 trials as exchangeable replicates of one shared kinetic
// process: each trial gets its own parameter vector theta_g, and those
// vectors are drawn from a common population-level log-normal distribution
// (following the partial-pooling approach in Rosenbaum, Raatz, Weithoff,
// Fussmann & Gaedke 2019, "Estimating Parameters From Multiple Time Series
// of Population Dynamics Using Bayesian Inference", Front. Ecol. Evol.,
// 10.3389/fevo.2018.00234, applied there to a predator-prey chemostat
// system fit with Stan/RK45/trajectory matching).
//
// State vector x[1:8] (unchanged from baseline_model.stan):
//   1: reducing sugars   2: sucrose        3: ethanol
//   4: lactic acid       5: acetic acid    6: yeast biomass
//   7: LAB biomass       8: AAB biomass
//
// theta[1:24] index order (unchanged from baseline_model.stan):
//   1-11:  yc1..yc11 (yield coefficients)
//   12-16: mu1..mu5  (max growth rates)
//   17-21: ks1..ks5  (half-saturation constants)
//   22-24: k1..k3    (biomass decay rate constants)
//
// Trials have different numbers of observation times (T[g] in [4,16]), so
// the time/observation arrays are padded out to Tmax and only the first
// T[g] entries of trial g are ever read (see the ragged loop in `model`
// and `generated quantities`); padding values are never touched.
//
// Data and initial states are scaled per trial by that trial's own
// column-wise max (`scl`, computed exactly as in baseline_model.stan),
// which keeps the RK45 integrator well-conditioned across trials whose
// absolute concentrations differ by orders of magnitude. `eps` is a small
// floor (in these scaled units) added before taking logs, since some
// trials have literal zero measurements (fully consumed substrate) that
// a log-normal likelihood cannot otherwise accommodate.

functions {
	vector cbf(real t, vector x, vector theta) {
		vector[8] dxdt;
		real v1 = theta[12]*x[1]*x[6]/(theta[17]+x[1]);
		real v2 = theta[13]*x[2]*x[6]/(theta[18]+x[2]);
		real v3 = theta[14]*x[1]*x[7]/(theta[19]+x[1]);
		real v4 = theta[15]*x[3]*x[8]/(theta[20]+x[3]);
		real v5 = theta[16]*x[4]*x[8]/(theta[21]*x[8]+x[4]);
		real v6 = theta[22]*x[6]*x[3];
		real v7 = theta[23]*x[7]*x[4];
		real v8 = theta[24]*x[8]*x[5]^2;

		dxdt[1] = - theta[1]*v1 - theta[2]*v3;
		dxdt[2] = - theta[3]*v2;
		dxdt[3] = theta[4]*v1 + theta[5]*v2 - theta[6]*v4;
		dxdt[4] = theta[7]*v3 - theta[8]*v5;
		dxdt[5] = theta[9]*v3 + theta[10]*v4 + theta[11]*v5;
		dxdt[6] = v1 + v2 - v6;
		dxdt[7] = v3 - v7;
		dxdt[8] = v4 + v5 - v8;

		return dxdt;
	}
}

data {
	int<lower=1> G;                          // number of trials (28)
	int<lower=1> Tmax;                       // longest trial (16)
	array[G] int<lower=1,upper=Tmax> T;      // observations per trial
	real t0;
	array[G,8] real<lower=0> x0;             // raw initial states
	array[G,Tmax] real ts;                   // observation times, padded
	array[G,Tmax,8] real<lower=0> x;         // raw observed states, padded
	array[G,8] real<lower=0> scl;            // per-trial per-variable scale
	real<lower=0> eps;                       // log-zero floor, scaled units
}

transformed data {
	array[G] vector[8] x0_1;                 // scaled initial states
	array[G,Tmax,8] real xn;                 // scaled observed states (only [1:T[g]] used)
	for (g in 1:G) {
		for (n in 1:8)
			x0_1[g,n] = x0[g,n] / scl[g,n];
		for (t in 1:T[g])
			for (n in 1:8)
				xn[g,t,n] = x[g,t,n] / scl[g,n];
	}
}

parameters {
	vector[24] mu_log_theta;                 // population mean, log scale
	vector<lower=0>[24] sigma_log_theta;     // population sd, log scale (half-normal(0,1))
	matrix[24,G] z_theta;                    // non-centered per-trial deviates
	real<lower=0> sigma;                     // shared residual sd, log scale
}

transformed parameters {
	matrix[24,G] log_theta = rep_matrix(mu_log_theta, G) + diag_pre_multiply(sigma_log_theta, z_theta);
	matrix<lower=0>[24,G] theta = exp(log_theta);
}

model {
	mu_log_theta ~ normal(log(0.5), 1);
	sigma_log_theta ~ std_normal();
	to_vector(z_theta) ~ std_normal();
	sigma ~ cauchy(0, 1);

	for (g in 1:G) {
		array[T[g]] vector[8] x_hat = ode_rk45_tol(cbf, x0_1[g], t0, ts[g,1:T[g]], 1e-6, 1e-6, 10000, theta[,g]);
		for (t in 1:T[g])
			for (n in 1:8)
				xn[g,t,n] + eps ~ lognormal(log(x_hat[t,n] + eps), sigma);
	}
}

generated quantities {
	array[G,Tmax,8] real log_lik;
	array[G,Tmax,8] real x_hat_rep;
	for (g in 1:G) {
		for (t in 1:Tmax)
			for (n in 1:8) {
				log_lik[g,t,n] = -1;
				x_hat_rep[g,t,n] = -1;
			}
		{
			array[T[g]] vector[8] x_hat = ode_rk45_tol(cbf, x0_1[g], t0, ts[g,1:T[g]], 1e-6, 1e-6, 10000, theta[,g]);
			for (t in 1:T[g])
				for (n in 1:8) {
					log_lik[g,t,n] = lognormal_lpdf(xn[g,t,n] + eps | log(x_hat[t,n] + eps), sigma);
					x_hat_rep[g,t,n] = x_hat[t,n] * scl[g,n];
				}
		}
	}
}

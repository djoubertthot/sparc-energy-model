# Milgrom's modified-inertia models on the two tests (2026-10-13)

Model: Milgrom 2022 (PRD 106, 064060), a(w) mu[A(w)/a0] = a_N(w) with A(w_in) = a_in + theta(w_ex/w_in) a_ex, theta(1) = 1, theta decreasing, theta(0) of order 2-3. The external acceleration enters with its OBSERVED value and with a weight set by the frequency ratio. Here mu is the inverse of the RAR nu.

## Test 1: wide binaries (a_ex = 1.8 a0, w_ex/w_in = 0.005-0.02, weight theta(0))

| g_N / a0 | measured ratio | Newton | QUMOND (field theory) | law without EFE | MI theta(0)=1 | MI theta(0)=2 | MI theta(0)=e | MI theta(0)=3 | MI theta(0)=5 |
|---|---|---|---|---|---|---|---|---|---|
| 30.0-1000.0 | 0.995 ± 0.006 | 1.000 | 1.000 | 1.000 | 0.999 | 0.999 | 0.999 | 0.999 | 0.999 |
| 10.0-30.0 | 1.015 ± 0.010 | 1.000 | 1.007 | 1.007 | 1.018 | 1.014 | 1.011 | 1.010 | 1.007 |
| 3.0-10.0 | 1.011 ± 0.011 | 1.000 | 1.047 | 1.048 | 1.064 | 1.043 | 1.033 | 1.030 | 1.016 |
| 1.0-3.0 | 1.046 ± 0.019 | 1.000 | 1.137 | 1.158 | 1.128 | 1.073 | 1.052 | 1.047 | 1.023 |
| 0.3-1.0 | 1.057 ± 0.022 | 1.000 | 1.167 | 1.356 | 1.179 | 1.091 | 1.063 | 1.055 | 1.026 |
| 0.1-0.3 | 1.055 ± 0.041 | 1.000 | 1.168 | 1.663 | 1.207 | 1.098 | 1.067 | 1.058 | 1.027 |
| 0.03-0.1 | 1.297 ± 0.061 | 1.000 | 1.168 | 2.018 | 1.217 | 1.101 | 1.068 | 1.059 | 1.027 |

chi2 on the three bins below 1 a0: Newton 32.4, QUMOND 37.5, law without EFE 549.7, MI best theta(0) = 2.3 (chi2 13.4, 1-sigma range 1.9-2.7); chi2 at theta(0) = 2: 13.8, at 3: 15.1, at 5: 22.1.
In the two clean bins (0.1-1 a0) the data give 1.057 ± 0.022 and 1.055 ± 0.041; MI with theta(0) = 2 predicts 1.091 and 1.098, with theta(0) = 3: 1.055 and 1.058.

## Test 2: lensing at 0.3-1 Mpc (a_ex = nu(e_N) e_N a0 = 0.0105 a0 for the isolated galaxies' e_N = 1.1e-4; weight theta_eff = theta(w_ex/w_in) free)

| theta_eff | chi2 bin 1 | chi2 bin 2 | chi2 bin 3 | chi2 bin 4 | total (48 points) |
|---|---|---|---|---|---|
| 0.00 | 4.5 | 14.8 | 17.0 | 40.6 | 76.8 |
| 0.03 | 4.5 | 14.4 | 17.4 | 40.7 | 77.0 |
| 0.08 | 4.7 | 14.0 | 18.2 | 40.8 | 77.6 |
| 0.20 | 5.5 | 13.2 | 20.4 | 41.2 | 80.3 |
| 0.50 | 9.7 | 12.6 | 27.6 | 43.7 | 93.6 |
| 1.26 | 26.0 | 18.5 | 53.3 | 57.0 | 154.8 |
| 3.16 | 65.9 | 50.3 | 135.4 | 119.7 | 371.2 |
| 7.94 | 130.1 | 129.3 | 315.3 | 325.9 | 900.6 |
| 10.00 | 148.6 | 156.2 | 374.9 | 411.1 | 1090.8 |

- Free amplitude per bin, R <= 1 Mpc: best theta_eff = 0.00, theta_eff < 0.20 at Delta chi2 = 4 (bin 1 alone: < 0.40); Delta chi2 = +53.8 at theta_eff = 1 (7.3 sigma) and +157.2 at theta_eff = 2 (12.5 sigma).
- Internal orbital periods at the lensing radii for bin 1 (V = 134 km/s): 100 kpc: 5 Gyr, 300 kpc: 14 Gyr, 500 kpc: 23 Gyr, 1000 kpc: 46 Gyr. Any large-scale field that varies within a Hubble time (13.8 Gyr) has w_ex > w_in beyond ~300 kpc, so theta_eff = theta(w_ex/w_in > 1) < 1 there.
- With Milgrom's example theta(y) = 2/(1+y^2), the bound theta_eff < 0.20 requires w_ex/w_in > 3.0, i.e. an external field varying faster than 15 Gyr at 1 Mpc: satisfied by the cosmic flow.

## Reading

- Wide binaries: with theta(0) = 2-3 the MI models predict a 5-9 % excess of the relative velocities between 1 and 0.1 a0, against 6 ± 2 % and 6 ± 4 % measured: compatible, where the field theory (QUMOND, +17 %) is 3-5 sigma above and Newton 2-3 sigma below. The reason is structural: in MI the Galactic acceleration enters the inertia of the binary with a weight theta(0) > 1 and isotropically, so the internal dynamics are closer to Newtonian than in modified gravity.
- Lensing at 0.3-1 Mpc: the data bound the external weight to theta_eff < 0.2 (2 sigma), and exclude the field-theory strength theta_eff = 1 at 7 sigma, more strongly than QUMOND was, because the algebraic MI form feels the external acceleration as soon as it is a fraction of the internal one, without the shielding of QUMOND's phantom mass inside r_t. In MI the weight is theta(w_ex/w_in): the model delivers theta_eff < 0.2 if w_ex/w_in > 3, i.e. if the cosmic field varies faster than 15 Gyr at 1 Mpc (orbital period 46 Gyr) and 5 Gyr at 300 kpc (14 Gyr): plausible for the first, marginal for the second. MI thus predicts a weak external-field effect where the lensing finds none, and a strong one in the binaries, where the field theories predict the reverse.
- What this does not establish: MI has no relativistic version, hence no lensing prediction of its own (the test above assumes light follows the circular test orbits); the model's Fourier construction presumes stationary trajectories, which the 50-Gyr orbits at 1 Mpc are not (Milgrom's 'initial conditions problem'); theta(y) is unknown apart from theta(1) = 1; and the binary excess of 5 % is also within the reach of hidden companions. The clusters do not reach MI either: neither lensing nor hydrostatic gas has a modified-inertia prediction. MI survives both tests, the first by predicting it, the second by predicting nothing firm there.

Figure: results/mi_milgrom_predictions.png
# Temporary mass from energy in transit (2026-10-10)

171 SPARC galaxies, inside R_last, every factor generous to the mechanism (see script docstring).

## Today, per galaxy (medians over the sample)

| source | m_temp / M_bary (median) | m_temp / M_missing (median) |
|---|---|---|
| starlight | 5.9e-09 | 2.3e-09 |
| sn_neutrinos | 3.0e-09 | 9.6e-10 |
| CMB | 9.1e-09 | 2.4e-09 |
| CnuB | 6.3e-09 | 1.7e-09 |
| EBL | 4.5e-10 | 1.2e-10 |
| cosmic_rays | 3.5e-08 | 9.4e-09 |
| magnetic | 3.5e-08 | 9.4e-09 |
| gas_thermal | 3.5e-08 | 9.4e-09 |
| m_temp | 1.4e-07 | 4.3e-08 |

Largest m_temp / M_bary in the sample: 2.8e-06 (UGC00128).

## Through cosmic time, median galaxy (UGC07125, M_bary = 6.6e+09 Msun, R_last = 18.7 kpc)

| z | M_bary(z) / Msun | m_temp(z) / Msun | m_temp / M_bary | dominant source |
|---|---|---|---|---|
| 0.0 | 6.6e+09 | 2.5e+03 | 3.8e-07 | cosmic_rays |
| 0.5 | 5.5e+09 | 2.5e+03 | 4.4e-07 | cosmic_rays |
| 1.0 | 4.1e+09 | 2.3e+03 | 5.6e-07 | cosmic_rays |
| 2.0 | 1.7e+09 | 1.7e+03 | 1.0e-06 | CMB |
| 3.0 | 5.9e+08 | 1.5e+03 | 2.5e-06 | CMB |
| 5.0 | 1.1e+08 | 1.9e+03 | 1.8e-05 | CMB |
| 8.0 | 1.8e+07 | 2.9e+03 | 1.6e-04 | CMB |

## If energy were retained for ever (100 % absorption, no re-emission)

- all starlight ever emitted: 1.52e+06 Msun = 2.3e-04 M_bary
- all CMB intercepted since z = 10, galaxy shrinking as 1/(1+z), opacity 1: 2.55e+08 Msun = 3.8e-02 M_bary; with the measured millimetre opacity 1e-5: 3.8e-07 M_bary
- required: M_missing / M_bary = 3.2 (median)

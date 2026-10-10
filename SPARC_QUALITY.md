# SPARC Quality Information

The SPARC catalog (Lelli et al. 2016) provides several indicators that can be used to assess data quality and deprojection reliability:

## Available Columns (Table 1)
- **Inc** : inclination (degrees)
- **e_Inc** : uncertainty on inclination
- **D** : distance (Mpc)
- **e_D** : distance uncertainty
- **Q** : quality flag (1 = high, 2 = medium, 3 = low)
- **Vflat** and its uncertainty

## How to Use Them
1. **Filter** : keep only Q = 1 or Q ≤ 2
2. **Weight** : weight galaxies by 1/e_Inc^{2} or by quality flag
3. **Uncertainty propagation** : Monte-Carlo on inclination and distance to estimate realistic RMS ranges

## Recommendation
Report results for:
- Full sample
- High-quality subsample (Q = 1 and 40° < i < 75°)

This gives more realistic error bars on the model performance.

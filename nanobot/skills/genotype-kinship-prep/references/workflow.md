# Workflow Notes

## Default software

- Format conversion and PCA: `plink`
- Kinship generation: `emmax-kin-intel64`
- Association engine: `emmax-intel64`
- Agents should resolve `emmax-intel64` directly before attempting any fallback binary names

## Key checkpoints

- Confirm that the sample order in the phenotype file matches `tfam`
- Confirm that PCA and kinship are calculated from the same set of samples
- Confirm whether the covariate table already includes an intercept column

## Delivery standard

- The generated `tped/tfam` files can be used directly by `emmax-intel64`
- The kinship file follows the same sample order as `tfam`
- All traits are collected into a clear trait list for batch execution

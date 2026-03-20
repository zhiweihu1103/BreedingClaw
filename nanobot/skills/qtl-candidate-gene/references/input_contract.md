# Input Contract

## Required files

- `gwas_summary.tsv`
- `gene_annotation.tsv`

## `gwas_summary.tsv` columns

- The header must include at least `trait	chrom	pos	snp_id	pvalue`

## `gene_annotation.tsv` columns

- The header must include at least `chrom	start	end	gene_id`

## Optional files

- `ld_table.tsv`
- genotype prefix for LD calculation:
  - `GENO_TFILE_PREFIX` pointing to `tped/tfam`
  - or `GENO_BFILE_PREFIX` pointing to `bed/bim/fam`

## Recommended outputs

- QTL interval table
- Lead SNP table
- Candidate gene table

# Script to create parameters for NUMBAT run with custom annotation


# Read the existing RDS file
a <- readRDS('/projectnb/paxlab/presh/projects/spatial_atac/Data/numbat/results/lowseq/atac_only/run_numbat_params.rds')

# Remove legacy params from source RDS not accepted by numbat 1.5.x
a$mode = NULL
# Add more parameters
a$t = 0.0001
a$ncores = 6        # was n_cores in older numbat versions
a$gamma = 5
a$n_cut = 5

# Save back to the same file
saveRDS(a, '/projectnb/paxlab/presh/projects/spatial_atac/Data/numbat/results/lowseq/atac_only_run2/par_numbat.rds')

# Verify
a <- readRDS('/projectnb/paxlab/presh/projects/spatial_atac/Data/numbat/results/lowseq/atac_only_run2/par_numbat.rds')
str(a)

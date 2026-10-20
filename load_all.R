# Run ONCE per R session from the repository folder:
#   source("load_all.R")
# Loads every function from all topic folders (no package install needed).
files <- list.files(pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
files <- files[grepl("/R/", files)]
for (f in files) source(f)
cat("Loaded", length(files), "files.\n")

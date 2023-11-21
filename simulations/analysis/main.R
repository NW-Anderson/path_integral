#!/usr/bin/env Rscript

library(data.table)
library(logger)
library(ggplot2)

results_base <- "/home/lkirk/Downloads/simulations/out"
## start_freqs_filename <- "start_freqs.csv.gz"
## end_freqs_filename <- "end_freqs.csv.gz"
allele_freqs_filename <- "allele_frequencies.csv.gz"
pop_stats_filename <- "pop_stats.csv.gz"

params <- fread(
    "/home/lkirk/repo/path_integral/simulations/params.txt",
    col.names = c("seed", "mut_rate", "effect_size")
)
setorder(params, mut_rate, effect_size)
params[, group_id := .GRP, by = .(mut_rate, effect_size)]
params[, group_id := as.factor(.SD$group_id)]

load_result_file <- function(base, seed, filename, ...) {
    cmd <- sprintf(
        "tar -zxOf %s/result_%d.tar.gz result_%d/%s", base, seed, seed, filename
    )
    fread(cmd = cmd, ...)
}

load_start_freqs <- function(seed) {
    log_info("loading start freqs for {seed}")
    dt <- load_result_file(results_base, seed, "start_freqs.csv", skip = 2)
    dt[, c("site", "seed") := list(as.integer(rownames(.SD)), seed)]
    dt
}

load_end_freqs <- function(seed) {
    log_info("loading end freqs for {seed}")
    dt <- load_result_file(results_base, seed, "end_freqs.csv", skip = 2)
    dt[, c("site", "seed") := list(as.integer(rownames(.SD)), seed)]
    dt
}

load_pop_stats <- function(seed) {
    log_info("loading pop stats for {seed}")
    dt <- load_result_file(results_base, seed, "popStats.csv")
    dt[, "seed" := seed]
    dt
}

load_effect_size <- function(seed) {
    log_info("loading effect sizes for {seed}")
    dt <- load_result_file(results_base, seed, "effect_sizes.csv")
    dt[, "seed" := seed]
    dt
}

merge_freqs <- function(start, end) {
    start <- melt(
        start,
        id.vars = c("site", "seed"),
        value.name = "start", variable.name = "deme"
    )
    end <- melt(
        end,
        id.vars = c("site", "seed"),
        value.name = "end", variable.name = "deme"
    )
    merged <- start[end, on = c("site", "seed", "deme")]
    merged[, deme := gsub("V", "", deme)]
    setcolorder(merged, c("seed", "deme", "site", "start", "end"))
    merged
}

if (!file.exists(allele_freqs_filename)) {
    start_freqs <- rbindlist(lapply(params$seed, load_start_freqs))
    end_freqs <- rbindlist(lapply(params$seed, load_end_freqs))
    allele_freqs <- merge_freqs(start_freqs, end_freqs)
    log_info("writing {allele_freqs_filename}")
    fwrite(allele_freqs, allele_freqs_filename)
} else {
    log_info("loading cached {allele_freqs_filename}")
    if (exists("allele_freqs")) {
        log_info("not loading {allele_freqs_filename}, allele_freqs exists")
    } else {
        allele_freqs <- fread(allele_freqs_filename)
    }
}

if (!file.exists(pop_stats_filename)) {
    pop_stats <- rbindlist(lapply(params$seed, load_pop_stats))
    log_info("writing {pop_stats_filename}")
    fwrite(pop_stats, pop_stats_filename)
} else {
    log_info("loading cached {pop_stats_filename}")
    if (exists("pop_stats")) {
        log_info("not loading {pop_stats_filename}, pop_stats exists")
    } else {
        pop_stats <- fread(pop_stats_filename)
    }
}


## if(!file.exists(start_freqs_filename)) {
##   start_freqs <- rbindlist(lapply(params[1:10]$seed, load_start_freqs))
##   log_info("writing start_freqs.csv.gz")
##   fwrite(start_freqs, start_freqs_filename)
## } else {
##   log_info("loading cached start_freqs.csv.gz")
##   start_freqs <- fread(start_freqs_filename)
## }

## if(!file.exists(end_freqs_filename)) {
##   end_freqs <- rbindlist(lapply(params[1:10]$seed, load_end_freqs))
##   log_info("writing end_freqs.csv.gz")
##   fwrite(end_freqs, end_freqs_filename)
## } else {
##   log_info("loading cached end_freqs.csv.gz")
##   end_freqs <- fread(end_freqs_filename)
## }


## pop_stats <- fread("tar -zxOf /home/lkirk/Downloads/simulations/out/result_5426.tar.gz result_5426/popStats.csv")
## ggplot(pop_stats, aes(gen, var_pheno)) + geom_line() + xlim(99900,100050)

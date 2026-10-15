# Nested space-filling subsampling from L_Q44_complete.csv ONLY.
# The 200-point response is DoC; the nested 50-point response is u.
# No regression fit or surrogate predictions are used.
# Repeated observations are retained; the design uses unique input settings.
# After design selection, draw one original row per selected setting.
# Shared settings use DoC and u from the SAME original row.
# Place L_Q44_complete.csv in the working directory, then run this script.
# Only base R and the standard stats package are needed.
# Definitions: https://arxiv.org/abs/2305.07202, Section 2.1.

# setwd("~/Downloads")

seed <- 10092026L
xcols <- c("T1", "T2", "t1", "t2")
n_low <- 200L
n_high <- 50L
n_trials_low <- 10000L
n_trials_high <- 10000L

# "maximin": maximize minimum pairwise distance; break ties by minimax.
# "minimax": minimize maximum nearest-design distance; break ties by maximin.
criterion <- "maximin"

# Keep the original responses by default. Set noise_sd_high <- 0.1
# to add independent N(0, 0.1^2) noise to the sampled u values.
noise_sd_low <- 0
noise_sd_high <- 0

read_filtered <- function(file) {
  dat <- read.csv(file, check.names = FALSE)
  required <- c(xcols, "DoC", "u")
  if (!all(required %in% names(dat))) {
    stop(file, ": missing columns: ",
         paste(setdiff(required, names(dat)), collapse = ", "))
  }
  if (!all(vapply(dat[required], is.numeric, logical(1)))) {
    stop(file, ": all input, filter, and response columns must be numeric.")
  }

  # Remove NA/NaN/Inf before applying the original DoC threshold.
  finite <- rowSums(!is.finite(as.matrix(dat[required]))) == 0L
  keep <- finite & !is.na(dat$DoC) & dat$DoC > 0.96
  out <- dat[keep, required, drop = FALSE]
  rownames(out) <- NULL
  message(file, ": retained ", nrow(out), " observations at ",
          nrow(unique(out[xcols])), " distinct input settings.")
  out
}

Generate44 <- read_filtered("L_Q44_complete.csv")

# The input pool comes from ONE table; no between-table matching is needed.
# All responses remain in Generate44 for random replicate selection later.
candidate_X <- unique(Generate44[xcols])
rownames(candidate_X) <- NULL
N <- nrow(candidate_X)
if (N < n_low) {
  stop("Only ", N, " distinct input settings remain in 44 after ",
       "DoC > 0.96 and finite-value checks. At least ", n_low,
       " are required for 200 distinct design points.")
}
stopifnot(n_high >= 2L, n_high <= n_low,
          n_trials_low >= 1L, n_trials_high >= 1L,
          noise_sd_low >= 0, noise_sd_high >= 0)
criterion <- match.arg(criterion, c("maximin", "minimax"))
message("Filtered 44 rows: ", nrow(Generate44),
        "; distinct candidate settings: ", N)

# Use the SAME [0,1] scaling for both designs. Physical output values
# stay unchanged. Constant columns contribute zero to all distances.
X <- as.matrix(candidate_X)
xmin <- apply(X, 2L, min)
xrange <- apply(X, 2L, max) - xmin
xrange[xrange == 0] <- 1
X_unit <- sweep(sweep(X, 2L, xmin, "-"), 2L, xrange, "/")

# Optional distance cache for speed. Avoid an N-by-N allocation for
# large candidate pools; the uncached branch computes distances as needed.
D2 <- if (N <= 3000L) as.matrix(stats::dist(X_unit))^2 else NULL

coverage_distance <- function(ids) {
  # Discrete minimax: ALL distinct filtered 44 inputs are the reference set,
  # including when assessing the nested 50-point design.
  nearest2 <- rep(Inf, N)
  for (j in ids) {
    d2 <- if (!is.null(D2)) D2[, j] else
      rowSums(sweep(X_unit, 2L, X_unit[j, ], "-")^2)
    nearest2 <- pmin(nearest2, d2)
  }
  sqrt(max(nearest2))
}

best_random_subset <- function(pool, size, n_trials,
                               criterion = "maximin") {
  criterion <- match.arg(criterion, c("maximin", "minimax"))
  stopifnot(size >= 2L, size <= length(pool), !anyDuplicated(pool),
            all(pool %in% seq_len(N)), n_trials >= 1L)
  best_ids <- NULL
  best_maximin <- -Inf
  best_minimax <- Inf
  tol <- 1e-12

  # A full-pool selection is unique; do not repeat the same trial.
  if (size == length(pool)) n_trials <- 1L
  for (b in seq_len(n_trials)) {
    ids <- pool[sample.int(length(pool), size, replace = FALSE)]
    current_maximin <- min(stats::dist(X_unit[ids, , drop = FALSE]))

    # Under maximin priority, inferior separation cannot win.
    # Skip the more expensive coverage calculation for these trials.
    if (criterion == "maximin" &&
        current_maximin < best_maximin - tol) next
    current_minimax <- coverage_distance(ids)

    better <- if (criterion == "maximin") {
      current_maximin > best_maximin + tol ||
        (abs(current_maximin - best_maximin) <= tol &&
           current_minimax < best_minimax - tol)
    } else {
      current_minimax < best_minimax - tol ||
        (abs(current_minimax - best_minimax) <= tol &&
           current_maximin > best_maximin + tol)
    }
    if (better) {
      best_ids <- ids
      best_maximin <- current_maximin
      best_minimax <- current_minimax
    }
  }
  list(ids = best_ids, maximin = best_maximin, minimax = best_minimax,
       n_trials = n_trials)
}

set.seed(seed)

# Stage 1: choose 200 points from the distinct 44 inputs (response: DoC).
low_design <- best_random_subset(seq_len(N), n_low,
                                 n_trials_low, criterion)
message("Selected the 200-point design; searching for its nested subset.")

# Stage 2: hold the 200-point design fixed; choose 50 WITHIN it (response: u).
# This is a conditional, two-stage random search, not joint optimization.
high_design <- best_random_subset(low_design$ids, n_high,
                                  n_trials_high, criterion)

# Put the shared 50 settings in the first 50 rows of the 200-row file,
# in exactly the same order as the high-fidelity file.
id_high <- high_design$ids
id_low <- c(id_high, setdiff(low_design$ids, id_high))

# For each selected setting, randomly choose ONE of its original observations.
# Sampling settings first prevents heavily replicated settings from getting
# extra weight in the space-filling search. Responses do not affect selection.
draw_observation_rows <- function(dat, selected_X) {
  vapply(seq_len(nrow(selected_X)), function(i) {
    same_X <- rep(TRUE, nrow(dat))
    for (nm in xcols) {
      same_X <- same_X & dat[[nm]] == selected_X[[nm]][i]
    }
    rows <- which(same_X)
    if (!length(rows)) stop("A selected setting has no original observation.")
    # Index with sample.int: sample(rows, 1) behaves differently when rows
    # contains just one numeric row index.
    rows[sample.int(length(rows), 1L)]
  }, integer(1))
}

# Draw ONE original 44 row per selected input. The first 50 row indices
# supply both responses, preserving their pairing at the shared inputs.
set.seed(seed + 1L)
selected_rows44 <- draw_observation_rows(Generate44,
                                        candidate_X[id_low, , drop = FALSE])
selected_rows_u <- selected_rows44[seq_len(n_high)]
ball4 <- data.frame(candidate_X[id_low, , drop = FALSE],
                    y = Generate44$DoC[selected_rows44])
ball2 <- data.frame(candidate_X[id_high, , drop = FALSE],
                    y = Generate44$u[selected_rows_u])
rownames(ball4) <- NULL
rownames(ball2) <- NULL

# Optional observation noise is separate from design selection.
set.seed(seed + 2L)
if (noise_sd_low > 0) {
  ball4$y <- ball4$y + rnorm(n_low, sd = noise_sd_low)
}
if (noise_sd_high > 0) {
  ball2$y <- ball2$y + rnorm(n_high, sd = noise_sd_high)
}

# Verify counts and exact nesting before saving.
stopifnot(nrow(ball4) == n_low, nrow(ball2) == n_high,
          !anyDuplicated(ball4[xcols]), !anyDuplicated(ball2[xcols]),
          identical(ball2[xcols], ball4[seq_len(n_high), xcols]))

write.csv(ball4, "Curing_44DoC_subsample200.csv", row.names = FALSE)
write.csv(ball2, "Curing_44u_subsample50_nested.csv", row.names = FALSE)

quality <- data.frame(
  design = c("44 DoC: 200 points", "44 u: 50 nested points"),
  maximin = c(low_design$maximin, high_design$maximin),
  minimax = c(low_design$minimax, high_design$minimax)
)
print(quality, row.names = FALSE)
message("Saved both CSV files. Higher maximin and lower minimax are better.")

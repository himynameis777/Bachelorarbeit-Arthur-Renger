# Habitat-loss experiment: analysis and plots
#
# Run from RStudio after setting the working directory to the CSV folder,
# or from Terminal:
# Rscript plot_habitat_results.R "/path/to/csv/files" "/path/to/output"

required_packages <- c(
  "readr", "dplyr", "tidyr", "ggplot2", "purrr", "scales"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
  library(scales)
})

# ----------------------------------------------------------------------
# File locations
# ----------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

data_dir <- if (length(args) >= 1) {
  normalizePath(args[[1]], mustWork = TRUE)
} else {
  getwd()
}

output_dir <- if (length(args) >= 2) {
  args[[2]]
} else {
  file.path(data_dir, "habitat_result_plots")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Selects the largest matching file. This handles names such as
# habitat_timeseries(2).csv.
find_file <- function(pattern, description) {
  candidates <- list.files(
    data_dir,
    pattern = pattern,
    full.names = TRUE
  )

  if (length(candidates) == 0) {
    stop(
      "Could not find ", description,
      " in: ", data_dir,
      call. = FALSE
    )
  }

  selected <- candidates[[which.max(file.info(candidates)$size)]]
  message("Using ", description, ": ", basename(selected))
  selected
}

run_file <- find_file(
  "^habitat_run_summary.*\\.csv$",
  "run summary"
)

timeseries_file <- find_file(
  "^habitat_timeseries.*\\.csv$",
  "time series"
)

treatment_file <- find_file(
  "^habitat_treatment_summary.*\\.csv$",
  "treatment summary"
)

seed_file <- find_file(
  "^habitat_seed_table.*\\.csv$",
  "seed table"
)

configuration_file <- find_file(
  "^habitat_configuration.*\\.csv$",
  "configuration"
)

runs <- read_csv(run_file, show_col_types = FALSE)
timeseries <- read_csv(timeseries_file, show_col_types = FALSE)
supplied_treatment_summary <- read_csv(
  treatment_file,
  show_col_types = FALSE
)
seed_table <- read_csv(seed_file, show_col_types = FALSE)
configuration <- read_csv(
  configuration_file,
  show_col_types = FALSE
)

# ----------------------------------------------------------------------
# Data preparation and validation
# ----------------------------------------------------------------------

as_boolean <- function(x) {
  if (is.logical(x)) {
    return(x)
  }

  tolower(trimws(as.character(x))) %in%
    c("true", "t", "1", "yes")
}

q_value <- function(x, probability) {
  if (all(is.na(x))) {
    return(NA_real_)
  }

  as.numeric(
    quantile(
      x,
      probability,
      na.rm = TRUE,
      names = FALSE
    )
  )
}

wilson_lower <- function(events, n, z = qnorm(0.975)) {
  probability <- events / n

  centre <- (
    probability + z^2 / (2 * n)
  ) / (
    1 + z^2 / n
  )

  half_width <- z * sqrt(
    probability * (1 - probability) / n +
      z^2 / (4 * n^2)
  ) / (
    1 + z^2 / n
  )

  pmax(0, centre - half_width)
}

wilson_upper <- function(events, n, z = qnorm(0.975)) {
  probability <- events / n

  centre <- (
    probability + z^2 / (2 * n)
  ) / (
    1 + z^2 / n
  )

  half_width <- z * sqrt(
    probability * (1 - probability) / n +
      z^2 / (4 * n^2)
  ) / (
    1 + z^2 / n
  )

  pmin(1, centre + half_width)
}

runs <- runs %>%
  mutate(
    prey_extinct = as_boolean(prey_extinct),
    predator_extinct = as_boolean(predator_extinct),
    survived_to_max_ticks = as_boolean(survived_to_max_ticks),
    stable_dynamics = as_boolean(stable_dynamics),
    loss = fixed_habitat_loss_percent,
    run_id = paste(loss, repetition, sep = "_")
  )

timeseries <- timeseries %>%
  mutate(
    loss = fixed_habitat_loss_percent,
    run_id = paste(loss, repetition, sep = "_")
  )

loss_levels <- sort(unique(runs$loss))
n_repetitions <- runs %>%
  count(loss) %>%
  summarise(n = min(n)) %>%
  pull(n)

treatment_tick <- unique(runs$treatment_tick)
maximum_tick <- max(runs$final_tick, na.rm = TRUE)

if (length(treatment_tick) != 1) {
  stop("More than one treatment tick was found.")
}

if (n_distinct(runs$config_id) != 1) {
  stop("The experiment contains more than one configuration.")
}

if (anyDuplicated(runs[c("loss", "repetition")]) > 0) {
  stop("Duplicate treatment/repetition records were found.")
}

if (
  anyDuplicated(
    timeseries[c("loss", "repetition", "tick")]
  ) > 0
) {
  stop("Duplicate treatment/repetition/tick records were found.")
}

runs_per_treatment <- runs %>%
  count(loss, name = "completed_runs")

if (n_distinct(runs_per_treatment$completed_runs) != 1) {
  stop("Treatments contain different numbers of repetitions.")
}

if (n_repetitions != nrow(seed_table)) {
  stop("The seed table does not match the repetitions.")
}

message(
  "Detected ",
  length(loss_levels),
  " habitat-loss treatments, ",
  n_repetitions,
  " paired repetitions per treatment, and ",
  nrow(runs),
  " total runs."
)

# ----------------------------------------------------------------------
# Plot appearance
# ----------------------------------------------------------------------

theme_set(
  theme_minimal(base_size = 11) +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(
        face = "bold",
        size = 13
      ),
      plot.subtitle = element_text(
        size = 9.5,
        colour = "grey30"
      ),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      plot.caption = element_text(
        size = 8,
        colour = "grey35",
        hjust = 0
      )
    )
)

population_colours <- c(
  "Prey" = "#0072B2",
  "Predators" = "#D55E00"
)

mortality_colours <- c(
  "Predation" = "#CC79A7",
  "Starvation" = "#E69F00",
  "Crowding" = "#009E73"
)

selected_levels <- intersect(
  c(0, 25, 50, 75, 100),
  loss_levels
)

selected_labels <- paste0(selected_levels, "%")

selected_colours <- setNames(
  c(
    "#202020",
    "#0072B2",
    "#009E73",
    "#E69F00",
    "#D55E00"
  )[seq_along(selected_levels)],
  selected_labels
)

plots <- list()

save_plot <- function(
  filename,
  plot,
  width = 11,
  height = 7
) {
  ggsave(
    filename = file.path(
      output_dir,
      paste0(filename, ".png")
    ),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    bg = "white"
  )

  plots[[filename]] <<- plot
}

# ----------------------------------------------------------------------
# 1. Intended versus realised habitat loss
# ----------------------------------------------------------------------

p <- ggplot(
  runs,
  aes(
    x = loss,
    y = actual_habitat_loss_percent
  )
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    colour = "grey50"
  ) +
  geom_point(
    alpha = 0.3,
    size = 1.4,
    colour = "#0072B2",
    position = position_jitter(width = 0.6)
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    colour = "#D55E00",
    size = 2.5
  ) +
  coord_equal() +
  scale_x_continuous(
    breaks = seq(0, 100, 10),
    limits = c(0, 100)
  ) +
  scale_y_continuous(
    breaks = seq(0, 100, 10),
    limits = c(0, 100)
  ) +
  labs(
    title = "Intended and realised habitat loss",
    subtitle = paste(
      "Orange points show treatment means;",
      "the dashed line indicates exact agreement"
    ),
    x = "Intended habitat loss (%)",
    y = "Realised habitat loss (%)"
  )

save_plot(
  "01_intended_vs_realised_habitat_loss",
  p
)

# ----------------------------------------------------------------------
# 2. Habitat-loss implementation error
# ----------------------------------------------------------------------

loss_error <- runs %>%
  mutate(
    error_pp =
      actual_habitat_loss_percent - loss
  )

p <- ggplot(
  loss_error,
  aes(
    x = factor(loss),
    y = error_pp
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey50"
  ) +
  geom_boxplot(
    fill = "#56B4E9",
    alpha = 0.65,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.15,
    alpha = 0.25,
    size = 0.8
  ) +
  labs(
    title = "Realised habitat-loss error",
    subtitle = paste(
      "Positive values indicate slightly more",
      "habitat loss than requested"
    ),
    x = "Intended habitat loss (%)",
    y = paste(
      "Realised minus intended loss",
      "(percentage points)"
    )
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

save_plot(
  "02_habitat_loss_error",
  p
)

# ----------------------------------------------------------------------
# 3. Terminal population distributions
# ----------------------------------------------------------------------

final_population_long <- runs %>%
  select(
    loss,
    repetition,
    seed,
    final_tick,
    final_preys,
    final_predators
  ) %>%
  pivot_longer(
    cols = c(
      final_preys,
      final_predators
    ),
    names_to = "population",
    values_to = "abundance"
  ) %>%
  mutate(
    population = recode(
      population,
      final_preys = "Prey",
      final_predators = "Predators"
    )
  )

p <- ggplot(
  final_population_long,
  aes(
    x = factor(loss),
    y = abundance,
    fill = population
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.65
  ) +
  geom_jitter(
    width = 0.16,
    alpha = 0.2,
    size = 0.65
  ) +
  facet_wrap(
    ~population,
    scales = "free_y",
    ncol = 1
  ) +
  scale_fill_manual(
    values = population_colours,
    guide = "none"
  ) +
  labs(
    title = "Terminal population abundance",
    subtitle = paste(
      "Values are recorded at tick 3000",
      "or at the first population extinction"
    ),
    x = "Habitat loss (%)",
    y = "Individuals"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

save_plot(
  "03_terminal_population_distributions",
  p
)

# ----------------------------------------------------------------------
# 4. Paired differences from the 0% control
# ----------------------------------------------------------------------

control_values <- runs %>%
  filter(loss == 0) %>%
  select(
    repetition,
    seed,
    control_preys = final_preys,
    control_predators = final_predators
  )

paired_response <- runs %>%
  left_join(
    control_values,
    by = c("repetition", "seed")
  ) %>%
  transmute(
    loss,
    repetition,
    seed,
    prey_difference =
      final_preys - control_preys,
    predator_difference =
      final_predators - control_predators,
    prey_relative_change = if_else(
      control_preys > 0,
      100 * prey_difference / control_preys,
      NA_real_
    ),
    predator_relative_change = if_else(
      control_predators > 0,
      100 * predator_difference /
        control_predators,
      NA_real_
    )
  )

paired_difference_long <- paired_response %>%
  select(
    loss,
    repetition,
    seed,
    prey_difference,
    predator_difference
  ) %>%
  pivot_longer(
    cols = c(
      prey_difference,
      predator_difference
    ),
    names_to = "population",
    values_to = "difference"
  ) %>%
  mutate(
    population = recode(
      population,
      prey_difference = "Prey",
      predator_difference = "Predators"
    )
  )

p <- ggplot(
  paired_difference_long,
  aes(
    x = loss,
    y = difference,
    group = seed
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey50"
  ) +
  geom_line(
    alpha = 0.12,
    colour = "grey35"
  ) +
  stat_summary(
    aes(
      group = 1,
      colour = population
    ),
    fun = median,
    geom = "line",
    linewidth = 1.2
  ) +
  stat_summary(
    aes(
      group = 1,
      colour = population
    ),
    fun = median,
    geom = "point",
    size = 1.8
  ) +
  facet_wrap(
    ~population,
    scales = "free_y",
    ncol = 1
  ) +
  scale_colour_manual(
    values = population_colours,
    guide = "none"
  ) +
  scale_x_continuous(
    breaks = seq(0, 100, 10)
  ) +
  labs(
    title = "Paired differences from the 0% control",
    subtitle = paste(
      "Thin lines connect runs using the same seed;",
      "coloured lines show medians"
    ),
    x = "Habitat loss (%)",
    y = "Difference from matched control"
  )

save_plot(
  "04_paired_population_differences",
  p
)

# ----------------------------------------------------------------------
# 5. Relative paired treatment response
# ----------------------------------------------------------------------

paired_relative_long <- paired_response %>%
  select(
    loss,
    repetition,
    seed,
    prey_relative_change,
    predator_relative_change
  ) %>%
  pivot_longer(
    cols = c(
      prey_relative_change,
      predator_relative_change
    ),
    names_to = "population",
    values_to = "relative_change"
  ) %>%
  mutate(
    population = recode(
      population,
      prey_relative_change = "Prey",
      predator_relative_change = "Predators"
    )
  )

p <- ggplot(
  paired_relative_long,
  aes(
    x = factor(loss),
    y = relative_change,
    fill = population
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey50"
  ) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.65
  ) +
  facet_wrap(
    ~population,
    scales = "free_y",
    ncol = 1
  ) +
  scale_fill_manual(
    values = population_colours,
    guide = "none"
  ) +
  labs(
    title = "Relative change from the matched control",
    subtitle = paste(
      "Each treatment is compared with the",
      "0% run using the same seed"
    ),
    x = "Habitat loss (%)",
    y = "Change from matched control (%)"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

save_plot(
  "05_paired_relative_population_change",
  p
)

# ----------------------------------------------------------------------
# 6. Terminal prey-predator relationship
# ----------------------------------------------------------------------

p <- ggplot(
  runs,
  aes(
    x = final_preys,
    y = final_predators,
    colour = loss
  )
) +
  geom_point(
    alpha = 0.55,
    size = 1.5
  ) +
  scale_colour_viridis_c(
    option = "C",
    end = 0.95,
    breaks = seq(0, 100, 20)
  ) +
  labs(
    title = "Terminal prey and predator abundance",
    subtitle = "Each point represents one run",
    x = "Terminal prey abundance",
    y = "Terminal predator abundance",
    colour = "Habitat loss (%)"
  )

save_plot(
  "06_terminal_prey_predator_relationship",
  p
)

# ----------------------------------------------------------------------
# 7. Extinction probabilities
# ----------------------------------------------------------------------

extinction_summary <- runs %>%
  group_by(loss) %>%
  summarise(
    n = n(),
    Prey = sum(prey_extinct),
    Predators = sum(predator_extinct),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(Prey, Predators),
    names_to = "population",
    values_to = "events"
  ) %>%
  mutate(
    probability = events / n,
    lower = wilson_lower(events, n),
    upper = wilson_upper(events, n)
  )

p <- ggplot(
  extinction_summary,
  aes(
    x = loss,
    y = probability,
    colour = population
  )
) +
  geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper,
      fill = population
    ),
    alpha = 0.14,
    colour = NA
  ) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_colour_manual(
    values = population_colours
  ) +
  scale_fill_manual(
    values = population_colours
  ) +
  scale_x_continuous(
    breaks = seq(0, 100, 10)
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1)
  ) +
  labs(
    title = "Extinction probability",
    subtitle = paste(
      "Shading shows 95% Wilson intervals;",
      "runs stop at the first population extinction"
    ),
    x = "Habitat loss (%)",
    y = "Extinction probability",
    colour = "Population",
    fill = "Population",
    caption = paste0(
      "n = ",
      n_repetitions,
      " paired runs per treatment."
    )
  )

save_plot(
  "07_extinction_probability",
  p
)

# ----------------------------------------------------------------------
# 8. Co-persistence probability
# ----------------------------------------------------------------------

persistence_summary <- runs %>%
  group_by(loss) %>%
  summarise(
    n = n(),
    events = sum(survived_to_max_ticks),
    .groups = "drop"
  ) %>%
  mutate(
    probability = events / n,
    lower = wilson_lower(events, n),
    upper = wilson_upper(events, n)
  )

p <- ggplot(
  persistence_summary,
  aes(
    x = loss,
    y = probability
  )
) +
  geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper
    ),
    fill = "#56B4E9",
    alpha = 0.22
  ) +
  geom_line(
    colour = "#0072B2",
    linewidth = 1
  ) +
  geom_point(
    colour = "#0072B2",
    size = 2
  ) +
  scale_x_continuous(
    breaks = seq(0, 100, 10)
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1)
  ) +
  labs(
    title = paste(
      "Probability that both populations",
      "persisted to tick 3000"
    ),
    subtitle = "Shading shows 95% Wilson intervals",
    x = "Habitat loss (%)",
    y = "Co-persistence probability"
  )

save_plot(
  "08_co_persistence_probability",
  p
)

# ----------------------------------------------------------------------
# 9. Time to extinction
# ----------------------------------------------------------------------

extinction_runs <- runs %>%
  filter(
    extinct_population %in%
      c("preys", "predators"),
    extinction_tick >= treatment_tick
  ) %>%
  mutate(
    extinct_population = recode(
      extinct_population,
      preys = "Prey",
      predators = "Predators"
    ),
    ticks_after_treatment =
      extinction_tick - treatment_tick
  )

if (nrow(extinction_runs) > 0) {
  p <- ggplot(
    extinction_runs,
    aes(
      x = factor(loss),
      y = ticks_after_treatment,
      colour = extinct_population
    )
  ) +
    geom_boxplot(
      aes(
        group = interaction(
          loss,
          extinct_population
        )
      ),
      outlier.shape = NA,
      alpha = 0.3
    ) +
    geom_jitter(
      width = 0.16,
      alpha = 0.55,
      size = 1.1
    ) +
    facet_wrap(
      ~extinct_population,
      scales = "free_x"
    ) +
    scale_colour_manual(
      values = population_colours,
      guide = "none"
    ) +
    labs(
      title = "Time from habitat loss to extinction",
      subtitle = paste(
        "Only runs with an observed extinction",
        "are included"
      ),
      x = "Habitat loss (%)",
      y = "Ticks after treatment"
    ) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      )
    )

  save_plot(
    "09_time_to_extinction",
    p
  )
}

# ----------------------------------------------------------------------
# 10. Co-persistence through time
# ----------------------------------------------------------------------

persistence_curve <- crossing(
  loss = selected_levels,
  tick = seq(
    treatment_tick,
    maximum_tick,
    by = 10
  )
) %>%
  left_join(
    runs %>%
      select(
        loss,
        repetition,
        final_tick
      ),
    by = "loss"
  ) %>%
  group_by(loss, tick) %>%
  summarise(
    proportion_active =
      mean(final_tick >= tick),
    .groups = "drop"
  ) %>%
  mutate(
    loss_label = factor(
      paste0(loss, "%"),
      levels = selected_labels
    )
  )

p <- ggplot(
  persistence_curve,
  aes(
    x = tick,
    y = proportion_active,
    colour = loss_label
  )
) +
  geom_line(linewidth = 1) +
  scale_colour_manual(
    values = selected_colours
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1)
  ) +
  labs(
    title = "Persistence after habitat loss",
    subtitle = paste(
      "A run remains active until the first",
      "population extinction"
    ),
    x = "Tick",
    y = "Runs with both populations present",
    colour = "Habitat loss"
  )

save_plot(
  "10_co_persistence_over_time",
  p
)

# ----------------------------------------------------------------------
# 11. Population trajectories
# ----------------------------------------------------------------------

population_long <- timeseries %>%
  select(
    loss,
    repetition,
    seed,
    tick,
    preys,
    predators
  ) %>%
  pivot_longer(
    cols = c(preys, predators),
    names_to = "population",
    values_to = "abundance"
  ) %>%
  mutate(
    population = recode(
      population,
      preys = "Prey",
      predators = "Predators"
    )
  )

population_summary <- population_long %>%
  group_by(
    loss,
    tick,
    population
  ) %>%
  summarise(
    n_active = n(),
    median = median(abundance),
    q10 = q_value(abundance, 0.10),
    q90 = q_value(abundance, 0.90),
    .groups = "drop"
  )

selected_population_summary <-
  population_summary %>%
  filter(loss %in% selected_levels) %>%
  mutate(
    loss_label = factor(
      paste0(loss, "%"),
      levels = selected_labels
    )
  )

p <- ggplot(
  selected_population_summary,
  aes(
    x = tick,
    y = median,
    colour = loss_label,
    fill = loss_label
  )
) +
  geom_ribbon(
    aes(
      ymin = q10,
      ymax = q90
    ),
    alpha = 0.10,
    colour = NA
  ) +
  geom_line(linewidth = 0.9) +
  geom_vline(
    xintercept = treatment_tick,
    linetype = "dashed",
    colour = "grey35"
  ) +
  facet_wrap(
    ~population,
    scales = "free_y",
    ncol = 1
  ) +
  scale_colour_manual(
    values = selected_colours
  ) +
  scale_fill_manual(
    values = selected_colours
  ) +
  labs(
    title = "Population trajectories",
    subtitle = paste(
      "Lines show medians and ribbons the",
      "10th–90th percentiles among active runs"
    ),
    x = "Tick",
    y = "Individuals",
    colour = "Habitat loss",
    fill = "Habitat loss",
    caption = paste(
      "The dashed line marks habitat loss.",
      "Runs end when the first population becomes extinct."
    )
  )

save_plot(
  "11_population_trajectories",
  p
)

# ----------------------------------------------------------------------
# 12. Population heatmaps
# ----------------------------------------------------------------------

regular_population_summary <-
  population_summary %>%
  filter(tick %% 10 == 0)

p <- regular_population_summary %>%
  filter(population == "Prey") %>%
  ggplot(
    aes(
      x = tick,
      y = loss,
      fill = median
    )
  ) +
  geom_tile() +
  geom_vline(
    xintercept = treatment_tick,
    linetype = "dashed",
    colour = "white"
  ) +
  scale_fill_viridis_c(
    option = "C",
    trans = "sqrt"
  ) +
  scale_y_continuous(
    breaks = seq(0, 100, 10)
  ) +
  labs(
    title = "Median prey abundance",
    subtitle = paste(
      "Medians use runs still active",
      "at each recorded tick"
    ),
    x = "Tick",
    y = "Habitat loss (%)",
    fill = "Median prey"
  )

save_plot(
  "12_prey_abundance_heatmap",
  p
)

p <- regular_population_summary %>%
  filter(population == "Predators") %>%
  ggplot(
    aes(
      x = tick,
      y = loss,
      fill = median
    )
  ) +
  geom_tile() +
  geom_vline(
    xintercept = treatment_tick,
    linetype = "dashed",
    colour = "white"
  ) +
  scale_fill_viridis_c(
    option = "D",
    trans = "sqrt"
  ) +
  scale_y_continuous(
    breaks = seq(0, 100, 10)
  ) +
  labs(
    title = "Median predator abundance",
    subtitle = paste(
      "Medians use runs still active",
      "at each recorded tick"
    ),
    x = "Tick",
    y = "Habitat loss (%)",
    fill = "Median predators"
  )

save_plot(
  "13_predator_abundance_heatmap",
  p
)

# ----------------------------------------------------------------------
# 13. Active-run heatmap
# ----------------------------------------------------------------------

active_fraction <- timeseries %>%
  filter(tick %% 10 == 0) %>%
  count(
    loss,
    tick,
    name = "active_runs"
  ) %>%
  mutate(
    active_fraction =
      active_runs / n_repetitions
  )

p <- ggplot(
  active_fraction,
  aes(
    x = tick,
    y = loss,
    fill = active_fraction
  )
) +
  geom_tile() +
  geom_vline(
    xintercept = treatment_tick,
    linetype = "dashed",
    colour = "white"
  ) +
  scale_fill_viridis_c(
    option = "B",
    limits = c(0, 1),
    labels = percent_format(accuracy = 1)
  ) +
  scale_y_continuous(
    breaks = seq(0, 100, 10)
  ) +
  labs(
    title = "Fraction of runs still active",
    subtitle = paste(
      "A run stops when the first",
      "population becomes extinct"
    ),
    x = "Tick",
    y = "Habitat loss (%)",
    fill = "Active runs"
  )

save_plot(
  "14_active_run_fraction",
  p
)

# ----------------------------------------------------------------------
# 14. Minimum post-treatment populations
# ----------------------------------------------------------------------

population_minima <- timeseries %>%
  filter(tick >= treatment_tick) %>%
  group_by(
    loss,
    repetition,
    seed
  ) %>%
  summarise(
    min_preys = min(preys),
    min_predators = min(predators),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(
      min_preys,
      min_predators
    ),
    names_to = "population",
    values_to = "minimum"
  ) %>%
  mutate(
    population = recode(
      population,
      min_preys = "Prey",
      min_predators = "Predators"
    )
  )

p <- ggplot(
  population_minima,
  aes(
    x = factor(loss),
    y = minimum,
    fill = population
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.65
  ) +
  facet_wrap(
    ~population,
    scales = "free_y",
    ncol = 1
  ) +
  scale_fill_manual(
    values = population_colours,
    guide = "none"
  ) +
  labs(
    title = "Minimum population after treatment",
    subtitle = paste(
      "Minimum observed value between",
      "tick 2000 and run termination"
    ),
    x = "Habitat loss (%)",
    y = "Minimum individuals"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

save_plot(
  "15_post_treatment_population_minima",
  p
)

# ----------------------------------------------------------------------
# 15. Hunger trajectories
# ----------------------------------------------------------------------

hunger_summary <- timeseries %>%
  select(
    loss,
    tick,
    mean_prey_hunger,
    mean_predator_hunger
  ) %>%
  pivot_longer(
    cols = c(
      mean_prey_hunger,
      mean_predator_hunger
    ),
    names_to = "population",
    values_to = "hunger"
  ) %>%
  mutate(
    population = recode(
      population,
      mean_prey_hunger = "Prey",
      mean_predator_hunger = "Predators"
    )
  ) %>%
  group_by(
    loss,
    tick,
    population
  ) %>%
  summarise(
    median = median(hunger, na.rm = TRUE),
    q10 = q_value(hunger, 0.10),
    q90 = q_value(hunger, 0.90),
    .groups = "drop"
  ) %>%
  filter(loss %in% selected_levels) %>%
  mutate(
    loss_label = factor(
      paste0(loss, "%"),
      levels = selected_labels
    )
  )

p <- ggplot(
  hunger_summary,
  aes(
    x = tick,
    y = median,
    colour = loss_label,
    fill = loss_label
  )
) +
  geom_ribbon(
    aes(
      ymin = q10,
      ymax = q90
    ),
    alpha = 0.10,
    colour = NA
  ) +
  geom_line(linewidth = 0.9) +
  geom_vline(
    xintercept = treatment_tick,
    linetype = "dashed",
    colour = "grey35"
  ) +
  facet_wrap(
    ~population,
    scales = "free_y",
    ncol = 1
  ) +
  scale_colour_manual(
    values = selected_colours
  ) +
  scale_fill_manual(
    values = selected_colours
  ) +
  labs(
    title = "Mean hunger through time",
    subtitle = paste(
      "Lines show medians and ribbons the",
      "10th–90th percentiles among active runs"
    ),
    x = "Tick",
    y = "Mean hunger",
    colour = "Habitat loss",
    fill = "Habitat loss"
  )

save_plot(
  "16_hunger_trajectories",
  p
)

# ----------------------------------------------------------------------
# 16. Prey on destroyed habitat
# ----------------------------------------------------------------------

trapped_summary <- timeseries %>%
  filter(
    loss %in% selected_levels,
    tick >= treatment_tick
  ) %>%
  mutate(
    trapped_share = if_else(
      preys > 0,
      100 * trapped_preys / preys,
      NA_real_
    )
  ) %>%
  select(
    loss,
    tick,
    trapped_preys,
    trapped_share
  ) %>%
  pivot_longer(
    cols = c(
      trapped_preys,
      trapped_share
    ),
    names_to = "metric",
    values_to = "value"
  ) %>%
  group_by(
    loss,
    tick,
    metric
  ) %>%
  summarise(
    median = median(value, na.rm = TRUE),
    q10 = q_value(value, 0.10),
    q90 = q_value(value, 0.90),
    .groups = "drop"
  ) %>%
  mutate(
    loss_label = factor(
      paste0(loss, "%"),
      levels = selected_labels
    ),
    metric = recode(
      metric,
      trapped_preys = "Number of prey",
      trapped_share = "Share of prey (%)"
    )
  )

p <- ggplot(
  trapped_summary,
  aes(
    x = tick,
    y = median,
    colour = loss_label,
    fill = loss_label
  )
) +
  geom_ribbon(
    aes(
      ymin = q10,
      ymax = q90
    ),
    alpha = 0.10,
    colour = NA
  ) +
  geom_line(linewidth = 0.9) +
  facet_wrap(
    ~metric,
    scales = "free_y",
    ncol = 1
  ) +
  scale_colour_manual(
    values = selected_colours
  ) +
  scale_fill_manual(
    values = selected_colours
  ) +
  labs(
    title = "Prey occupying destroyed habitat",
    subtitle = paste(
      "Lines show medians and ribbons the",
      "10th–90th percentiles among active runs"
    ),
    x = "Tick",
    y = NULL,
    colour = "Habitat loss",
    fill = "Habitat loss"
  )

save_plot(
  "17_trapped_prey",
  p
)

# ----------------------------------------------------------------------
# 17. Post-treatment prey mortality
# ----------------------------------------------------------------------

deaths_at_treatment <- timeseries %>%
  filter(tick == treatment_tick) %>%
  select(
    loss,
    repetition,
    seed,
    baseline_predation =
      prey_deaths_predation,
    baseline_starvation =
      prey_deaths_starvation,
    baseline_crowding =
      prey_deaths_crowding
  )

post_treatment_deaths <- runs %>%
  select(
    loss,
    repetition,
    seed,
    final_predation =
      prey_deaths_predation,
    final_starvation =
      prey_deaths_starvation,
    final_crowding =
      prey_deaths_crowding
  ) %>%
  left_join(
    deaths_at_treatment,
    by = c(
      "loss",
      "repetition",
      "seed"
    )
  ) %>%
  mutate(
    Predation = pmax(
      0,
      final_predation - baseline_predation
    ),
    Starvation = pmax(
      0,
      final_starvation - baseline_starvation
    ),
    Crowding = pmax(
      0,
      final_crowding - baseline_crowding
    )
  ) %>%
  select(
    loss,
    repetition,
    seed,
    Predation,
    Starvation,
    Crowding
  ) %>%
  pivot_longer(
    cols = c(
      Predation,
      Starvation,
      Crowding
    ),
    names_to = "cause",
    values_to = "deaths"
  )

mortality_summary <- post_treatment_deaths %>%
  group_by(
    loss,
    cause
  ) %>%
  summarise(
    mean_deaths = mean(deaths),
    median_deaths = median(deaths),
    total_deaths = sum(deaths),
    .groups = "drop"
  ) %>%
  group_by(loss) %>%
  mutate(
    proportion =
      total_deaths / sum(total_deaths)
  ) %>%
  ungroup()

p <- ggplot(
  mortality_summary,
  aes(
    x = factor(loss),
    y = mean_deaths,
    fill = cause
  )
) +
  geom_col(width = 0.82) +
  scale_fill_manual(
    values = mortality_colours
  ) +
  labs(
    title = "Mean post-treatment prey deaths",
    subtitle = paste(
      "Cumulative deaths after tick 2000,",
      "averaged across runs"
    ),
    x = "Habitat loss (%)",
    y = "Mean deaths per run",
    fill = "Cause"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

save_plot(
  "18_post_treatment_deaths",
  p
)

p <- ggplot(
  mortality_summary,
  aes(
    x = factor(loss),
    y = proportion,
    fill = cause
  )
) +
  geom_col(width = 0.82) +
  scale_fill_manual(
    values = mortality_colours
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1)
  ) +
  labs(
    title = "Composition of post-treatment mortality",
    subtitle = paste(
      "Shares use deaths pooled across all runs",
      "within each treatment"
    ),
    x = "Habitat loss (%)",
    y = "Share of recorded prey deaths",
    fill = "Cause"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

save_plot(
  "19_post_treatment_mortality_composition",
  p
)

# ----------------------------------------------------------------------
# 18. Stability diagnostic
# ----------------------------------------------------------------------

stability_summary <- runs %>%
  group_by(loss) %>%
  summarise(
    n = n(),
    events = sum(stable_dynamics),
    .groups = "drop"
  ) %>%
  mutate(
    probability = events / n,
    lower = wilson_lower(events, n),
    upper = wilson_upper(events, n)
  )

p <- ggplot(
  stability_summary,
  aes(
    x = loss,
    y = probability
  )
) +
  geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper
    ),
    fill = "#009E73",
    alpha = 0.20
  ) +
  geom_line(
    colour = "#009E73",
    linewidth = 1
  ) +
  geom_point(
    colour = "#009E73",
    size = 2
  ) +
  scale_x_continuous(
    breaks = seq(0, 100, 10)
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1)
  ) +
  labs(
    title = "Runs meeting the stability criteria",
    subtitle = paste(
      "This is a model diagnostic,",
      "not the primary habitat-loss endpoint"
    ),
    x = "Habitat loss (%)",
    y = "Stable runs"
  )

save_plot(
  "20_stability_diagnostic",
  p
)

# ----------------------------------------------------------------------
# Export derived data
# ----------------------------------------------------------------------

calculated_treatment_summary <- runs %>%
  group_by(loss) %>%
  summarise(
    completed_runs = n(),
    prey_extinction_rate =
      mean(prey_extinct),
    predator_extinction_rate =
      mean(predator_extinct),
    co_persistence_rate =
      mean(survived_to_max_ticks),
    stable_dynamics_rate =
      mean(stable_dynamics),
    mean_actual_habitat_loss =
      mean(actual_habitat_loss_percent),
    median_final_preys =
      median(final_preys),
    median_final_predators =
      median(final_predators),
    median_final_tick =
      median(final_tick),
    .groups = "drop"
  )

write_csv(
  calculated_treatment_summary,
  file.path(
    output_dir,
    "derived_treatment_summary.csv"
  )
)

write_csv(
  paired_response,
  file.path(
    output_dir,
    "derived_paired_responses.csv"
  )
)

write_csv(
  mortality_summary,
  file.path(
    output_dir,
    "derived_mortality_summary.csv"
  )
)

write_csv(
  population_summary,
  file.path(
    output_dir,
    "derived_population_timeseries.csv"
  )
)

# ----------------------------------------------------------------------
# Combined multipage PDF
# ----------------------------------------------------------------------

pdf(
  file.path(
    output_dir,
    "all_habitat_result_plots.pdf"
  ),
  width = 11,
  height = 7.5,
  onefile = TRUE
)

walk(plots, print)

dev.off()

message(
  "Finished. Created ",
  length(plots),
  " figures in: ",
  normalizePath(output_dir)
)
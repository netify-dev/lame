# Helper function to update progress bars in lame
update_lame_progress <- function() {
  # Replace old progress bar code with cli versions
  
  # Example replacements:
  # OLD: pbBurn <- txtProgressBar(min=1,max=burn,style=3)
  # NEW: 
  create_burn_progress <- function(burn) {
    if(burn > 0) {
      cli::cli_progress_bar(
        format = "{cli::pb_spin} Burn-in: {cli::pb_current}/{cli::pb_total} [{cli::pb_bar}] {cli::pb_percent}",
        total = burn,
        clear = FALSE
      )
    }
  }
  
  create_main_progress <- function(nscan) {
    cli::cli_progress_bar(
      format = "{cli::pb_spin} Sampling: {cli::pb_current}/{cli::pb_total} [{cli::pb_bar}] {cli::pb_percent} | ETA: {cli::pb_eta}",
      total = nscan,
      clear = FALSE
    )
  }
}

# Message replacements
show_lame_messages <- function() {
  # Replace cat() with cli functions
  
  # WARNING: row effects are not estimable
  cli::cli_warn("Row effects are not estimable using this procedure")
  
  # WARNING: intercept is not estimable  
  cli::cli_warn("An intercept is not estimable using this procedure")
  
  # WARNING: Random reordering for ties
  cli::cli_warn("Random reordering used to break ties in ranks")
  
  # Starting burn-in
  cli::cli_h2("Starting burn-in period")
  
  # Burn-in complete
  cli::cli_alert_success("Burn-in period complete")
  
  # Model fitting complete
  cli::cli_alert_success("MCMC sampling complete")
}

# Parameter display with cli formatting
show_parameters <- function(s, beta_means, eff_sizes = NULL) {
  cli::cli_div(theme = list(span.param = list(color = "blue")))
  
  param_str <- paste(round(beta_means, 2), collapse = ", ")
  cli::cli_text("Iteration {.val {s}}: beta = [{.param {param_str}}]")
  
  if(!is.null(eff_sizes)) {
    eff_str <- paste(round(eff_sizes), collapse = ", ")
    cli::cli_text("  Effective sizes: [{.emph {eff_str}}]")
  }
  
  cli::cli_end()
}
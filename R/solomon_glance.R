#' @importFrom generics glance
#' @export
glance.solomon_glm <- function(x, ...) {
  eff <- x$effects
  ate <- eff$contrast == "ATE (avg over pretest)"
  tibble::tibble(
    n            = stats::nobs(x$model),
    df_resid     = x$model$df.residual,
    robust       = x$robust,
    ate_est      = eff$estimate[ate],
    ate_conf.low = if (is.null(eff$conf.low)) NA_real_ else eff$conf.low[ate],
    ate_conf.high = if (is.null(eff$conf.high)) NA_real_ else eff$conf.high[ate],
    ate_p        = eff$p.value[ate],
    inter_p      = eff$p.value[eff$contrast == "Pretest x Treatment"],
    pre_p        = eff$p.value[eff$contrast == "Treatment | pretested"],
    un_p         = eff$p.value[eff$contrast == "Treatment | unpretested"]
  )
}

#' @export
glance.solomon_glm <- function(x, ...) {
  eff <- x$effects
  tibble::tibble(
    n = nobs(x$model),
    df_resid = x$model$df.residual,
    ate_est = eff$estimate[eff$contrast=="ATE (avg over pretest)"],
    ate_p   = eff$p.value[eff$contrast=="ATE (avg over pretest)"],
    inter_p = eff$p.value[eff$contrast=="Pretest x Treatment"],
    pre_p   = eff$p.value[eff$contrast=="Treatment | pretested"],
    un_p    = eff$p.value[eff$contrast=="Treatment | unpretested"]
  )
}

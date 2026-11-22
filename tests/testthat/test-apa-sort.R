# Issue #191: report_solomon() and analysis_plan_solomon() order their
# references with .apa_sort(), by the rule of tools/check-references.R, in
# an order that does not depend on the session's collation. (That script
# also checks that .apa_sort() puts the whole reference table in its order.)

test_that(".apa_sort() orders letter by letter, the end of a name first", {
  refs <- c(
    "Browning, A. (2001). Title.",
    "Little, R. J. A. (1990). Title.",
    "Brown, J. (2010). Title.",
    "Little, R. J. (1995). Title.",
    "Smith, A. (2020). Title.",
    "Smith, A. (n.d.). Title.",
    "Smith, A., & Jones, B. (2001). Title.",
    "Dorf, B. (2000). Title.",
    "Dor\u00e9n, A. (2000). Title."
  )
  # "Brown" before "Browning", and "Little, R. J." before "Little, R. J. A.":
  # the end of a name comes before any letter. An accented letter counts as
  # its base letter (Doren before Dorf). No date before a date, and one
  # author before the same author with others.
  expect_identical(.apa_sort(refs), refs[c(3, 1, 9, 8, 4, 2, 6, 5, 7)])
  expect_identical(.apa_sort(rev(refs)), .apa_sort(refs))

  table <- unname(.solomon_reference_text)
  expect_identical(.apa_sort(rev(table)), .apa_sort(table))
  expect_identical(sort(.apa_sort(table), method = "radix"), sort(table, method = "radix"))
})

test_that(".apa_sort() gives one order under the C and Czech collations", {
  old <- Sys.getlocale("LC_COLLATE")
  on.exit(Sys.setlocale("LC_COLLATE", old), add = TRUE)

  table <- unname(.solomon_reference_text)
  small <- c("Hill, B. (2000). Title.", "Chan, A. (2000). Title.",
             "Carpenter, J. (2000). Title.")
  plan_refs <- function() {
    p <- analysis_plan_solomon(sensitization = "equivalence", equivalence_bound = 2)
    refs <- p$text[seq(which(p$text == "## References") + 2L, length(p$text))]
    refs[nzchar(refs)]
  }

  expect_true(nzchar(Sys.setlocale("LC_COLLATE", "C")))
  table_c <- .apa_sort(table)
  plan_c <- plan_refs()
  expect_identical(.apa_sort(small), small[c(3, 2, 1)])

  tried <- c("Czech_Czechia.utf8", "cs_CZ.UTF-8", "cs_CZ.utf8", "Czech_Czechia.1250", "cs_CZ")
  czech <- ""
  for (loc in tried) {
    if (nzchar(suppressWarnings(Sys.setlocale("LC_COLLATE", loc)))) {
      czech <- loc
      break
    }
  }
  if (!nzchar(czech)) {
    skip(paste("no Czech collation can be set on this machine; tried", paste(tried, collapse = ", ")))
  }
  # The test proves something only if this collation puts "ch" after "h".
  if (!identical(order(c("chan", "hill")), c(2L, 1L))) {
    skip(paste0("the collation ", czech, " does not sort 'ch' after 'h' on this machine"))
  }
  expect_identical(.apa_sort(small), small[c(3, 2, 1)])
  expect_identical(.apa_sort(table), table_c)
  expect_identical(plan_refs(), plan_c)
})

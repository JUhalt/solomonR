# Text set in the package's figures (issue #108): captions, subtitles, and
# keys are broken into lines that fit a figure 7 inches wide.

test_that("the line length follows the width and the point size", {

  # An average character width of 0.55 em on a line of 6.8 inches.
  expect_equal(.figure_text_chars(.figure_caption_size), 92L)
  expect_equal(.figure_text_chars(.figure_subtitle_size), 74L)
  expect_equal(.figure_text_chars(12, width = 5), 54L)
})


test_that("text is wrapped between words, and its line breaks are kept", {

  caution <- paste(
    "Teaching and replication aid, not a recommended workflow: the conditional",
    "sequence ending in Test I inflates experiment-wise Type I error",
    "(Sawilowsky et al., 1994)."
  )
  wrapped <- .wrap_figure_text(caution, .figure_caption_size)
  lines <- strsplit(wrapped, "\n", fixed = TRUE)[[1]]
  expect_equal(lines, c(
    "Teaching and replication aid, not a recommended workflow: the conditional sequence ending in",
    "Test I inflates experiment-wise Type I error (Sawilowsky et al., 1994)."
  ))
  expect_true(all(nchar(lines) <= 92L))
  expect_equal(gsub("\n", " ", wrapped, fixed = TRUE), caution)

  # Short text is unchanged, and each element or line starts a new line.
  expect_equal(.wrap_figure_text("First line.\nSecond line.", 9.6), "First line.\nSecond line.")
  expect_equal(.wrap_figure_text(c("First.", "Second."), 9.6), "First.\nSecond.")
  expect_equal(.wrap_figure_text("Extra   spaces  go.", 9.6), "Extra spaces go.")
  expect_null(.wrap_figure_text(NULL, 9.6))
})


test_that("filled text is broken only between its pieces", {

  pieces <- c("Path taken (1988 flow): A -> D -> E -> H -> I",
              "(alpha allocation: method2_conservative (Sawilowsky, 1996))")
  expect_equal(.fill_figure_text(pieces, 12), paste(pieces, collapse = "\n"))
  expect_equal(.fill_figure_text(c("Path taken: A -> B -> C", "(alpha = 0.05)"), 12),
               "Path taken: A -> B -> C (alpha = 0.05)")

  # A piece longer than a line gets a line of its own.
  long <- strrep("x", 100)
  expect_equal(.fill_figure_text(c("a", long, "b"), 12), paste("a", long, "b", sep = "\n"))
  expect_equal(.fill_figure_text(c("x", "y"), 12, sep = "; "), "x; y")
  expect_equal(.fill_figure_text(character(), 12), "")
})

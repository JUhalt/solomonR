# Text set in the package's figures (issue #108).
#
# A ggplot does not know the size it will be drawn at, so captions, subtitles,
# and keys are broken into lines when the figure is made. They are broken for
# a figure 7 inches wide: the narrowest width at which the vignettes show
# these figures, and the default in knitr, R Markdown, and Quarto. pkgdown
# draws its figures 7.29 inches wide. A line of text has 6.8 inches: the 7
# inches less the plot margins of 5.5 points on each side.
.figure_text_width <- 6.8

# Point sizes of the subtitle and the caption. The figures use a 12-point base
# size; ggplot2 sets the subtitle at the base size and the caption at 0.8 of
# it.
.figure_subtitle_size <- 12
.figure_caption_size <- 9.6

# Characters that fit on a line of figure text `width` inches wide at `size`
# points. The count assumes an average character width of 0.55 em, so that a
# line fits in any of the usual default fonts: on the package's captions the
# average is about 0.45 em in Arial and Helvetica (the defaults on Windows and
# macOS, and the metrics of the pdf() device) and about 0.51 em in DejaVu
# Sans (the default on most Linux systems).
.figure_text_chars <- function(size, width = .figure_text_width) {
  as.integer(floor(width * 72 / (0.55 * size)))
}

# Sets `pieces` on lines that fit `width` inches at `size` points, joined by
# `sep` and broken only between pieces. A piece longer than a line gets a
# line of its own.
.fill_figure_text <- function(pieces, size, width = .figure_text_width, sep = " ") {
  chars <- .figure_text_chars(size, width)
  lines <- character()
  current <- NULL
  for (piece in pieces) {
    if (is.null(current)) {
      current <- piece
    } else if (nchar(current) + nchar(sep) + nchar(piece) <= chars) {
      current <- paste0(current, sep, piece)
    } else {
      lines <- c(lines, current)
      current <- piece
    }
  }
  paste(c(lines, current), collapse = "\n")
}

# Wraps figure text, such as a caption or a subtitle, between words so that
# each line fits `width` inches at `size` points. Line breaks already in the
# text are kept.
.wrap_figure_text <- function(text, size, width = .figure_text_width) {
  if (is.null(text)) {
    return(NULL)
  }
  paragraphs <- strsplit(paste(text, collapse = "\n"), "\n", fixed = TRUE)[[1]]
  wrapped <- vapply(paragraphs, function(paragraph) {
    words <- strsplit(paragraph, " ", fixed = TRUE)[[1]]
    .fill_figure_text(words[nzchar(words)], size, width)
  }, character(1), USE.NAMES = FALSE)
  paste(wrapped, collapse = "\n")
}

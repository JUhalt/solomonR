# Check that every reference entry in the package matches the canonical
# APA 7 list in vignettes/articles/references.Rmd, and that each reference
# list is in APA order. Run from the package root:
#
#   Rscript tools/check-references.R
#
# Reference entries are read from roxygen `@references` blocks in R/ and from
# "References" sections of the vignettes, articles, and README.Rmd. Entries
# are compared after collapsing whitespace, so line wrapping does not matter.
# Exits with status 1 if any entry differs from the canonical list.

# A reference starts with an author, not a blockquote marker or list bullet.
ref_pattern <- "^[^()>*-][^()]*? \\(\\d{4}[a-z]?(, [^)]*)?\\)\\."

paragraphs <- function(lines) {
  groups <- cumsum(!nzchar(trimws(lines)))
  paras <- split(trimws(lines), groups)
  paras <- vapply(paras, function(p) paste(p[nzchar(p)], collapse = " "), "")
  unname(paras[nzchar(paras)])
}

sort_key <- function(entry) {
  authors <- sub(" \\(\\d{4}.*$", "", entry)
  year <- sub("^[^()]+? \\((\\d{4}[a-z]?).*$", "\\1", entry)
  paste(gsub("[^a-z]", "", tolower(gsub("&", "", authors, fixed = TRUE))), year)
}

canonical_file <- file.path("vignettes", "articles", "references.Rmd")
canonical <- paragraphs(readLines(canonical_file, encoding = "UTF-8"))
canonical <- canonical[grepl(ref_pattern, canonical, perl = TRUE)]

problems <- character()
note <- function(where, what) problems <<- c(problems, sprintf("%s: %s", where, what))

check_list <- function(entries, where) {
  entries <- entries[grepl(ref_pattern, entries, perl = TRUE)]
  for (e in setdiff(entries, canonical)) {
    note(where, paste("not in the canonical list:", substr(e, 1, 90)))
  }
  if (anyDuplicated(entries)) note(where, "duplicate entry")
  keys <- vapply(entries, sort_key, "")
  if (is.unsorted(keys)) note(where, "entries are not in APA order")
  length(entries)
}

n_checked <- 0L

for (path in list.files("R", pattern = "[.]R$", full.names = TRUE)) {
  lines <- readLines(path, encoding = "UTF-8")
  starts <- grep("^#'\\s*@references\\s*$", lines)
  for (s in starts) {
    end <- s
    while (end < length(lines) && grepl("^#'", lines[end + 1]) &&
           !grepl("^#'\\s*@\\w", lines[end + 1])) {
      end <- end + 1
    }
    block <- sub("^#' ?", "", lines[seq_len(end - s) + s])
    n_checked <- n_checked + check_list(paragraphs(block), sprintf("%s:%d", path, s))
  }
}

md_files <- c(
  list.files("vignettes", pattern = "[.]Rmd$", full.names = TRUE),
  setdiff(list.files(file.path("vignettes", "articles"), pattern = "[.]Rmd$", full.names = TRUE),
          canonical_file),
  "README.Rmd"
)
for (path in md_files) {
  lines <- readLines(path, encoding = "UTF-8")
  starts <- grep("^#{2,3} References\\s*$", lines)
  for (s in starts) {
    end <- s
    while (end < length(lines) && !grepl("^(#{1,3} |```|-{3,}\\s*$|<)", lines[end + 1])) {
      end <- end + 1
    }
    block <- lines[seq_len(end - s) + s]
    n_checked <- n_checked + check_list(paragraphs(block), sprintf("%s:%d", path, s))
  }
}

canonical_lines <- readLines(canonical_file, encoding = "UTF-8")
sections <- split(canonical_lines, cumsum(grepl("^## ", canonical_lines)))
for (section in sections[-1]) {
  invisible(check_list(paragraphs(section[-1]), paste(canonical_file, section[1])))
}
if (anyDuplicated(canonical)) note(canonical_file, "duplicate entry")

cat(sprintf("Checked %d reference entries against %d canonical references.\n",
            n_checked, length(canonical)))
if (length(problems)) {
  cat(problems, sep = "\n")
  quit(status = 1)
}
cat("All reference entries match the canonical list.\n")

# Check that every reference entry in the package matches the canonical
# APA 7 list in vignettes/articles/references.Rmd, and that each reference
# list is in APA order. Run from the package root:
#
#   Rscript tools/check-references.R
#
# Reference entries are read from roxygen `@references` blocks in R/ and from
# "References" sections of the vignettes, articles, and README.Rmd. Entries
# are compared after collapsing whitespace, so line wrapping does not matter,
# and after removing the angle brackets of explicit Markdown links such as
# <https://doi.org/...>, which DOIs containing "::" need in order to render as
# one link. It also checks that every work in the Solomon-design section of
# the canonical list appears in the coverage article. Exits with status 1 if
# any entry differs from the canonical list or lacks a coverage entry.

# A reference starts with an author, not a blockquote marker or list bullet.
ref_pattern <- "^[^()>*-][^()]*? \\(\\d{4}[a-z]?(, [^)]*)?\\)\\."

paragraphs <- function(lines) {
  groups <- cumsum(!nzchar(trimws(lines)))
  paras <- split(trimws(lines), groups)
  paras <- vapply(paras, function(p) paste(p[nzchar(p)], collapse = " "), "")
  paras <- gsub("<(https?://[^>[:space:]]+)>", "\\1", paras)
  unname(paras[nzchar(paras)])
}

# APA 7 orders entries by the first author's surname and initials, then by
# the other authors and the year, so "Little, R. J." precedes "Little, R. J.
# A."; the first author is therefore compared on its own.
sort_key <- function(entry) {
  authors <- sub(" \\(\\d{4}.*$", "", entry)
  year <- sub("^[^()]+? \\((\\d{4}[a-z]?).*$", "\\1", entry)
  parts <- strsplit(authors, ", ", fixed = TRUE)[[1]]
  first <- paste(parts[seq_len(min(2L, length(parts)))], collapse = " ")
  rest <- paste(parts[-seq_len(min(2L, length(parts)))], collapse = " ")
  letters_only <- function(x) gsub("[^a-z]", "", tolower(gsub("&", "", x, fixed = TRUE)))
  paste(letters_only(first), letters_only(rest), year)
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
    # Roxygen markdown needs literal brackets escaped, as in "\[Paper
    # presentation\]"; compare the text as it renders.
    block <- gsub("\\\\([][])", "\\1", block)
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
    while (end < length(lines) &&
           !grepl("^(#{1,3} |```|-{3,}\\s*$|<(?!https?://))", lines[end + 1], perl = TRUE)) {
      end <- end + 1
    }
    block <- lines[seq_len(end - s) + s]
    n_checked <- n_checked + check_list(paragraphs(block), sprintf("%s:%d", path, s))
  }
}

# The references report_solomon() cites (R/solomon_report_refs.R).
registry <- new.env()
sys.source(file.path("R", "solomon_report_refs.R"), envir = registry)
for (key in names(registry$.solomon_reference_text)) {
  entry <- enc2utf8(registry$.solomon_reference_text[[key]])
  if (!entry %in% canonical) {
    note("R/solomon_report_refs.R", paste("not in the canonical list:", key))
  }
}
n_checked <- n_checked + length(registry$.solomon_reference_text)

canonical_lines <- readLines(canonical_file, encoding = "UTF-8")
sections <- split(canonical_lines, cumsum(grepl("^## ", canonical_lines)))
for (section in sections[-1]) {
  invisible(check_list(paragraphs(section[-1]), paste(canonical_file, section[1])))
}
if (anyDuplicated(canonical)) note(canonical_file, "duplicate entry")

# Every work in the Solomon-design section of the bibliography must have an
# entry in the coverage article, whose reference list names each one (#80).
coverage_file <- file.path("vignettes", "articles", "coverage.Rmd")
is_solomon <- vapply(sections, function(s) grepl("^## The Solomon design literature", s[1]), NA)
solomon_refs <- paragraphs(sections[[which(is_solomon)]][-1])
solomon_refs <- solomon_refs[grepl(ref_pattern, solomon_refs, perl = TRUE)]
coverage_lines <- readLines(coverage_file, encoding = "UTF-8")
coverage_start <- grep("^## References\\s*$", coverage_lines)
coverage_refs <- paragraphs(coverage_lines[seq(coverage_start + 1L, length(coverage_lines))])
for (e in setdiff(solomon_refs, coverage_refs)) {
  note(coverage_file, paste("no coverage entry for:", substr(e, 1, 90)))
}

cat(sprintf("Checked %d reference entries against %d canonical references.\n",
            n_checked, length(canonical)))
if (length(problems)) {
  cat(problems, sep = "\n")
  quit(status = 1)
}
cat("All reference entries match the canonical list.\n")

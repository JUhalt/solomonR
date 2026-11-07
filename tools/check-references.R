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

# A reference starts with an author, not a blockquote marker or list bullet,
# and gives a year or "n.d." (no date) in parentheses.
ref_pattern <- "^[^()>*-][^()]*? \\((\\d{4}[a-z]?|n\\.d\\.)(, [^)]*)?\\)\\."

paragraphs <- function(lines) {
  groups <- cumsum(!nzchar(trimws(lines)))
  paras <- split(trimws(lines), groups)
  paras <- vapply(paras, function(p) paste(p[nzchar(p)], collapse = " "), "")
  paras <- gsub("<(https?://[^>[:space:]]+)>", "\\1", paras)
  unname(paras[nzchar(paras)])
}

# APA 7 orders entries letter by letter by the first author's surname and
# initials, then by the other authors and the year. Spaces and punctuation
# within a name are ignored, but the end of a surname or of a set of
# initials comes before any letter (nothing precedes something), so
# "Brown, J." precedes "Browning, A.", and "Little, R. J." precedes
# "Little, R. J. A.". Each name and set of initials is therefore a field of
# its own, ended by a tab. Accented Latin letters count as their base
# letters ("Ethier" with an acute accent sorts under E), from a fixed table
# rather than the platform's iconv; the table is written in escapes, so the
# file stays ASCII and reads the same in any locale. A work with no date
# ("n.d.") precedes the same authors' dated works. The keys are compared
# byte by byte (method = "radix"), so the order does not depend on the
# locale.
accented <- paste0(
  "\u00e0\u00e1\u00e2\u00e3\u00e4\u00e5\u0101\u0103\u0105\u01ce",       # a
  "\u00e7\u0107\u0109\u010b\u010d",                                     # c
  "\u010f\u0111\u00f0",                                                 # d
  "\u00e8\u00e9\u00ea\u00eb\u0113\u0115\u0117\u0119\u011b",             # e
  "\u011d\u011f\u0121\u0123",                                           # g
  "\u0125\u0127",                                                       # h
  "\u00ec\u00ed\u00ee\u00ef\u0129\u012b\u012d\u012f\u0131\u01d0",       # i
  "\u0135",                                                             # j
  "\u0137",                                                             # k
  "\u013a\u013c\u013e\u0140\u0142",                                     # l
  "\u00f1\u0144\u0146\u0148",                                           # n
  "\u00f2\u00f3\u00f4\u00f5\u00f6\u00f8\u014d\u014f\u0151\u01d2",       # o
  "\u0155\u0157\u0159",                                                 # r
  "\u015b\u015d\u015f\u0161\u0219",                                     # s
  "\u0163\u0165\u0167\u021b",                                           # t
  "\u00f9\u00fa\u00fb\u00fc\u0169\u016b\u016d\u016f\u0171\u0173\u01d4", # u
  "\u0175",                                                             # w
  "\u00fd\u00ff\u0177",                                                 # y
  "\u017a\u017c\u017e"                                                  # z
)
unaccented <- paste(rep(c("a", "c", "d", "e", "g", "h", "i", "j", "k", "l", "n", "o",
                          "r", "s", "t", "u", "w", "y", "z"),
                        c(10, 5, 3, 9, 4, 2, 10, 1, 1, 5, 4, 10, 3, 5, 4, 11, 1, 3, 3)),
                    collapse = "")
stopifnot(nchar(accented) == nchar(unaccented))
letters_only <- function(x) {
  x <- tolower(gsub("&", "", x, fixed = TRUE))
  for (pair in list(c("\u00df", "ss"), c("\u00e6", "ae"), c("\u0153", "oe"), c("\u00fe", "th"))) {
    x <- gsub(pair[1], pair[2], x, fixed = TRUE)
  }
  gsub("[^a-z]", "", chartr(accented, unaccented, x))
}
sort_key <- function(entry) {
  authors <- sub(" \\((\\d{4}|n\\.d\\.).*$", "", entry)
  year <- sub("^[^()]+? \\((\\d{4}[a-z]?|n\\.d\\.).*$", "\\1", entry)
  if (identical(year, "n.d.")) year <- "0000"
  parts <- strsplit(authors, ", ", fixed = TRUE)[[1]]
  paste(c(letters_only(parts), year), collapse = "\t")
}
in_order <- function(keys) identical(order(keys, method = "radix"), seq_along(keys))

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
  if (!in_order(keys)) note(where, "entries are not in APA order")
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

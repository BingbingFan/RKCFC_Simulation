# Caption = table name only; everything else goes to the note under the table.
# Author's instruction, 2026-10-01. The rules live in R/table_captions.txt (plain
# text, so LaTeX and regular expressions need no R-string escaping); this file
# reads them and is sourced by R/make_tables.R and R/make_realdata_tables.R,
# whose emit() calls split_caption() before writing a table.

read_caption_spec <- function(path = file.path("R", "table_captions.txt")) {
  x <- readLines(path, encoding = "UTF-8", warn = FALSE)
  x <- x[!grepl("^[[:space:]]*(#|$)", x)]
  spec <- list(); cur <- NULL
  for (ln in x) {
    if (startsWith(ln, "@@ ")) {
      cur <- sub("^@@ ", "", ln)
      spec[[cur]] <- list(title = NA_character_, match = NA_character_, lead = "")
    } else if (is.null(cur)) {
      stop("R/table_captions.txt: a rule line appears before the first @@ block: ", ln)
    } else if (startsWith(ln, "TITLE ")) {
      spec[[cur]]$title <- sub("^TITLE ", "", ln)
    } else if (startsWith(ln, "MATCH ")) {
      spec[[cur]]$match <- sub("^MATCH ", "", ln)
    } else if (startsWith(ln, "LEAD")) {
      spec[[cur]]$lead <- sub("^LEAD ?", "", ln)
    } else stop("R/table_captions.txt: unrecognized line in block ", cur, ": ", ln)
  }
  for (id in names(spec))
    if (is.na(spec[[id]]$title) || is.na(spec[[id]]$match))
      stop("R/table_captions.txt: block ", id, " needs both TITLE and MATCH")
  spec
}
CAPTION_SPEC <- read_caption_spec()

# Returns list(caption = <table name>, notes = <moved text, then the old notes>).
split_caption <- function(id, caption, notes) {
  sp <- CAPTION_SPEC[[id]]
  if (is.null(sp)) stop("no caption rule for ", id, " in R/table_captions.txt")
  m <- regexpr(sp$match, caption, perl = TRUE)
  if (m != 1)
    stop(id, ": the caption no longer starts the way R/table_captions.txt expects.\n",
         "  caption: ", caption, "\n  MATCH:   ", sp$match)
  fill <- function(s) {
    cs <- attr(m, "capture.start"); cl <- attr(m, "capture.length")
    if (!is.null(cs) && length(cs))
      for (k in seq_along(cs))
        s <- gsub(paste0("{{", k, "}}"), substr(caption, cs[k], cs[k] + cl[k] - 1), s, fixed = TRUE)
    if (grepl("{{", s, fixed = TRUE)) stop(id, ": unfilled placeholder in ", s)
    s
  }
  rest <- substr(caption, attr(m, "match.length") + 1L, nchar(caption))
  lead <- trimws(paste0(fill(sp$lead), rest))
  old  <- if (is.null(notes)) character(0) else paste(notes, collapse = " ")
  list(caption = fill(sp$title),
       notes = paste(c(if (nzchar(lead)) lead, old), collapse = " "))
}

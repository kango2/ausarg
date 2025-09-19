# install.packages(c("DBI","RMariaDB"))
library(DBI)
library(RMariaDB)

#' Fetch rows by RefCode from a MySQL/MariaDB table
#'
#' @param refcodes  character vector of RefCode values to retrieve
#' @param db        database name (default "Team_AusARG")
#' @param table     table name    (default "G_genomics")
#' @param host      hostname or IP (default "georges.biomatix.org")
#' @param user      username (default "anon")
#' @param password  password (default "dragon")
#' @param port      TCP port (default 3306)
#' @param ssl.ca    optional path to CA cert, if required
#' @param chunk_size query in batches of this many refcodes (default 1000)
#' @return data.frame with the same column names as in the table
#' @example
#' refvec <- c("123","121")
#' df <- gn.pull.records(refvec)
#' df

gn.pull.records <- function(refcodes,
                             db = "Team_AusARG",
                             table = "G_genomics",
                             host = "georges.biomatix.org",
                             user = "anon",
                             password = "dragon",
                             port = 3306,
                             ssl.ca = NULL,
                             chunk_size = 1000) {
  
  # basic checks / cleaning
  refcodes <- unique(na.omit(refcodes))
  if (length(refcodes) == 0) {
    con <- dbConnect(MariaDB(), host = host, user = user, password = password,
                     dbname = db, port = port, ssl.ca = ssl.ca)
    on.exit(try(dbDisconnect(con), silent = TRUE), add = TRUE)
    return(dbGetQuery(con, sprintf("SELECT * FROM `%s`.`%s` WHERE 1=0", db, table)))
  }
  
  con <- dbConnect(MariaDB(), host = host, user = user, password = password,
                   dbname = db, port = port, ssl.ca = ssl.ca)
  on.exit(try(dbDisconnect(con), silent = TRUE), add = TRUE)
  
  # helper: one chunk query using safe quoting
  run_chunk <- function(vals) {
    quoted <- dbQuoteString(con, vals)
    sql <- sprintf("SELECT * FROM `%s`.`%s` WHERE `RefCode` IN (%s)",
                   db, table, paste(quoted, collapse = ","))
    dbGetQuery(con, sql)
  }
  
  # split into chunks (avoids overly long IN() lists)
  idx <- ceiling(seq_along(refcodes) / chunk_size)
  chunks <- split(refcodes, idx)
  
  out_list <- lapply(chunks, run_chunk)
  df <- do.call(rbind, out_list)
  
  return(df)
}

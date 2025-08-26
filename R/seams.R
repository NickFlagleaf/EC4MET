# these are non exported, to help with mocking

.download_to <- function(url, dest, method = "libcurl", quiet = TRUE, mode = "wb") {
	  mapply(function(u, d) {
			     utils::download.file(url = u, destfile = d, method = method, quiet = quiet, mode = mode)
			         return(d)
			       }, url, dest, SIMPLIFY = TRUE, USE.NAMES = FALSE)
}

# this is just a helper to make writing tests easier
.copy_fixture <- function(src, dest) {
  mapply(function(s, d) { file.copy(s, d, overwrite = TRUE); d },
         src, dest, SIMPLIFY = TRUE, USE.NAMES = FALSE)
}


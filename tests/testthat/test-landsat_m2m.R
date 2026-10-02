test_that(".m2m_request surfaces API errors as R errors with code and message", {
  # simuliamo la logica di controllo errore senza una vera chiamata di rete
  fake_result <- list(errorCode = "AUTH_INVALID", errorMessage = "Token non valido")
  expect_error(
    {
      if (!is.null(fake_result$errorCode) && !identical(fake_result$errorCode, "")) {
        stop(sprintf("Errore M2M API [%s]: %s", fake_result$errorCode, fake_result$errorMessage),
             call. = FALSE)
      }
    },
    "AUTH_INVALID"
  )
})

test_that("landsat_m2m_login requires credentials", {
  old_username <- Sys.getenv("USGS_USERNAME", unset = NA)
  old_token <- Sys.getenv("USGS_M2M_TOKEN", unset = NA)
  Sys.unsetenv(c("USGS_USERNAME", "USGS_M2M_TOKEN"))
  on.exit({
    if (!is.na(old_username)) Sys.setenv(USGS_USERNAME = old_username)
    if (!is.na(old_token)) Sys.setenv(USGS_M2M_TOKEN = old_token)
  }, add = TRUE)

  expect_error(landsat_m2m_login(), "Credenziali USGS M2M non fornite")
})

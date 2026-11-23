# Tests for inst/shiny/server.R, driven with shiny::testServer(). An error
# in an observer ends a real session, and closes the test session, so
# session$isClosed() shows whether an observer failed.

test_that("combining one table says two or more are needed (#277)", {
  env <- shiny_server_env()
  shiny::testServer(env$server, {
    vals$musica <- dbs_musica
    session$setInputs(combined_table_name = "comb", combine_tables = "SBS96")
    session$setInputs(Combine = 1)
    expect_false(session$isClosed())
    expect_error(output$combine_warning, "two or more tables")
  })
})

test_that("comparing with no threshold asks for one (#277)", {
  env <- shiny_server_env()
  shiny::testServer(env$server, {
    session$setInputs(select_result_a = "result", select_model_a = "model",
                      Threshold = "")
    session$setInputs(compare_results = 1)
    expect_false(session$isClosed())
    expect_error(output$compare_validate, "similarity threshold")
  })
})

test_that("comparing before the model input exists doesn't fail (#278)", {
  env <- shiny_server_env()
  shiny::testServer(env$server, {
    session$setInputs(select_result_a = "result", Threshold = ".5")
    session$setInputs(compare_results = 1)
    expect_false(session$isClosed())
    expect_error(output$compare_validate, "model ID")
  })
})

test_that("predicting before its inputs exist doesn't fail (#278)", {
  env <- shiny_server_env()
  shiny::testServer(env$server, {
    vals$musica <- musica_sbs96_tiny
    session$setInputs(predict_sigs = 1)
    expect_false(session$isClosed())
    expect_error(output$predict_warning, "result list name")
  })
})

test_that("discovering with a cleared number box doesn't fail (#278)", {
  env <- shiny_server_env()
  shiny::testServer(env$server, {
    vals$musica <- musica_sbs96_tiny
    session$setInputs(discover_result_name = "result",
                      discover_model_name = "model",
                      number_of_signatures = NA, n_start = 1,
                      select_discover_table = "SBS96")
    session$setInputs(discover_signatures = 1)
    expect_false(session$isClosed())
    expect_error(output$discover_warning, "number of expected signatures")
  })
})

test_that("group and color by annotation work before the other exists (#278)", {
  env <- shiny_server_env()
  alerts <- record_shinyalerts()
  shiny::testServer(env$server, {
    vals$musica <- musica_sbs96_tiny
    session$setInputs(group1 = "annotation")
    expect_false(session$isClosed())
  })
  shiny::testServer(env$server, {
    vals$musica <- musica_sbs96_tiny
    session$setInputs(color = "annotation")
    expect_false(session$isClosed())
  })
  # musica_sbs96_tiny has no annotations, so each says so
  texts <- vapply(alerts(), `[[`, "", "text")
  expect_length(texts, 2)
  expect_match(texts, "Annotation not found")
})

test_that("a failed combine doesn't say the table was created (#279)", {
  notes <- character()
  local_mocked_bindings(
    showNotification = function(ui, ...) {
      notes <<- c(notes, as.character(ui))
      invisible(NULL)
    },
    .package = "shiny"
  )
  alerts <- record_shinyalerts()
  env <- shiny_server_env()
  shiny::testServer(env$server, {
    vals$musica <- dbs_musica
    # "SBS96" exists already, and combine_count_tables() won't overwrite it
    session$setInputs(combined_table_name = "SBS96",
                      combine_tables = c("SBS96", "DBS78"))
    session$setInputs(Combine = 1)
    expect_false(session$isClosed())
  })
  expect_length(alerts(), 1)
  expect_false("Table created." %in% notes)
})

test_that("a successful combine says the table was created (#279)", {
  notes <- character()
  local_mocked_bindings(
    showNotification = function(ui, ...) {
      notes <<- c(notes, as.character(ui))
      invisible(NULL)
    },
    .package = "shiny"
  )
  env <- shiny_server_env()
  shiny::testServer(env$server, {
    vals$musica <- dbs_musica
    session$setInputs(combined_table_name = "comb",
                      combine_tables = c("SBS96", "DBS78"))
    session$setInputs(Combine = 1)
    expect_true("comb" %in% names(tables(vals$musica)))
  })
  expect_true("Table created." %in% notes)
})

test_that("tab-navigation alerts are errors with well-formed text (#280)", {
  alerts <- record_shinyalerts()
  env <- shiny_server_env()
  shiny::testServer(env$server, {
    # Nothing uploaded
    for (tab in c("musica", "annotations", "download", "tables", "discover",
                  "visualization")) {
      session$setInputs(menu = tab)
    }
    # Variants uploaded, but no musica object created
    vals$var <- variants(musica_sbs96_tiny)
    for (tab in c("tables", "predict", "compare")) {
      session$setInputs(menu = tab)
    }
    # musica object created, but no count tables
    musica_no_tables <- musica_sbs96_tiny
    tables(musica_no_tables) <- list()
    vals$musica <- musica_no_tables
    for (tab in c("discover", "heatmap")) {
      session$setInputs(menu = tab)
    }
    expect_false(session$isClosed())
  })
  shown <- alerts()
  expect_length(shown, 11)
  for (alert in shown) {
    expect_identical(alert$type, "error")
    expect_no_match(alert$text, "muscia")
    # Each sentence and clause is followed by a space
    expect_no_match(alert$text, "[.,]\\S")
    expect_no_match(alert$text, "error$")
    # Quotes come in pairs
    expect_true(lengths(regmatches(alert$text,
                                   gregexpr("\"", alert$text))) %% 2 == 0)
  }
})

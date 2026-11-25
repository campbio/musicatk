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

test_that("comparing with a threshold outside 0 to 1 asks for one", {
  env <- shiny_server_env()
  shiny::testServer(env$server, {
    session$setInputs(select_result_a = "result", select_model_a = "model")
    for (threshold in c("abc", "5")) {
      session$setInputs(Threshold = threshold)
      session$setInputs(compare_results = threshold)
      expect_error(output$compare_validate, "similarity threshold")
    }
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

test_that("discovering before its inputs exist doesn't fail (#278)", {
  env <- shiny_server_env()
  shiny::testServer(env$server, {
    vals$musica <- musica_sbs96_tiny
    session$setInputs(discover_signatures = 1)
    expect_false(session$isClosed())
    expect_error(output$discover_warning, "result list name")
  })
})

test_that("discovery accepts 10 or more signatures (#278)", {
  env <- shiny_server_env()
  shiny::testServer(env$server, {
    vals$musica <- musica_sbs96_tiny
    # The model ID is missing, so the check stops before discovery runs.
    # number_of_signatures is text, and "10" >= 2 is FALSE as strings.
    session$setInputs(discover_result_name = "result",
                      discover_model_name = "",
                      number_of_signatures = "10", n_start = "1",
                      select_discover_table = "SBS96")
    session$setInputs(discover_signatures = 1)
    msg <- tryCatch(output$discover_warning,
                    error = function(e) conditionMessage(e))
    expect_match(msg, "model ID")
    expect_no_match(msg, "2 or more signatures")
  })
})

test_that("discovery asks for random starts that are a number of 1 or more", {
  env <- shiny_server_env()
  shiny::testServer(env$server, {
    vals$musica <- musica_sbs96_tiny
    session$setInputs(discover_result_name = "result",
                      discover_model_name = "model",
                      number_of_signatures = "2",
                      select_discover_table = "SBS96")
    for (n_start in c("ten", "0")) {
      session$setInputs(n_start = n_start)
      session$setInputs(discover_signatures = n_start)
      expect_false(session$isClosed())
      expect_error(output$discover_warning, "random starts")
    }
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

# Sets the Discover inputs for a 2-signature lda run with model ID 'model'
set_discover_inputs <- function(session, model) {
  session$setInputs(discover_result_name = "result",
                    discover_model_name = model, Method = "lda",
                    number_of_signatures = "2", n_start = "1",
                    select_discover_table = "SBS96")
}

# Sets the Predict inputs for an lda prediction with model ID 'model', using
# the first two signatures of 'sigs'
set_predict_inputs <- function(session, model, sigs) {
  session$setInputs(predict_result_name = "result",
                    predict_model_name = model, predict_table = "SBS96",
                    predict_algorithm = "lda",
                    pred_sigs = colnames(signatures(sigs))[1:2])
}

test_that("discovery with a model ID already used asks for another (#282)", {
  env <- shiny_server_env()
  shiny::testServer(env$server, {
    vals$musica <- discover_signatures(musica_sbs96_tiny, "SBS96", 2,
                                       nstart = 1, model_id = "m1",
                                       make_copy = TRUE)
    vals$res_names <- "result/SBS96/m1"
    set_discover_inputs(session, "m1")
    session$setInputs(discover_signatures = 1)
    expect_false(session$isClosed())
    expect_error(output$discover_warning, "already exists")
    expect_identical(names(get_modality(vals$musica, "result", "SBS96")),
                     "m1")
    expect_identical(vals$res_names, "result/SBS96/m1")
  })
})

test_that("prediction with a model ID already used asks for another (#282)", {
  env <- shiny_server_env()
  shiny::testServer(env$server, {
    vals$musica <- predict_exposure(musica_sbs96_tiny, "SBS96",
                                    cosmic_v2_sigs, algorithm = "lda",
                                    model_id = "p1", make_copy = TRUE)
    vals$res_names <- "result/SBS96/p1"
    vals$p_res <- cosmic_v2_sigs
    set_predict_inputs(session, "p1", cosmic_v2_sigs)
    session$setInputs(predict_sigs = 1)
    expect_false(session$isClosed())
    expect_error(output$predict_warning, "already exists")
    expect_identical(names(get_modality(vals$musica, "result", "SBS96")),
                     "p1")
    expect_identical(vals$res_names, "result/SBS96/p1")
  })
})

test_that("a failed prediction shows the error and stores nothing (#283)", {
  notes <- record_notifications()
  alerts <- record_shinyalerts()
  env <- shiny_server_env()
  shiny::testServer(env$server, {
    vals$musica <- musica_sbs96_tiny
    # DBS signatures don't fit the SBS96 table
    vals$p_res <- cosmic_v3_dbs_sigs
    set_predict_inputs(session, "p1", cosmic_v3_dbs_sigs)
    session$setInputs(predict_sigs = 1)
    expect_false(session$isClosed())
    expect_null(vals$res_names)
    expect_length(result_list(vals$musica), 0)
  })
  expect_length(alerts(), 1)
  expect_match(alerts()[[1]]$text, "motifs of 'signature_res' don't match")
  expect_length(notes(), 0)
})

test_that("a successful discovery is recorded and announced", {
  notes <- record_notifications()
  env <- shiny_server_env()
  shiny::testServer(env$server, {
    vals$musica <- musica_sbs96_tiny
    set_discover_inputs(session, "m1")
    session$setInputs(discover_signatures = 1)
    expect_identical(vals$res_names, "result/SBS96/m1")
    expect_identical(vals$res_names_full$result$SBS96, "m1")
    expect_identical(names(get_modality(vals$musica, "result", "SBS96")),
                     "m1")
  })
  expect_match(notes(), "Discovery, m1, completed")
})

test_that("the model ID check uses the current musica object", {
  notes <- record_notifications()
  env <- shiny_server_env()
  shiny::testServer(env$server, {
    # IDs left from an earlier object, which the current one doesn't have
    vals$res_names <- "result/SBS96/m1"
    vals$musica <- musica_sbs96_tiny
    set_discover_inputs(session, "m1")
    session$setInputs(discover_signatures = 1)
    expect_false(session$isClosed())
    expect_identical(names(get_modality(vals$musica, "result", "SBS96")),
                     "m1")
  })
  expect_match(notes(), "Discovery, m1, completed")

  alerts <- record_shinyalerts()
  shiny::testServer(env$server, {
    # A model the object has, though the app's list doesn't name it
    vals$musica <- predict_exposure(musica_sbs96_tiny, "SBS96",
                                    cosmic_v2_sigs, algorithm = "lda",
                                    model_id = "p1", make_copy = TRUE)
    vals$p_res <- cosmic_v2_sigs
    set_predict_inputs(session, "p1", cosmic_v2_sigs)
    session$setInputs(predict_sigs = 1)
    expect_false(session$isClosed())
    expect_error(output$predict_warning, "already exists")
  })
  expect_length(alerts(), 0)
})

test_that("uploaded annotations are matched to samples by name (#211)", {
  notes <- record_notifications()
  alerts <- record_shinyalerts()
  annot <- read.table(system.file("extdata", "sample_annotations.txt",
                                  package = "musicatk"),
                      sep = "\t", header = TRUE)
  env <- shiny_server_env()
  shiny::testServer(env$server, {
    vals$musica <- musica_sbs96
    vals$annotations <- annot
    session$setInputs(annot_sample_column = "Sample_Names")
    session$setInputs(add_annotation = 1)
    expect_false(session$isClosed())
    samples <- sample_names(musica_sbs96)
    expect_identical(sample_names(vals$musica), samples)
    expected <- annot$Tumor_Subtypes[match(samples, annot$Sample_Names)]
    expect_identical(samp_annot(vals$musica)$Tumor_Subtypes, expected)
  })
  expect_length(alerts(), 0)
  expect_identical(notes(), "Annotations have been added")
})

test_that("uploaded annotations keep number-like sample IDs (#287)", {
  notes <- record_notifications()
  alerts <- record_shinyalerts()
  counts <- get_count_table(tables(musica_sbs96)$SBS96)[, 1:3]
  colnames(counts) <- c("001", "002", "100000")
  musica <- create_musica_from_counts(counts, "SBS96")
  annot_file <- withr::local_tempfile(fileext = ".txt")
  # Spaces after the commas, which number-like IDs read as numbers lost
  writeLines(c("group, id, score", "A, 100000, 3", "B, 001, 1.5",
               "C, 002, 2"), annot_file)
  env <- shiny_server_env()
  shiny::testServer(env$server, {
    vals$musica <- musica
    session$setInputs(annotation_header = TRUE, annotation_delimiter = ",",
                      annotations_file = list(name = "annot.txt",
                                              datapath = annot_file))
    output$annotations
    session$setInputs(annot_sample_column = "id")
    session$setInputs(add_annotation = 1)
    expect_false(session$isClosed())
    expect_identical(samp_annot(vals$musica)$group, c("B", "C", "A"))
    expect_identical(samp_annot(vals$musica)$score, c(1.5, 2, 3))
  })
  expect_length(alerts(), 0)
  expect_identical(notes(), "Annotations have been added")
})

test_that("annotation upload shows missing samples in a notification", {
  notes <- record_notifications()
  alerts <- record_shinyalerts()
  samples <- sample_names(musica_sbs96)
  env <- shiny_server_env()
  shiny::testServer(env$server, {
    vals$musica <- musica_sbs96
    vals$annotations <- data.frame(id = samples[-1], group = "a")
    session$setInputs(annot_sample_column = "id")
    session$setInputs(add_annotation = 1)
    expect_false(session$isClosed())
    expect_identical(samp_annot(vals$musica)$group,
                     c(NA, rep("a", length(samples) - 1)))
  })
  expect_length(alerts(), 0)
  expect_length(notes(), 2)
  expect_match(notes()[1], "1 sample has no row")
  expect_identical(notes()[2], "Annotations have been added")
})

test_that("unsupported strand genomes are rejected before overwriting (#227)", {
  for (choice in c("SBS192 - Transcript_Strand",
                   "SBS192 - Replication_Strand")) {
    modals <- record_modals()
    alerts <- record_shinyalerts()
    musica_strand <- musica_sbs96
    # Stands in for an existing SBS192 table, so the overwrite prompt applies
    tables(musica_strand)[["SBS192_Trans"]] <- tables(musica_strand)$SBS96
    tables(musica_strand)[["SBS192_Rep"]] <- tables(musica_strand)$SBS96
    env <- shiny_server_env()
    shiny::testServer(env$server, {
      vals$musica <- musica_strand
      session$setInputs(select_table = choice, table_genome_list = "mm10")
      session$setInputs(add_table = 1)
      expect_false(session$isClosed())
      expect_identical(vals$musica, musica_strand)
    })
    expect_identical(modals(), 0)
    expect_length(alerts(), 1)
    expect_match(alerts()[[1]]$text, "genomes: hg19, hg38\\.$")
  }
})

test_that("an unsupported strand genome is rejected for a new table (#223)", {
  alerts <- record_shinyalerts()
  env <- shiny_server_env()
  shiny::testServer(env$server, {
    vals$musica <- musica_sbs96
    session$setInputs(select_table = "SBS192 - Transcript_Strand",
                      table_genome_list = "mm9")
    session$setInputs(add_table = 1)
    expect_false(session$isClosed())
    expect_identical(vals$musica, musica_sbs96)
  })
  expect_length(alerts(), 1)
  expect_match(alerts()[[1]]$text, "genomes: hg19, hg38\\.$")
})

test_that("a supported strand genome still builds the table", {
  notes <- record_notifications()
  alerts <- record_shinyalerts()
  env <- shiny_server_env()
  shiny::testServer(env$server, {
    vals$musica <- musica_sbs96
    session$setInputs(select_table = "SBS192 - Replication_Strand",
                      table_genome_list = "hg38")
    session$setInputs(add_table = 1)
    expect_false(session$isClosed())
    expect_true("SBS192_Rep" %in% names(tables(vals$musica)))
  })
  expect_length(alerts(), 0)
  expect_true("Table created." %in% notes())
})

test_that("warnings from building a table are shown as notifications (#223)", {
  local_mocked_bindings(
    build_standard_table = function(...) {
      warning("2 of 9 variants have motifs not in the SBS96 schema")
    },
    .package = "musicatk"
  )
  notes <- record_notifications()
  alerts <- record_shinyalerts()
  env <- shiny_server_env()
  shiny::testServer(env$server, {
    vals$musica <- musica_sbs96
    session$setInputs(select_table = "DBS", table_genome_list = "hg38")
    session$setInputs(add_table = 1)
    expect_false(session$isClosed())
  })
  expect_length(alerts(), 0)
  expect_identical(
    notes(),
    c("2 of 9 variants have motifs not in the SBS96 schema", "Table created.")
  )
})

test_that("Create Table without a genome shows an alert", {
  alerts <- record_shinyalerts()
  env <- shiny_server_env()
  shiny::testServer(env$server, {
    vals$musica <- musica_sbs96
    session$setInputs(select_table = "SBS96")
    session$setInputs(add_table = 1)
    expect_false(session$isClosed())
  })
  expect_length(alerts(), 1)
  expect_match(alerts()[[1]]$text, "genomes: hg19, hg38, mm9, mm10\\.$")
})

context("Test Plotting Functions")
library("musicatk")

test_that(desc = "Testing Exposures Plotting", {
  data(res)
  p <- plot_exposures(res, model_name = "res", plot_type = "bar")
  expect_s3_class(p, "ggplot")
  p <- plot_exposures(res, model_name = "res", plot_type = "violin")
  expect_s3_class(p, "ggplot")
  p <- plot_exposures(res, model_name = "res", plot_type = "box")
  expect_s3_class(p, "ggplot")
})

test_that(desc = "Box and violin plots work with group_by = 'none'", {
  data(res)
  p <- plot_exposures(res, model_name = "res", plot_type = "box")
  expect_s3_class(p$layers[[1]]$geom, "GeomBoxplot")
  expect_identical(rlang::as_label(p$mapping$x), "signature")
  p <- plot_exposures(res, model_name = "res", plot_type = "violin")
  expect_s3_class(p$layers[[1]]$geom, "GeomViolin")
  expect_identical(rlang::as_label(p$mapping$x), "signature")
})

test_that(desc = "threshold sets exposures below it to 0", {
  data(res)
  exposure_matrix <- exposures(res, "result", "SBS96", "res")
  full <- plot_exposures(res, model_name = "res")$data
  p <- plot_exposures(res, model_name = "res", threshold = 100)
  expect_equal(nrow(p$data), length(exposure_matrix))
  observed <- p$data[order(p$data$sample, p$data$signature), ]
  full <- full[order(full$sample, full$signature), ]
  expected <- ifelse(full$exposure < 100, 0, full$exposure)
  expect_equal(observed$exposure, expected)
})

test_that(desc = "Proportional exposures name samples with zero exposure", {
  data(res)
  exposure_matrix <- exposures(res, "result", "SBS96", "res")
  exposure_matrix[, 2] <- 0
  exposures(res, "result", "SBS96", "res") <- exposure_matrix
  expect_error(plot_exposures(res, model_name = "res", proportional = TRUE),
               regexp = colnames(exposure_matrix)[2])
})

test_that(desc = ".pivot_signatures matches a color_variable vector by motif", {
  data(res)
  sigs <- signatures(res, "result", "SBS96", "res")
  tab <- tables(res)[["SBS96"]]
  annot <- get_annot_tab(tab)
  rownames(annot) <- annot$motif
  color_vector <- as.character(annot[rownames(sigs), "mutation"])

  out <- .pivot_signatures(sigs, tab, color_variable = color_vector)
  expected <- color_vector[match(out$df$motif, rownames(sigs))]
  expect_equal(as.character(out$df$mutation_color), expected)
})

test_that(desc = ".pivot_signatures uses a color_variable column name", {
  data(res)
  sigs <- signatures(res, "result", "SBS96", "res")
  tab <- tables(res)[["SBS96"]]
  annot <- get_annot_tab(tab)
  rownames(annot) <- annot$motif

  out <- .pivot_signatures(sigs, tab, color_variable = "context")
  expect_equal(as.character(out$df$mutation_color),
               as.character(annot[out$df$motif, "context"]))
})

# Capture the two panels plot_signatures() passes to ggpubr::ggarrange()
.capture_signature_panels <- function(...) {
  captured <- NULL
  local_mocked_bindings(
    ggarrange = function(...) {
      captured <<- list(...)
      NULL
    },
    .package = "ggpubr"
  )
  plot_signatures(...)
  captured
}

.layers_of <- function(p, geom_class) {
  # Exact class: GeomBar also inherits from GeomRect
  Filter(function(layer) class(layer$geom)[1] == geom_class, p$layers)
}

test_that(desc = "Signature header colours are matched by name", {
  data(res)
  color_mapping <- get_color_mapping(tables(res)[["SBS96"]])
  # Swap the C>A and C>G colours, so C>A is the dark one, and reverse
  color_mapping[c("C>A", "C>G")] <- color_mapping[c("C>G", "C>A")]
  color_mapping <- rev(color_mapping)
  panels <- .capture_signature_panels(res, model_id = "res",
                                      color_mapping = color_mapping)
  header <- panels[[1]]
  rect <- .layers_of(header, "GeomRect")[[1]]
  boxes <- rect$data
  expect_equal(unname(as.character(rect$aes_params$fill)),
               unname(color_mapping[boxes$mutation_color]))

  # Label text is white on the dark C>A box only
  text <- .layers_of(header, "GeomText")[[1]]
  expected <- ifelse(text$data$mutation_color == "C>A", "white", "black")
  expect_equal(unname(text$aes_params$colour), expected)
})

test_that(desc = "Signature header and bars agree with a custom variable", {
  data(res)
  panels <- .capture_signature_panels(res, model_id = "res",
                                      color_variable = "context")
  header <- panels[[1]]
  rect <- .layers_of(header, "GeomRect")[[1]]
  expect_false(anyNA(rect$aes_params$fill))
  expect_equal(nrow(rect$data), 16)
})

test_that(desc = "Unnamed or incomplete color_mapping", {
  data(res)
  color_mapping <- get_color_mapping(tables(res)[["SBS96"]])
  panels <- .capture_signature_panels(res, model_id = "res",
                                      color_mapping = unname(color_mapping))
  rect <- .layers_of(panels[[1]], "GeomRect")[[1]]
  expect_equal(unname(as.character(rect$aes_params$fill)),
               unname(color_mapping[rect$data$mutation_color]))

  expect_error(plot_signatures(res, model_id = "res",
                               color_mapping = color_mapping[-2]),
               regexp = "C>G")
})

test_that(desc = "plot_signatures ignores a global 'annotation_text'", {
  data(res)
  assign("annotation_text",
         data.frame(label = "global", signature = "x",
                    mutation_color = "C>A"),
         envir = globalenv())
  withr::defer(rm("annotation_text", envir = globalenv()))
  panels <- .capture_signature_panels(res, model_id = "res")
  expect_length(.layers_of(panels[[2]], "GeomText"), 1)
})

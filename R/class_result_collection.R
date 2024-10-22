# Result Collection object -------------------------------

#' The Result Collection object that contains modality, input parameters, prior hyperparameters
#' @slot modality a list contains model results for different modality
#' @slot parameter a list contains input parameters
#' @slot hyperparameter a list contains prior and tuning parameters
#' @importFrom S4Vectors SimpleList
#' @exportClass result_collection
#' @export
setClass("result_collection", slots = c(modality = "SimpleList", 
                                        parameter = "list", 
                                        hyperparameter= "list"),
         prototype = list(modality = SimpleList(), 
                          parameter = list(), 
                          hyperparameter = list())
)

#' @title Retrieve a specific modality entry from a musica or result_collection object
#' @description  \code{modality} list contains model results for a modality
#' @param x  A \code{\linkS4class{result_model}} or
#' \code{\linkS4class{result_collection}} object
#' @param result The name of the result_list entry.
#' @param modality The modality.
#' @param ... Other inputs
#' @rdname get_modality
#' @return A list of modality which contains result_model objects
#' @export
#' @examples
#' data(res)
#' get_modality(res, "result", "SBS96")
setGeneric(
  name = "get_modality",
  def = function(x, ...)
  {
    standardGeneric("get_modality")
  }
)

#' @rdname get_modality
setMethod(
  f = "get_modality",
  signature = "musica",
  definition = function(x, result, modality) {
    return(x@result_list[[result]]@modality[[modality]])
  }
)

#' @rdname get_modality
setMethod(
  f = "get_modality",
  signature = "result_collection",
  definition = function(x, modality) {
    return(x@modality[[modality]])
  }
)

#' @title Retrieve parameter from a musica or result_collection object
#' @description  The \code{parameter} contains input parameters used in the model 
#' @param x  A \code{\linkS4class{result_model}} or
#' \code{\linkS4class{result_collection}} object
#' @param result The name of the result_list entry.
#' @param ... Other inputs
#' @rdname parameter
#' @return a list of parameters
#' @export
#' @examples
#' data(res)
#' parameter(res, "result")
setGeneric(
  name = "parameter",
  def = function(x, ...)
  {
    standardGeneric("parameter")
  }
)

#' @rdname parameter 
setMethod(
  f = "parameter",
  signature = "musica",
  definition = function(x, result) {
    return(x@result_list[[result]]@parameter)
  }
)

#' @rdname parameter 
setMethod(
  f = "parameter",
  signature = "result_collection",
  definition = function(x) {
    return(x@parameter)
  }
)

#' @param x  A \code{\linkS4class{result_model}} or
#' \code{\linkS4class{result_collection}} object
#' @param value a list of input parameters
#' @rdname parameter
#' @export

setGeneric(
  name = "parameter<-",
  def = function(x, ..., value)
  {
    standardGeneric("parameter<-")
  }
)

#' @rdname parameter
setReplaceMethod(
  f = "parameter",
  signature = c("result_collection", "list"),
  definition = function(x, value)
  {
    x@parameter <- value
    return(x)
  }
) 

#' @rdname parameter
setReplaceMethod(
  f = "parameter",
  signature = c("musica", "list"),
  definition = function(x, result, value)
  {
    x@result_list[[result]]@parameter <- value
    return(x)
  }
) 

#' @title Retrieve hyperparameter from a musica or result_collection object
#' @description The \code{hyperparameter} contain list of prior and tuning parameters 
#' @param x  A \code{\linkS4class{result_model}} or
#' \code{\linkS4class{result_collection}} object
#' @param result The name of the result_list entry.
#' @param ... Other inputs
#' @rdname hyperparameter
#' @return A list of hyperparameters
#' @export
#' @examples
#' data(res)
#' hyperparameter(res, "result")
setGeneric(
  name = "hyperparameter",
  def = function(x, ...)
  {
    standardGeneric("hyperparameter")
  }
)

#' @rdname hyperparameter
setMethod(
  f = "hyperparameter",
  signature = "musica",
  definition = function(x, result) {
    return(x@result_list[[result]]@hyperparameter)
  }
)

#' @rdname hyperparameter
setMethod(
  f = "hyperparameter",
  signature = "result_collection",
  definition = function(x) {
    return(x@hyperparameter)
  }
)

#' @rdname hyperparameter
#' @param x  A \code{\linkS4class{result_model}} or
#' \code{\linkS4class{result_collection}} object
#' @param value A \code{\linkS4class{list}} of hyperparameters for model
#' @export
setGeneric(
  name = "hyperparameter<-",
  def = function(x, ..., value)
  {
    standardGeneric("hyperparameter<-")
  }
)

#' @rdname hyperparameter
setReplaceMethod(
  f = "hyperparameter",
  signature = c("musica", "list"),
  definition = function(x, result, value)
  {
    x@result_list[[result]]@hperparameter
    return(x)
  }
)

#' @rdname hyperparameter
setReplaceMethod(
  f = "hyperparameter",
  signature = c("result_collection", "list"),
  definition = function(x, value)
  {
    x@hyperparameter <- value
    return(x)
  }
) 


#' @title Retrieve model from a musica or result collection object
#' @description Extract the \code{\linkS4class{result_model}} object from the
#' \code{\linkS4class{musica}} or \code{\linkS4class{result_collection}} object
#' that contains the model.
#' @param x A \code{\linkS4class{musica}} or
#' \code{\linkS4class{result_collection}} object.
#' @param result The name of the result_list entry.
#' @param modality The modality.
#' @param model The name of the model.
#' @param ... Other inputs
#' @rdname get_model
#' @return A \code{\linkS4class{result_model}} object
#' @export
#' @examples
#' data(res)
#' get_model(res, "result", "SBS96", "res")
setGeneric(
  name = "get_model",
  def = function(x, ...)
  {
    standardGeneric("get_model")
  }
)

#' @rdname get_model
setMethod(
  f = "get_model",
  signature = "musica",
  definition = function(x, result, modality, model) {
    return(x@result_list[[result]]@modality[[modality]][[model]])
  }
)

#' @rdname get_model
setMethod(
  f = "get_model",
  signature = "result_collection",
  definition = function(x, modality, model) {
    return(x@modality[[modality]][[model]])
  }
)






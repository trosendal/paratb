##' ncalf
##'
##' A method to calculate the number of calves.
##'
##' @title ncalf
##' @return A character string
ncalf <- function() {
    "S1 + T1H + T1L"
}

##' icalf
##'
##' A method to calculate the number of infected calves.
##'
##' @title icalf
##' @return A character string
icalf <- function() {
    "T1H + T1L"
}

##' nheifer
##'
##' A method to calculate the number of heifers.
##'
##' @title nheifer
##' @return A character string
nheifer <- function() {
   "S2 + T2L + T2H + E2L + E2H"
}

##' iheifer
##'
##' A method to calculate the number of infected heifers.
##'
##' @title iheifer
##' @return A character string
iheifer <- function() {
       "T2L + T2H + E2L + E2H"
}

##' ncow
##'
##' A method to calculate the the number of cows.
##'
##' @title ncow
##' @return A character string
ncow <- function() {
    "S3 + E3H + E3L + L3H + L3L + H3"
}

##' icow
##'
##' A method to calculate the number of infected cows.
##'
##' @title icow
##' @return A character string
icow <- function() {
    "E3H + E3L + L3H + L3L + H3"
}

##' imilk
##'
##' A method to calculate the number of positive milk tests.
##'
##' @title imilk
##' @return A character string
imilk <- function() {
    "BMp"
}

##' nmilk
##'
##' A method to calculate the number of milk tests.
##'
##' @title nmilk
##' @return A character string
nmilk <- function() {
    "BMp + BMn"
}
##' prev_class
##'
##' @title prevclass
##' @param class Age class of the prevalence
##' @return A formula
##' @importFrom stats formula
##' @export
##' @author Thomas Rosendal
prev_class <- function(class = c("true", "calf", "heifer", "cow", "apparent", "milk"))
{
    class <- match.arg(class)
    class <- switch(class,
                    true = paste(icalf(),   "+",
                                 iheifer(), "+",
                                 icow(),    "~",
                                 ncalf(),   "+",
                                 nheifer(), "+",
                                 ncow())
                                   ,
                    calf = paste(icalf(),   "~",
                                 ncalf())
                                   ,
                    heifer = paste(iheifer(),   "~",
                                   nheifer())
                                     ,
                    cow = paste(icow(),   "~",
                                ncow())
                                  ,
                    apparent = paste(imilk(), "~",
                                     nmilk())
                                       ,
                    milk = paste(imilk(),     "~",
                                 nmilk())
                    )
    return(formula(class))
}

##' A method for calculating prevalance in a paratb model
##'
##' This just assembles the pieces to pass onto the SimInf::prevalance
##' function.
##'
##' @title prev
##' @param model A paratb model object
##' @param class The age category or sample type that you want the
##'     prevalance from
##' @param type the type pf prevalance
##' @param ... Other arguments
##' @return A data.frame from the SimInf::prevalance method
##' @export
##' @author Thomas Rosendal
prev <- function(model,
                 class = c("true", "calf", "heifer", "cow", "apparent", "milk"),
                 type = c("nop", "wnp"),
                 ...)
{
    type <- match.arg(type)
    class <- match.arg(class)
    class <- prev_class(class)
    SimInf::prevalence(model, class, type = type, ...)
}

##' A method to count the number in some compartments.
##'
##'
##' Similar to prev but without the denomenator
##'
##' @title counti
##' @param model A paratb model object
##' @param class The age category or sample type that you want the
##'     prevalance from
##' @param type The type of prevalence
##' @return A data.frame from the SimInf::prevalance method
##' @return A data.frame from the SimInf::prevalance method
##' @importFrom stats formula
##' @export
##' @author Thomas Rosendal
counti <- function(model,
                   class = c("true", "calf", "heifer", "cow", "apparent", "milk", "nmilk"),
                   type = c("nop", "wnp"))
{
    type <- match.arg(type)
    class <- match.arg(class)
    compartments <- switch(class,
                           true = paste("~", icalf(),   "+",
                                             iheifer(), "+",
                                             icow()
                                                ),
                           calf = paste("~", icalf()),
                           heifer = paste("~", iheifer()),
                           cow = paste("~", icow()),
                           apparent = paste("~", imilk()),
                           milk = paste("~", imilk()),
                           nmilk = paste("~", nmilk())
                           )
    res <- trajectory(model, compartments = formula(compartments))
    if(type == "wnp") return(res)
    resb <- res[,!(names(res) %in% c("node", "time"))]
    if(is.null(dim(resb))) resb <- matrix(resb, ncol = 1)
    data.frame(time = model@tspan,
               count = tapply(rowSums(resb) > 0, res$time, "sum"))
}

##' A method to count the number of herds missed by the surveillance
##'
##' @title missed
##' @param model A paratb model object
##' @return A data.frame
##' @export
##' @author Thomas Rosendal
missed <- function(model)
{
    true_state <- counti(model, "true", "wnp")
    true_state$pos <- rowSums(true_state[,!(names(true_state) %in% c("node", "time"))])
    test_state <- counti(model, "apparent", "wnp")
    test_state$pos <- rowSums(test_state[,!(names(test_state) %in% c("node", "time"))])
    true_state$missed <- true_state$pos > 0 & test_state$pos == 0
    data.frame(time = model@tspan,
               count = tapply(true_state$missed > 0, true_state$time, "sum"), stringsAsFactors = FALSE)
}

##' A method to list the positive herds for each tstep
##'
##' @title which_pos
##' @param model A paratb model object
##' @return A data.frame
##' @export
##' @author Thomas Rosendal
which_pos <- function(model)
{
    true_state <- counti(model, "true", "wnp")
    true_state$pos <- rowSums(true_state[,!(names(true_state) %in% c("node", "time"))])
    data.frame(time = model@tspan,
               pos = do.call("c", lapply(unique(true_state$time), function(x){
                   paste(c("herds:", true_state$node[true_state$time == x & true_state$pos > 0]), collapse = ",")
               })), stringsAsFactors = FALSE
               )
}

##' A method to list the detected herds for each tstep
##'
##' @title which_test_pos
##' @param model A paratb model object
##' @return A data.frame
##' @export
##' @author Thomas Rosendal
which_test_pos <- function(model)
{
    test_state <- counti(model, "apparent", "wnp")
    test_state$pos <- rowSums(test_state[,!(names(test_state) %in% c("node", "time"))])
    data.frame(time = model@tspan,
               pos = do.call("c", lapply(unique(test_state$time), function(x){
                   paste(c("herds:", test_state$node[test_state$time == x & test_state$pos > 0]), collapse = ",")
               })), stringsAsFactors = FALSE
               )
}

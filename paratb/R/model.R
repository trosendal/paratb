##' Class \code{paratb}
##'
##' Class to handle the \code{paratb} \code{SimInf_model}.
##' @export
setClass("paratb", contains = "SimInf_model")

##' Create a model for the SimInf framework
##'
##' Create a model to be used by the SimInf framework.
##' @param u0 A data.frame with the initial state in each node.
##' @param tspan A vector (length >= 2) of increasing time points
##'     where the state of each node is to be returned.
##' @import SimInf
##' @import methods
##' @import Matrix
##' @export
##' @examples
##' ## Please add example(s) how to use the model
paratb <- function(u0 = NULL, tspan = NULL) {

    compartments <- c("S0", "S1", "T1H", "T1L", "S2", "T2H", "T2L",
                      "E2H", "E2L", "S3", "E3H", "E3L", "L3H", "L3L",
                      "H3", "BM", "BMp", "BMn")

    ## Check tspan (if null, make it 1 year)
    if(is.null(tspan)) {
        tspan <- 1:365
    }

    ## Check u0
    if(!is.null(u0)) {
        if (!is.data.frame(u0))
            u0 <- as.data.frame(u0)
        if (!all(compartments %in% names(u0)))
            stop("Missing columns in u0")
        u0 <- u0[, compartments, drop = FALSE]
    } else {

    ## Start with 10 farms with 200 animals each without any infection
    u0 <- structure(list(S0 = c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
                         S1 = c(50, 50, 50, 50, 50, 50, 50, 50, 50, 50),
                         T1H = c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
                         T1L = c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
                         S2 = c(50, 50, 50, 50, 50, 50, 50, 50, 50, 50),
                         T2H = c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
                         T2L = c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
                         E2H = c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
                         E2L = c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
                         S3 = c(100, 100, 100, 100, 100, 100, 100, 100, 100, 100),
                         E3H = c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
                         E3L = c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
                         L3H = c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
                         L3L = c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
                         H3 = c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
                         BM = c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
                         BMp = c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
                         BMn = c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L)),
                    row.names = c(NA, 10L), class = "data.frame")
    }

    G <- new("dgCMatrix", i = c(1L, 2L, 3L, 1L, 2L, 3L, 4L, 5L, 6L,
                                9L, 12L, 13L, 26L, 1L, 2L, 3L, 4L, 5L,
                                7L, 10L, 12L, 13L, 26L, 1L, 2L, 3L,
                                4L, 5L, 8L, 11L, 12L, 13L, 26L, 4L,
                                5L, 6L, 7L, 9L, 10L, 12L, 13L, 26L,
                                4L, 5L, 6L, 8L, 9L, 11L, 12L, 13L,
                                26L, 4L, 5L, 6L, 9L, 12L, 13L, 16L,
                                21L, 26L, 4L, 5L, 7L, 10L, 12L, 13L,
                                14L, 19L, 22L, 26L, 4L, 5L, 8L, 11L,
                                12L, 13L, 15L, 20L, 23L, 26L, 4L, 5L,
                                6L, 9L, 12L, 13L, 26L, 4L, 5L, 7L,
                                10L, 12L, 13L, 26L, 4L, 5L, 8L, 11L,
                                12L, 13L, 26L, 4L, 5L, 12L, 13L, 14L,
                                16L, 19L, 21L, 22L, 26L, 4L, 5L, 12L,
                                13L, 15L, 16L, 20L, 21L, 23L, 26L, 4L,
                                5L, 12L, 13L, 14L, 17L, 19L, 22L, 24L,
                                26L, 4L, 5L, 12L, 13L, 15L, 18L, 20L,
                                23L, 25L, 26L, 0L, 1L, 2L, 3L, 4L, 5L,
                                12L, 13L, 16L, 21L, 26L, 30L, 39L,
                                40L, 41L, 0L, 1L, 2L, 3L, 4L, 5L, 12L,
                                13L, 17L, 24L, 26L, 27L, 31L, 39L,
                                40L, 41L, 0L, 1L, 2L, 3L, 4L, 5L, 12L,
                                13L, 18L, 25L, 26L, 28L, 32L, 39L,
                                40L, 41L, 0L, 1L, 2L, 3L, 4L, 5L, 12L,
                                13L, 14L, 19L, 22L, 26L, 27L, 31L,
                                39L, 40L, 41L, 0L, 1L, 2L, 3L, 4L, 5L,
                                12L, 13L, 15L, 20L, 23L, 26L, 28L,
                                32L, 39L, 40L, 41L, 4L, 5L, 12L, 13L,
                                16L, 21L, 26L, 4L, 5L, 12L, 13L, 14L,
                                19L, 22L, 26L, 4L, 5L, 12L, 13L, 15L,
                                20L, 23L, 26L, 4L, 5L, 12L, 13L, 17L,
                                24L, 26L, 4L, 5L, 12L, 13L, 18L, 25L,
                                26L, 0L, 1L, 2L, 3L, 4L, 5L, 12L, 13L,
                                26L, 28L, 30L, 32L, 39L, 40L, 41L, 0L,
                                1L, 2L, 3L, 4L, 5L, 12L, 13L, 26L,
                                27L, 29L, 31L, 33L, 36L, 39L, 40L,
                                41L, 0L, 1L, 2L, 3L, 4L, 5L, 12L, 13L,
                                26L, 28L, 32L, 34L, 37L, 39L, 40L,
                                41L, 0L, 1L, 2L, 3L, 4L, 5L, 12L, 13L,
                                26L, 29L, 33L, 35L, 36L, 38L, 39L,
                                40L, 41L, 0L, 1L, 2L, 3L, 4L, 5L, 12L,
                                13L, 26L, 30L, 39L, 40L, 41L, 0L, 1L,
                                2L, 3L, 4L, 5L, 12L, 13L, 26L, 27L,
                                31L, 39L, 40L, 41L, 0L, 1L, 2L, 3L,
                                4L, 5L, 12L, 13L, 26L, 28L, 32L, 39L,
                                40L, 41L, 0L, 1L, 2L, 3L, 4L, 5L, 12L,
                                13L, 26L, 29L, 33L, 36L, 39L, 40L,
                                41L, 0L, 1L, 2L, 3L, 4L, 5L, 12L, 13L,
                                26L, 34L, 37L, 39L, 40L, 41L, 0L, 1L,
                                2L, 3L, 4L, 5L, 12L, 13L, 26L, 35L,
                                38L, 39L, 40L, 41L, 0L, 1L, 2L, 3L,
                                4L, 5L, 12L, 13L, 26L, 29L, 30L, 33L,
                                36L, 39L, 40L, 41L, 0L, 1L, 2L, 3L,
                                4L, 5L, 12L, 13L, 26L, 30L, 34L, 37L,
                                39L, 40L, 41L, 0L, 1L, 2L, 3L, 4L, 5L,
                                12L, 13L, 26L, 30L, 35L, 38L, 39L,
                                40L, 41L, 40L, 41L, 40L, 41L, 40L,
                                41L),
             p = c(0L, 3L, 13L, 23L, 33L, 42L, 51L, 60L, 70L, 80L,
                   87L, 94L, 101L, 111L, 121L, 131L, 141L, 156L, 172L,
                   188L, 205L, 222L, 229L, 237L, 245L, 252L, 259L,
                   274L, 291L, 307L, 324L, 337L, 351L, 365L, 380L,
                   394L, 408L, 424L, 439L, 454L, 456L, 458L, 460L), Dim = c(42L, 42L),
             Dimnames = list(c("@ -> S0", "S0 -> S1", "S0 -> T1H",
                               "S0 -> T1L", "S1 -> T1H", "S1 -> T1L",
                               "S1 -> S2", "T1H -> T2H", "T1L -> T2L",
                               "S1 -> @", "T1H -> @", "T1L -> @",
                               "S2 -> T2H", "S2 -> T2L", "T2H -> E2H",
                               "T2L -> E2L", "S2 -> S3", "E2H -> E3H",
                               "E2L -> E3L", "T2H -> E3H",
                               "T2L -> E3L", "S2 -> @", "T2H -> @",
                               "T2L -> @", "E2H -> @", "E2L -> @",
                               "S3 -> E3L", "E3H -> L3H",
                               "E3L -> L3L", "L3H -> H3", "S3 -> @",
                               "E3H -> @", "E3L -> @", "L3H -> @",
                               "L3L -> @", "H3 -> @", "L3H -> S3",
                               "L3L -> S3", "H3 -> S3", "@ -> BM",
                               "BM -> BMp", "BM -> BMn"),
                             c("1", "2", "3", "4", "5", "6", "7", "8",
                               "9", "10", "11", "12", "13", "14",
                               "15", "16", "17", "18", "19", "20",
                               "21", "22", "23", "24", "25", "26",
                               "27", "28", "29", "30", "31", "32",
                               "33", "34", "35", "36", "37", "38",
                               "39", "43", "44", "45")),
             x = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                   1), factors = list())

    S <- new("dgCMatrix", i = c(0L, 0L, 1L, 0L, 2L, 0L, 3L, 1L, 2L,
                                1L, 3L, 1L, 4L, 2L, 5L, 3L, 6L, 1L,
                                2L, 3L, 4L, 5L, 4L, 6L, 5L, 7L, 6L,
                                8L, 4L, 9L, 7L, 10L, 8L, 11L, 5L, 10L,
                                6L, 11L, 4L, 5L, 6L, 7L, 8L, 9L, 11L,
                                10L, 12L, 11L, 13L, 12L, 14L, 9L, 10L,
                                11L, 12L, 13L, 14L, 9L, 12L, 9L, 13L,
                                9L, 14L, 15L, 15L, 16L, 15L, 17L),
             p = c(0L, 1L, 3L, 5L, 7L, 9L, 11L, 13L, 15L, 17L, 18L,
                   19L, 20L, 22L, 24L, 26L, 28L, 30L, 32L, 34L, 36L,
                   38L, 39L, 40L, 41L, 42L, 43L, 45L, 47L, 49L, 51L,
                   52L, 53L, 54L, 55L, 56L, 57L, 59L, 61L, 63L, 64L,
                   66L, 68L), Dim = c(18L, 42L),
             Dimnames = list(
                 c("S0", "S1", "T1H", "T1L", "S2", "T2H", "T2L",
                   "E2H", "E2L", "S3", "E3H", "E3L", "L3H", "L3L",
                   "H3", "BM", "BMp", "BMn"),
                 c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                   "11", "12", "13", "14", "15", "16", "17", "18",
                   "19", "20", "21", "22", "23", "24", "25", "26",
                   "27", "28", "29", "30", "31", "32", "33", "34",
                   "35", "36", "37", "38", "39", "43", "44", "45")),
             x = c(1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1,
                   -1, 1, -1, -1, -1, -1, 1, -1, 1, -1, 1, -1, 1, -1,
                   1, -1, 1, -1, 1, -1, 1, -1, 1, -1, -1, -1, -1, -1,
                   -1, 1, -1, 1, -1, 1, -1, 1, -1, -1, -1, -1, -1, -1,
                   1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1), factors = list())

    E <- methods::as(matrix(integer(0), nrow = 0, ncol = 0), "dgCMatrix")

    N <- matrix(integer(0), nrow = 0, ncol = 0)

    v0 <- matrix(numeric(0), nrow  = 0, ncol = nrow(u0))
    storage.mode(v0) <- "double"

    ldata <- matrix(numeric(0), nrow = 0, ncol = nrow(u0))
    storage.mode(ldata) <- "double"

    ## By default don't add any events
    events <- NULL

    gdata <- structure(c(1000, 0.01, 0.04, 0.04, 0.00273973,
                         0.00273973, 0.00547945, 0.335, 0.0589041,
                         0.00145205, 0.00468493, 0.00123288,
                         0.00273973, 0.000246575, 0.00273973,
                         2.73973e-05, 0.00112329, 0.000136986,
                         7.980511, -2.048039, 0, 1),
                       .Names = c("p1", "gamma_E", "gamma_L",
                                  "gamma_H", "beta", "beta_a", "phi",
                                  "eta", "sigma_H", "sigma_L", "nu",
                                  "mu_b", "rho_1", "mu_1", "rho_2",
                                  "mu_2", "mu_3", "mu_4", "intercept",
                                  "coef", "delta2", "Sp2"))

    model <- SimInf_model(G      = G,
                          S      = S,
                          E      = E,
                          N      = N,
                          tspan  = tspan,
                          events = events,
                          ldata  = ldata,
                          gdata  = gdata,
                          u0     = u0,
                          v0     = v0)

    methods::as(model, "paratb")
}

##' Run the model
##'
##' @rdname run-methods
##' @param model The model to run.
##' @param threads Number of threads. Default is NULL, i.e. to use all
##'     available processors.
##' @param solver Which numerical solver to utilize. Default is 'ssa'.
##' @return model with result from simulation.
##' @export
##' @useDynLib paratb, .registration=TRUE
setMethod("run",
    signature(model = "paratb"),
    function(model, threads = NULL, solver = NULL)
    {
        methods::validObject(model)
        .Call(SimInf_model_run, model, threads, solver,
              PACKAGE = "paratb")
    })

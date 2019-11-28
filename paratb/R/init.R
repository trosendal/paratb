##' random_events
##'
##' A method to add random external transfer events to a model. This
##' is to illustrate the model working with events despite that we
##' cannot share the recording events used in the study.
##'
##' @title random_events
##' @param model The model that the random event should be added to
##' @param n The number of random events
##' @import Matrix
##' @export
##' @return A model object
random_events <- function(model, n) {

    ## All events are external transfers
    event <- 3L

    ## Sample the vector from minimum to maximum of tspan
    time <- sample(seq(min(model@tspan), max(model@tspan)), n, replace = TRUE)

    ## sample pairs of nodes for node and dest
    nodes <- sapply(1:n, function(i){
        sample(seq_len(ncol(model@u0)), 2, replace = FALSE)
    })
    node <- nodes[1,]
    dest <- nodes[2,]

    ## n is just any number because we will use the proportion
    number <- 1L

    ## Sample the proportion from between 1% and 10% of the animals in
    ## the compartment(s) for the "select".
    proportion <- round(runif(n, 0.01, 0.05), 2)

    ## The select can be one of the 3 external transfers 2, 4, 5
    select <- sample(c(2, 4, 5), n, replace = TRUE)

    ## Apply no shift
    shift <- 0

    ## Stick it in a data.frame
    events <- data.frame(event = event,
                         time = as.integer(time),
                         node = as.numeric(node),
                         dest = as.numeric(dest),
                         n = number,
                         proportion = as.numeric(proportion),
                         select = as.numeric(select),
                         shift = as.numeric(shift))

    ## select matrix
    E <- Matrix(c(1, 1, 0, 0, 0, 0, 0, 0, 0, #S0
                  0, 1, 1, 0, 0, 0, 0, 0, 0, #S1
                  0, 1, 1, 0, 0, 1, 0, 0, 0, #T1H
                  0, 1, 1, 0, 0, 1, 0, 0, 0, #T1L
                  0, 0, 0, 1, 0, 0, 0, 0, 0, #S2
                  0, 0, 0, 1, 0, 0, 1, 0, 0, #T2H
                  0, 0, 0, 1, 0, 0, 1, 0, 0, #T2L
                  0, 0, 0, 1, 0, 0, 1, 0, 0, #E2H
                  0, 0, 0, 1, 0, 0, 1, 0, 0, #E2L
                  0, 0, 0, 0, 1, 0, 0, 0, 0, #S3
                  0, 0, 0, 0, 1, 0, 0, 1, 0, #E3H
                  0, 0, 0, 0, 1, 0, 0, 1, 0, #E3L
                  0, 0, 0, 0, 1, 0, 0, 1, 0, #L3H
                  0, 0, 0, 0, 1, 0, 0, 1, 0, #L3L
                  0, 0, 0, 0, 1, 0, 0, 1, 0, #H3
                  0, 0, 0, 0, 0, 0, 0, 0, 1, #BM
                  0, 0, 0, 0, 0, 0, 0, 0, 0, #BMp
                  0, 0, 0, 0, 0, 0, 0, 0, 0),#BMn
                nrow   = 18,
                ncol   = 9,
                byrow  = TRUE,
                sparse = TRUE)
    row.names(E) <- c("S0", "S1", "T1H", "T1L", "S2", "T2H", "T2L", "E2H",
                      "E2L", "S3", "E3H", "E3L", "L3H", "L3L", "H3",
                      "BM", "BMp", "BMn")
    colnames(E) <- c("Calves born", "Calves Die/Move",
                     "Calves Age", "Heifers Age/Die/Move",
                     "Cows Age/Die/Move", "Seed a calf", "Seed a heifer",
                     "Seed a cow", "Milk sampling")

    ## shift matrix
    N <- matrix(c(0, 0, #S0
                  3, 0, #S1
                  3, 0, #T1H
                  3, 0, #T1L
                  0, 5, #S2
                  0, 5, #T2H
                  0, 5, #T2L
                  0, 3, #E2H
                  0, 3, #E2L
                  0, 0, #S3
                  0, 0, #E3H
                  0, 0, #E3L
                  0, 0, #L3H
                  0, 0, #L3L
                  0, 0, #H3
                  0, 0, #BM
                  0, 0, #BMp
                  0, 0), #BMn
                nrow   = 18,
                ncol   = 2,
                byrow  = TRUE)
    row.names(N) <- c("S0", "S1", "T1H", "T1L", "S2", "T2H", "T2L", "E2H",
                      "E2L", "S3", "E3H", "E3L", "L3H", "L3L", "H3",
                      "BM", "BMp", "BMn")
    colnames(N) <- c("Calves Age", "Heifers Age")

    ## return the model with added events
    model@events <- SimInf_events(E = E, N = N, events = events)
    model
}

##' seed_herd
##'
##' A method for seeding infected animals in a herd.
##'
##' @title seed_herd
##' @param model The model object
##' @param i the node to seed in
##' @param age the age category to seed in. One of c("calf", "heifer", "cow").
##' @param n the number of animals to seed
##' @param seed The seed
##' @return a paratb model object
##' @export
##' @author Thomas Rosendal
seed_herd <- function(model,
                      i,
                      age = c("calf", "heifer", "cow"),
                      n,
                      seed = NULL)
{
    if(!is.null(seed)) {
        if(exists(".Random.seed")) {
            savedSeed = .Random.seed
            on.exit(.Random.seed <<- savedSeed)
        }
        set.seed(seed)
    }
    age <- switch(age,
                  calf = c("S1", "T1H", "T1L"),
                  heifer = c("S2", "T2H", "T2L", "E2L", "E2H"),
                  cow = c("S3", "E3H", "E3L", "L3H", "L3L", "H3")
                  )
    n <- as.integer(n)
    ## Move animals to a random infected compartment until you are
    ## done or have no more susceptibles available in that age
    ## category:
    while(model@u0[age[1], i] > 0 & n > 0) {
        dest <- sample(age[-1], 1)
        model@u0[dest, i] <- model@u0[dest, i] + 1L
        model@u0[age[1], i] <- model@u0[age[1], i] - 1L
        n <- n - 1L
    }
    return(model)
}

##' add_sampling
##'
##' Add a couple of columns to the select matrix to do blood and milk
##' sampling. This will require scheduling birth events with select 9
##' and 10
##'
##' @title add_sampling
##' @param model A paratb model object
##' @return A modified model object
##' @author Thomas Rosendal
add_sampling <- function(model)
{
    model@events@E <- cbind(model@events@E,
                            Matrix(c(rep(0, 15), 1, rep(0, 23), 1, rep(0, 2)),
                                   ncol = 2,
                                   sparse = TRUE,
                                   dimnames = list(NULL, c("Blood sampling", "Milk sampling"))))
    model
}

##' Add surveillance events
##'
##' A method for adding surveillance events to the model. These are
##' events of Event type == 1 that increment the BL or BM compartments
##' at a specific time point. These require select = 9 for blood
##' sampling of adults and select = 10 for bulk milk sampling. When
##' these units are added to the model they are redistributed to the
##' appropriate test status by the simulator.
##'
##' @title add_surveillance_events
##' @param model The model object
##' @param i the node to add a surveillance event to
##' @param t the timepoint that you want to do surveillance
##' @param type One of either c("blood", "milk")
##' @param n the number of animals to seed
##' @param events A dataframe of the events the default is just to
##'     extract them from the model object.
##' @return a paratb model object with some extra events
##' @export
##' @author Thomas Rosendal
add_surveillance_event <- function(model,
                            i,
                            t = 1,
                            type = c("blood", "milk"),
                            n = 1,
                            events = as(model@events, "data.frame"))
{
    ## check that you added the appropriate select columns to the E
    ## matrix:
    model <- add_sampling(model)
    ## And turn off the stochastic samplings events
    model@gdata["delta2"] <- 0
    type <- switch(type,
                   blood = 9,
                   milk = 10)
    seed_events <- data.frame(event = 1L,
                              time = as.integer(t),
                              node = as.integer(i),
                              dest = 1L,
                              n = as.integer(n),
                              proportion = 0,
                              select = as.integer(type),
                              shift = 0L)
    events <- rbind(events, seed_events)
    model@events <- SimInf_events(model@events@E, model@events@N, events)
    model
}

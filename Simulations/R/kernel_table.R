

rej_f_tbl <- read.csv("Simulations/Products/crank_bs_ear1_sims_res.csv")
ar1_tbl <- read.csv("Simulations/Products/crank_bs_ear1_sims_ar1.csv")

id <- rbind(
    data.frame(Var1 = "No Splines", Var2 = ""),
    expand.grid(
        c(
        "Order 1",
        "Order 2"
        ),
        c(
        "4 splines",
        "12 splines",
        "24 splines",
        "36 splines"
        )
    )
)
ar1_tbl<-cbind(id, ar1_tbl)

rep(id[1:3,], 3)

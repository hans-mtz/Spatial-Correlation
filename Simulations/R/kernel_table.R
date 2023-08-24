

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
### ------


kernel_results <- read.csv("Simulations/Products/grid_sims_res.csv", na.strings="NaN")

e_ar1_tbl <- read.csv("Simulations/Products/grid_sims_ar1.csv")

row_names_kernel <- c(
    # " ",
    # "SCPC",
    rep(
        c(
            "HR",
            "Kernel-0.015",
            "Kernel-0.030",
            "Kernel-0.050"
            # "0.07"
        ),
        19
    )
)

row_names_ar1_kernel <- c(
    rep(
        "No splines",
        4
    ),
    rep(
        c(
        "4 splines",
        "12 splines",
        "24 splines",
        "36 splines",
        "48 splines",
        "64 splines",
        "88 splines",
        "102 splines",
        "124 splines"
        ),
        each=4,
        times =2
    )
)
aux <- c(
    rep(
        "",
        4
    ),
    rep(
        c(
            "Order 1",
            "Order 2"
        ),
        each = 36
    )
)
kernel_id<-paste0(aux, row_names_ar1_kernel)

table_kernel <- cbind(
    aux,
    row_names_ar1_kernel,
    row_names_kernel, 
    kernel_results, 
    kernel_id
)
names(table_kernel)<- c(
    "Spline order", "Spline qty.",
    "Estimator",
    "Cons.",
    "$\\beta_1$",
    "id"
)

row_names_ar1 <- c(
    "No splines",
    rep(
        c(
        "4 splines",
        "12 splines",
        "24 splines",
        "36 splines",
        "48 splines",
        "64 splines",
        "88 splines",
        "102 splines",
        "124 splines"
        ),
        2
    )
)
aux2 <- c(
    rep(
        "",
        1
    ),
    rep(
        c(
            "Order 1",
            "Order 2"
        ),
        each = 9
    )
)
ar1_id<-paste0(aux2, row_names_ar1)
table_ar1 <-cbind(
    ar1_id,
    e_ar1_tbl
)
names(table_ar1)<-c(
    "id",
    "Cons2",
    "$\\gamma$",
    "BIC"
)

table<-merge(table_kernel,table_ar1, by="id",sort=FALSE)

col_names <- c(
    "Spline order",
    "No. of splines",
    "Estimator",
    "Cons.",
    "$\\beta_1$",
    # "Cons.",
    "$\\gamma$",
    "BIC"
)

table[, -c(1,7)] |> head()

##


kernel_results <- read.csv("Simulations/Products/crank_bs_ear1_sims_res_dd_2x.csv")

e_ar1_tbl <- read.csv("Simulations/Products/crank_bs_ear1_sims_ar1_dd_2x.csv")

dist_v_tbl <- read.csv("Simulations/Products/crank_bs_ear1_sims_dist_v_dd_2x.csv")

aux0<- c(
    # " ",
    "HR",
    "SCPC",
    "C-SCPC",
    rep(
        c(
            "Kernel ",
            "Kernel ",
            "Kernel "#,
            # "Damian's"
            # "0.07"
        ),
        1
    ),
        rep(
        c(
            "HR",
            "Kernel ",
            "Kernel ",
            "Kernel ",
            "Damian's "
            # "0.07"
        ),
        18
    )
)

dist_tbl <- c(
    rep(
        "",
        3
    ),
    "0.015",
    "0.030",
    "0.080",
    gsub("NaN"," ", round(dist_v_tbl[[1]],3))
)
row_names_kernel <- paste0(aux0,dist_tbl)

row_names_ar1_kernel <- c(
    rep(
        "No splines",
        4+2
    ),
    rep(
        c(
        "4 splines",
        "12 splines",
        "24 splines",
        "36 splines",
        "48 splines",
        "64 splines",
        "88 splines",
        "102 splines",
        "124 splines"
        ),
        each=5,
        times =2
    )
)
aux <- c(
    rep(
        "",
        4+2
    ),
    rep(
        c(
            "Order 1",
            "Order 2"
        ),
        each = 45
    )
)
kernel_id<-paste0(aux, row_names_ar1_kernel)

table_kernel <- cbind(
    aux,
    row_names_ar1_kernel,
    row_names_kernel, 
    kernel_results, 
    kernel_id
)
names(table_kernel)<- c(
    "Order", "Qty.",
    "Estimator",
    "Cons.",
    "$\\beta_1$",
    "id"
)


row_names_ar1 <- c(
    "No splines",
    rep(
        c(
        "4 splines",
        "12 splines",
        "24 splines",
        "36 splines",
        "48 splines",
        "64 splines",
        "88 splines",
        "102 splines",
        "124 splines"
        ),
        2
    )
)
aux2 <- c(
    rep(
        "",
        1
    ),
    rep(
        c(
            "Order 1",
            "Order 2"
        ),
        each = 9
    )
)
ar1_id<-paste0(aux2, row_names_ar1)
table_ar1 <-cbind(
    ar1_id,
    e_ar1_tbl
)
names(table_ar1)<-c(
    "id",
    "Cons2",
    "$\\gamma_1$",
    "Hansen's",
    "Damian's",
    "Dropped"
)
table<-merge(table_kernel,table_ar1, by="id",sort=FALSE)


##


kernel_results <- read.csv("Simulations/Products/bic_spline_sims_res_dd_2x.csv")

e_ar1_tbl <- read.csv("Simulations/Products/bic_spline_sims_ar1_dd_2x.csv")


row_names_kernel<- c(
    # " ",
    "HR",
    "SCPC",
    "C-SCPC",
    rep(
        c(
            "Kernel 0.015",
            "Kernel 0.030",
            "Kernel 0.080"#,
            # "Damian's"
            # "0.07"
        ),
        1
    ),
        rep(
        c(
            "HR",
            "Kernel 0.015",
            "Kernel 0.030",
            "Kernel 0.080",
            "Damian's "
            # "0.07"
        ),
        2
    )
)

# dist_tbl <- c(
#     rep(
#         "",
#         3
#     ),
#     "0.015",
#     "0.030",
#     "0.080",
#     gsub("NaN"," ", round(dist_v_tbl[[1]],3))
# )
# row_names_kernel <- paste0(aux0,dist_tbl)

# row_names_ar1_kernel <- c(
#     rep(
#         "No splines",
#         4+2
#     ),
#     rep(
#         c(
#         "4 splines",
#         "12 splines",
#         "24 splines",
#         "36 splines",
#         "48 splines",
#         "64 splines",
#         "88 splines",
#         "102 splines",
#         "124 splines"
#         ),
#         each=1,
#         times =2
#     )
# )
aux <- c(
    rep(
        "",
        4+2
    ),
    rep(
        c(
            "Order 1",
            "Order 2"
        ),
        each = 5
    )
)
# kernel_id<-paste0(aux, row_names_ar1_kernel)

table_kernel <- cbind(
    aux,
    # row_names_ar1_kernel,
    row_names_kernel, 
    kernel_results#, 
    # kernel_id
)
names(table_kernel)<- c(
    "Splines",
    "Estimator",
    "Cons.",
    "$\\beta_1$"#,
    # "id"
)


# row_names_ar1 <- c(
#     "No splines",
#     rep(
#         c(
#         "4 splines",
#         "12 splines",
#         "24 splines",
#         "36 splines",
#         "48 splines",
#         "64 splines",
#         "88 splines",
#         "102 splines",
#         "124 splines"
#         ),
#         2
#     )
# )
aux2 <- c(
    rep(
        "",
        1
    ),
    rep(
        c(
            "Order 1",
            "Order 2"
        ),
        1
    )
)
# ar1_id<-paste0(aux2, row_names_ar1)
table_ar1 <-cbind(
    aux2,
    e_ar1_tbl
)
names(table_ar1)<-c(
    "Splines",
    "Cons2",
    "$\\gamma_1$",
    "Hansen's",
    "Damian's",
    "Dropped"
)
table<-merge(table_kernel,table_ar1, by="Splines",sort=FALSE)

##

kernel_results <- read.csv("Simulations/Products/grid_sims_res_dd_2x.csv")

e_ar1_tbl <- read.csv("Simulations/Products/grid_sims_ar1_dd_2x.csv")

dist_v_tbl <- read.csv("Simulations/Products/grid_sims_dist_v_dd_2x.csv")

aux0 <- c(
    # " ",
    # "SCPC",
    rep(
        c(
            "HR",
            "Kernel ",
            "Kernel ",
            "Kernel "
            # "0.07"
        ),
        1
    ),
    rep(
        c(
            "HR",
            "Kernel ",
            "Kernel ",
            "Kernel ",
            "Damian's "
            # "0.07"
        ),
        16
    )
)
dist_tbl <- c(
    rep(
        "",
        1
    ),
    "0.015",
    "0.030",
    "0.080",
    gsub("NaN"," ", round(dist_v_tbl[[1]],3))
)
row_names_kernel <- paste0(aux0,dist_tbl)
row_names_ar1_kernel <- c(
    rep(
        "No splines",
        4
    ),
    rep(
        c(
        "4 splines",
        "8 splines",
        "12 splines",
        "24 splines",
        "36 splines",
        "48 splines",
        "60 splines",
        "72 splines"#,
        # "88 splines",
        # "124 splines"
        ),
        each=5,
        times =2
    )
)
aux <- c(
    rep(
        "",
        4
    ),
    rep(
        c(
            "Order 1",
            "Order 2"
        ),
        each = 40
    )
)
kernel_id<-paste0(aux, row_names_ar1_kernel)

table_kernel <- cbind(
    aux,
    row_names_ar1_kernel,
    row_names_kernel, 
    kernel_results, 
    kernel_id
)
names(table_kernel)<- c(
    "Order", "Qty.",
    "Estimator",
    "Cons.",
    "$\\beta_1$",
    "id"
)

row_names_ar1 <- c(
    "No splines",
    rep(
        c(
        "4 splines",
        "8 splines",
        "12 splines",
        "24 splines",
        "36 splines",
        "48 splines",
        "60 splines",
        "72 splines"#,
        # "88 splines",
        # "124 splines"
        ),
        2
    )
)
aux2 <- c(
    rep(
        "",
        1
    ),
    rep(
        c(
            "Order 1",
            "Order 2"
        ),
        each = 8
    )
)
ar1_id<-paste0(aux2, row_names_ar1)
table_ar1 <-cbind(
    ar1_id,
    e_ar1_tbl
)
names(table_ar1)<-c(
    "id",
    "Cons2",
    "$\\gamma$",
    "Hansen's",
    "Damian's"
)

table<-merge(table_kernel,table_ar1, by="id",sort=FALSE)

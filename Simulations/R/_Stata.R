Stata <- function (src = stop("At least 'src' must be specified"), data.in = NULL, 
    data.out = FALSE, stata.path = getOption("RStata.StataPath", 
        stop("You need to set up a Stata path; ?chooseStataBin")), 
    stata.version = getOption("RStata.StataVersion", stop("You need to specify your Stata version")), 
    stata.echo = getOption("RStata.StataEcho", TRUE), ...) 
{
    if (!is.character(src)) 
        stop("src must be a character")
    if (!(is.null(data.in) | is.data.frame(data.in))) 
        stop("data.in must be NULL or a data.frame")
    if (!is.logical(data.out)) 
        stop("data.out must be logical")
    if (!is.numeric(stata.version)) 
        stop("stata.version must be numeric")
    if (!is.logical(stata.echo)) 
        stop("stata.echo must be logical")
    OS <- Sys.info()["sysname"]
    OS.type <- .Platform$OS.type
    SRC <- unlist(lapply(src, strsplit, "\n"))
    dataIn <- is.data.frame(data.in)
    dataOut <- data.out[1L]
    stataVersion <- stata.version[1L]
    stataEcho <- stata.echo[1L]
    if (OS %in% c("Windows","Darwin")) {
        winRStataLog <- "RStata.log"
        on.exit(unlink(winRStataLog))
    } #else if (OS %in% "Darwin") {
    #     macRStataLog <- "RStata.log"
    #     on.exit(unlink(macRStataLog))
    # }
    doFile <- "RStata.do"
    on.exit(unlink(doFile), add = TRUE)
    if (dataIn) {
        dtaInFile <- "RStataDataIn.dta"
        on.exit(unlink(dtaInFile), add = TRUE)
        foreign::write.dta(data.in, file = dtaInFile, version = if (stataVersion >= 
            7) 
            7L
        else 6L, ...)
    }
    if (dataOut) {
        dtaOutFile <- "RStataDataOut.dta"
        on.exit(unlink(dtaOutFile), add = TRUE)
    }
    if (file.exists(SRC[1L])) {
        SRC <- readLines(SRC[1L])
    }
    cut_me_here <- "RSTATA: cut me here"
    cut_me_comment <- paste0("/*", cut_me_here, "*/")
    SRC <- c({
        if (dataIn) sprintf("use %s", tools::file_path_sans_ext(dtaInFile)) else ""
    }, "capture noisily {", cut_me_comment, SRC, cut_me_comment, 
        "} /* end capture noisily */")
    SRC <- c("set more off", SRC)
    if (dataOut) {
        save_cmd <- sprintf("%s %s%s", if (stataVersion >= 13) 
            "saveold"
        else "save", tools::file_path_sans_ext(dtaOutFile), if (stataVersion >= 
            14) 
            ", version(12)"
        else "")
        SRC <- c(SRC, save_cmd)
    }
    SRC <- c(SRC, "exit, clear STATA")
    stata.batch <- if (OS == "Windows") {
            "/e"
        } else if (OS == "Darwin") {
            "-e"
        } else {
            ""
    }
    stataCmd <- paste(
        stata.path,
        stata.batch,
        "do",
        doFile
    ) #Adding -e  to run stata in batch  mode
    con <- file(doFile, "w")
    writeLines(SRC, con)
    close(con)
    rdl <- pipe(stataCmd, "r")
    stataLog <- readLines(rdl)
    close(rdl)
    if (stataEcho) {
        if (OS %in% c("Windows","Darwin")) {
            stataLog <- readLines(winRStataLog)
        } #else if (OS %in% "Darwin") {
            # stataLog <- readLines(macRStataLog)
        # }
        cutpoints <- grep(cut_me_here, stataLog)
        stataLog <- stataLog[seq.int(cutpoints[1] + 1, cutpoints[2] - 
            1)]
        cat(stataLog, sep = "\n")
    }
    if (dataOut) {
        res <- foreign::read.dta(dtaOutFile, ...)
        invisible(res)
    }
}

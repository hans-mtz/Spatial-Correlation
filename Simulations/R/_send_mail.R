# install.packages("mailR")
library(mailR)
send.mail(
    from = "hans.yoel.martinez@gmail.com",
    to = c("hans.yoel.martinez@gmail.com", "Hans UWO <hmarti33@uwo.ca>"),
    # replyTo = c("Reply to someone else <someone.else@gmail.com>"),
    subject = "Update",
    body = "Hi, here's an update \nBest,\nHans",
    smtp = list(host.name = "smtp.gmail.com", port = 465, user.name = "hans.yoel.martinez", passwd = "wtrfpvxjzxxasacr", ssl = TRUE),
    authenticate = TRUE,
    send = TRUE,
    attach.files = "./Document/SCPC.pdf",
)

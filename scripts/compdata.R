compdata <- sas7bdat::read.sas7bdat("/media/mtr/A5C2-009E/SCL/c9732_demographic.sas7bdat")
aedata <- sas7bdat::read.sas7bdat("E:/SCL/c9732_ae.sas7bdat")

sclc_demo$STATUS[sclc_demo$STATUS == 1] <- 0
sclc_demo$STATUS[sclc_demo$STATUS == 3] <- 1
sclc_demo$STATUS[sclc_demo$STATUS == 2] <- 1


compdata$TIMEDIFF <- compdata$OS_TIME - compdata$PFS_TIME
unique(compdata$TRT_ARM_LABEL)

compdata <- compdata[compdata$OS_TIME > 0, ]

compdata$PFS_STATUS[(compdata$STATUS == 1) & compdata$TIMEDIFF == 0] <- 0

compdata$TIMEDIFF[(compdata$PFS_STATUS == 1) & (compdata$TIMEDIFF == 0)] <- 1e-5

saveRDS(compdata, "data-raw/sclc.RDS")

compdata <- readRDS("data-raw/sclc.RDS")

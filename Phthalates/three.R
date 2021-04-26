# third trimester urinary phthalates 

# data: 3rdtrim_phthalates.xlsx, tides-fsex-3u-phthalates-plots.xlsx (df "three")

# AAPrev and AAPmonosg calculation

three$AAPrev_3U<-three$T3_mBP + (0.24*three$T3_miBP) + (0.26*three$T3_mBzP) + (0.61*three$T3_mEHP) + (0.024*three$T3_mEP)
three$lnAAPrev_3U<-log(three$AAPrev)

three$AAPmonosg_3U<- three$T3_MBP_SGadj + (0.24*three$T3_MiBP_SGadj) + (0.26*three$T3_MBzP_SGadj) + (0.61*three$T3_MEHP_SGadj) + (0.024*three$T3_MEP_SGadj)
three$lnAAPmonosg_3U<-log(three$AAPmonosg)

# log transformation of provided SG adjusted data

three$MBP_lnT3 <- log(three$T3_MBP_SGadj)
three$MBzP_lnT3 <- log(three$T3_MBzP_SGadj)
three$MCNP_lnT3 <- log(three$T3_MCNP_SGadj)
three$MCOP_lnT3 <- log(three$T3_MCOP_SGadj)
three$MCPP_lnT3 <- log(three$T3_MCPP_SGadj)
three$MECPP_lnT3 <- log(three$T3_MECPP_SGadj)
three$MECPTP_lnT3 <- log(three$T3_MECPTP_SGadj)
three$MEHHP_lnT3 <- log(three$T3_MEHHP_SGadj)
three$MEHP_lnT3 <- log(three$T3_MEHP_SGadj)
three$MEOHP_lnT3 <- log(three$T3_MEOHP_SGadj)
three$MEP_lnT3 <- log(three$T3_MEP_SGadj)
three$MHiBP_lnT3 <- log(three$T3_MHiBP_SGadj)
three$MiBP_lnT3 <- log(three$T3_MiBP_SGadj)
three$MONP_lnT3 <- log(three$T3_MONP_SGadj)

# saved as tides-fsex-3u-phthalates-plots.xlsx
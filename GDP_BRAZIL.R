###########Econometrics II#################
###########Cristhian Ortiz#################

# Clear plots
if(!is.null(dev.list())) dev.off()
#Tamanio de Memoria
memory.size(max=F)
# Garbage collection
gc()
# Clear console
cat("\014") 
# Clean workspace
rm(list=ls())
# Restart R
.rs.restartR()

#Obteniendo la data del Banco mundial

library("WDI") 
new_wdi_cache <- WDIcache() 

wdi_dat <- WDI(indicator = c("PA.NUS.FCRF", "NY.GDP.MKTP.KN", "NE.CON.PRVT.KN", 
                             "NE.CON.GOVT.KN", "NE.EXP.GNFS.KN", "NE.IMP.GNFS.KN", 
                             "NE.GDI.FTOT.KN", "NE.GDI.STKB.KN"), 
               country = "BRA", start = 1975, end = 2013, extra = FALSE)

wdi_dat.ts <- ts(wdi_dat)

#	GDP (constant LCU)
gdp <- wdi_dat.ts[,"NY.GDP.MKTP.KN"] 
#	Household final consumption expenditure (constant LCU)
con <- wdi_dat.ts[,"NE.CON.PRVT.KN"]
#General government final consumption expenditure (constant LCU)
gov <- wdi_dat.ts[,"NE.CON.GOVT.KN"]
#Exports of goods and services (constant LCU)
exp <- wdi_dat.ts[,"NE.EXP.GNFS.KN"]
#Imports of goods and services (constant LCU)
imp <- wdi_dat.ts[,"NE.IMP.GNFS.KN"]
#Changes in inventories (constant LCU)
invch <- wdi_dat.ts[,"NE.GDI.STKB.KN"]
#Gross fixed capital formation (constant LCU)
inv <- wdi_dat.ts[,"NE.GDI.FTOT.KN"]

#Aplicando la tasa de cambio con diff-log
g_gdp <- diff(log(gdp))
g_con <- diff(log(con))
g_gov <- diff(log(gov))
g_exp <- diff(log(exp))
g_imp <- diff(log(imp))
g_invch <- diff(log(invch))
g_inv <- diff(log(inv))

#Suavizamiento solo con logaritmo
lgdp <- log(gdp)
lcon <- log(con)
lgov <- log(gov)
lexp <- log(exp)
limp <- log(imp)
linvch <- log(invch)
linv <- log(inv)

#Aplicando los lags

lag_gdp <- lag(lgdp,1)
lag_con <- lag(lcon,1)
lag_gov <- lag(lgov,1)
lag_exp <- lag(lexp,1)
lag_imp <- lag(limp,1)
lag_invch <- lag(linvch,1)
lag_inv <- lag(linv,1)

library(tseries)
#Aplicaremos la prueba de cointegracion de Phillips Oularis

cointeg_test <- as.matrix(cbind(g_gdp, g_con, g_gov, g_exp ))
po.test(cointeg_test)
#p-value = 0.01 las series son cointegradas

cointeg_test2 <- as.matrix(cbind(g_gdp, g_imp, g_invch, g_inv))
po.test(cointeg_test2)
#p-value = 0.02 las series son cointegradas

cointeg_test3 <- as.matrix(cbind(g_gdp, lag_gdp, lag_con, lag_gov, lag_inv))
po.test(cointeg_test3)
#p-value = 0.02 las series son cointegradas

cointeg_test4 <- as.matrix(cbind(g_gdp, lag_exp, lag_imp, lag_invch))
po.test(cointeg_test4)
#p-value = 0.01 las series son cointegradas

#Analizando la Estacionariedad en el Residuo

aux3.1 <- lm(g_gdp~g_con+g_gov+g_inv+g_exp+g_imp+g_invch+
               lag_gdp[-1]+lag_con[-1]+lag_gov[-1]+lag_exp[-1]+lag_imp[-1]+
               lag_invch[-1]+lag_inv[-1])

aux3.2 <- resid(aux3.1)
pp.test(aux3.2)
#p-value=0.01 el residuo es estacionario

m_aux <- lm(g_gdp ~ lag_gdp[-1]+lag_con[-1]+lag_gov[-1]+lag_inv[-1]+lag_exp[-1]+
              lag_imp[-1]+lag_invch[-1]+g_con+g_gov+g_exp+g_imp+g_invch+g_inv)
stage_1 <- resid(m_aux)
pp.test(stage_1)
#p-value=0.01 el residuo es estacionario

aux <- lm(g_gdp ~ lag_gdp[-1]+lag_con[-1]+lag_gov[-1]+lag_inv[-1]+lag_exp[-1]+
            lag_imp[-1]+g_con+g_gov+g_exp+g_imp+g_inv)

aux2 <- resid(aux)
pp.test(aux2)
#p-value=0.01 el residuo es estacionario

#Replicando el paper obteniendo lag 2

diff_ggdp <- diff(g_gdp)
diff_gcon <- diff(g_con)
diff_ggov <- diff(g_gov)
diff_gexp <- diff(g_exp)
diff_gimp <- diff(g_imp)
diff_ginvch <- diff(g_invch)
diff_ginv <- diff(g_inv)

l2_gdp <- lag(g_gdp,1)
l2_con <- lag(g_con,1)
l2_gov <- lag(g_gov,1)
l2_exp <- lag(g_exp,1)
l2_imp <- lag(g_imp,1)
l2_invch <- lag(g_invch,1)
l2_inv <- lag(g_inv,1)

model_1 <- lm(diff_ggdp ~ l2_gdp[-1] + l2_con[-1] + l2_gov[-1] + l2_inv[-1] +
                l2_exp[-1] + l2_imp[-1] + l2_invch[-1] + diff_gcon + diff_ggov + 
                diff_gexp + diff_gimp + diff_ginvch + diff_ginv + aux2[-1])
summary(model_1)

model_2 <- lm(diff_ggdp ~ l2_gdp[-1] + l2_con[-1] + l2_gov[-1] + l2_inv[-1] +
                l2_exp[-1] + l2_imp[-1] + diff_gcon + diff_ggov + 
                diff_gexp + diff_gimp + diff_ginv + aux2[-1])
summary(model_2)

model_3 <- lm(diff_ggdp ~ diff_gcon + diff_ggov + diff_gexp + diff_gimp 
              + diff_ginv + aux2[-1])
summary(model_3)

model_def <- lm(diff_ggdp ~ l2_gdp[-1] + l2_con[-1] + l2_gov[-1] + l2_inv[-1] +
                  l2_exp[-1] + l2_imp[-1] + diff_gcon + diff_ggov + diff_gexp + 
                  diff_gimp  + diff_ginv + aux2[-1])
summary(model_def)

gr_gdp <- predict(model_def)

growth_rate_gdp <- ts(gr_gdp, start=c(1975,1), end=c(2013,1),
                      frequency=1)
obs_gr_gdp <- ts(diff_ggdp, start=c(1975,1), end=c(2013,1),
                 frequency=1)

grate_gdp_brasil <- cbind(growth_rate_gdp, obs_gr_gdp)

plot.ts(grate_gdp_brasil, plot.type=c("single"), 
        main="Brazil annual GDP Growth 1975-2013",xlab="Year", ylab="Diff Tasas de crecimiento", col=c("black", "blue"),lty=2:1, lwd=1:1)
legend("topright", legend = c("ECM estimation","Observed"),col=c("black", "blue"),lty=2:1, cex=.65)


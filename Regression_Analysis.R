##----Import the libraries
library(vars)
library(tsDyn)
library(tseries)
library(urca)
library(dynlm)

########################################################################################
#### PREPARATION STEPS
##----To Import the Cleaned Dataset in Python
df <- read.csv('./millennium/cleaned_sample.csv')
vecm_df_1 = data.frame(df['Adj_Close_log'], df['Sig_log'])

##----TO perform Augmented Dicky FUller test for the stationarity information
adf_ac = adf.test(vecm_df_1$Adj_Close_log)
adf_ac_diff = adf.test(diff(vecm_df_1$Adj_Close_log))
adf_sig = adf.test(vecm_df_1$Sig_log)
adf_sig_diff = adf.test(diff(vecm_df_1$Sig_log))
#--Result: the log value of adjusted close and that of the signal as well as their first order differences are all staionary

##----To estimate the optimal lag length, firs estimate a VAR model
var_md = VARselect(vecm_df_1, lag.max = 20, type = c("const"), season = NULL)
var_md
#--Result: Optimal lag length being 3 or 17

##----Implementing Johnason test to investigating the cointegration relationships
jt = ca.jo(vecm_df_1, ecdet = "const", type = "eigen", K = 17)
summary(jt)
#--Result: 'Adj_close_log' and 'Sig_log' are cointegrated
########################################################################################

########################################################################################
#### CAUSALITY TEST FOLLOWING TODA-YAMAMTO PRCEDURE
##----We set up 2 VAR models with lag = 3 and lag = 17 and do serial tests on them
var_3<- VAR(vecm_df_1, p=3, type = 'const')
serial.test(var_3, type = "ES")
var_17<- VAR(vecm_df_1, p=17, type = 'const')
serial.test(var_17, type = "ES")
#--Result: It turns out that the lag length of 17 result in much smaller serial correlation in the residuals
##----Then we set up a VAR model with the lag length of 18, and do the subsequent causality test
var_18 <- VAR(vecm_df_1, p=18, type = 'const')
causality(var_18, cause = 'Adj_Close_log')
causality(var_18, cause = 'Sig_log')
causality(var_18, cause = 'Signal')
#--Result: Results clearly show the fact that the causation relationship between 'Adj_Close_log' and 'Sig_log'
#--is one-way, in particular, it is that 'Sig_log' Granger causes 'Adj_Close_log'.
########################################################################################

########################################################################################
#### RUNNING VECM MODEL
##----Implementing VECM model with the lag to be 17-1=16
vecm = VECM(vecm_df_1, lag = 16, LRinclude = "const", r = 1, estim = 'ML')
#vecm_3_jo = cajorls(jt_3, r=1) ## This is exactly the same model as the above one
summary(vecm)
trsfm_var = vec2var(jt, r=1)
##----Testing the autocorrelation, heteroscedasticity and normality of the residuals
serial.test(trsfm_var, type = "ES")
arch.test(trsfm_var)
normality.test(trsfm_var, multivariate.only = TRUE)
##----There is no autocorrelation problem found in this vecm model. But heteroskedasticity is detected.
##----Heteroskedasticity won't affect the unbiasness or consistency of the OLS, but causes inefficiency as well as the statistics.
##----We thus use IRF and FEVD analyze the estimated VECM
plot(irf(trsfm_var, n.ahead = 20), plot.type = "single")
plot(fevd(trsfm_var, n.ahead = 20), xlab = "Period")

##----We use the NewyWest-corrected variance for the correction of the statistics of the OLS.
adj_log <- ts(vecm_df_1$Adj_Close_log)
sig_log <- ts(vecm_df_1$Sig_log)
vecm_ect = ts(adj_log - 0.8744696*sig_log - 2.451243)

VECM_EQ1 <- dynlm(d(adj_log) ~ L(d(adj_log), 1:16) + L(d(sig_log), 1:16) + L(vecm_ect))
VECM_EQ2 <- dynlm(d(sig_log) ~ L(d(adj_log), 1:16) + L(d(sig_log), 1:16) + L(vecm_ect))

coeftest(VECM_EQ1, vcov. = NeweyWest(VECM_EQ1, prewhite = F, adjust = T))
coeftest(VECM_EQ2, vcov. = NeweyWest(VECM_EQ2, prewhite = F, adjust = T))
#--Result: The statistics of the coefficients change little.

########################################################################################
#### Doing the spliting test
df_train <- read.csv('./millennium/cleaned_sample_train.csv')
df_test <- read.csv('./millennium/cleaned_sample_test.csv')
train_df = data.frame(df_train['Adj_Close_log'], df_train['Sig_log'])
test_df = data.frame(df_test['Adj_Close_log'], df_test['Sig_log'])


#vecm = VECM(train_df, lag = 16, LRinclude = "const", r = 1, estim = 'ML')
jt_train = ca.jo(train_df, ecdet = "const", type = "eigen", K = 17)
var = vec2var(jt_train, r=1)

forecast = predict(var, n.ahead = 17)
forecast
fanchart(forecast, plot.type = "single", xlim=c(1000,1038), xlab = "Period")


fit <- arima(train_df$Adj_Close_log)
predict(fit, n.ahead = 10)

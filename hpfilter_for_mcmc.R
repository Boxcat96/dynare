library(mFilter)
data <- new.env()

# 推計開始時期と終了時期
date.start <- "1997-01-01"
date.end <- "2007-12-31"

tickers <- c("GDPC1", "PCECC96", "GPDIC1", "LES1252881600Q", "LREM64TTUSQ156S")

# FREDのデータベースにアクセス
library("quantmod")
getSymbols( tickers
            , src = "FRED"
            , from = date.start 
            , to = date.end  
            , env = data
)

#階層分け
dtx1 <- data$GDPC1
x1 <- dtx1[paste(date.start,date.end,sep="/")]
dtx2 <- data$PCECC96
x2 <- dtx2[paste(date.start,date.end,sep="/")]
dtx3 <- data$GPDIC1
x3 <- dtx3[paste(date.start,date.end,sep="/")]
dtx4 <- data$LES1252881600Q
x4 <- dtx4[paste(date.start,date.end,sep="/")]
dtx5 <- data$LREM64TTUSQ156S
x5 <- dtx5[paste(date.start,date.end,sep="/")]

#HPフィルター
hpf1 <- hpfilter(log(x1),freq = 1600)
hpf2 <- hpfilter(log(x2),freq = 1600)
hpf3 <- hpfilter(log(x3),freq = 1600)
hpf4 <- hpfilter(log(x4),freq = 1600)
hpf5 <- hpfilter(log(x5),freq = 1600)


# 結果
out <- xts(cbind(hpf1$cycle, hpf2$cycle, hpf3$cycle, hpf4$cycle, hpf5$cycle)
           , index(x1, x2, x3, x4, x5))
colnames(out) <- c("y_obs", "c_obs", "i_obs", "w_obs", "n_obs")

# 結果の表示
head(out)

# 結果の出力
write.csv(out, "C:/cat/mcmc.csv", row.names = F)

# 結果の描画
par(mfrow = c(2, 1), mar = c(3, 2, 2, 1))
colnames(out) <- c("real GDP", "real consumption", "real investment", 
                   "real wage", "employment rate")
plot(out[,"real GDP"], t= "n", main = "real GDP", col = "steelblue")
plot(out[,"real consumption"], t= "n", main = "real consumption", col = "steelblue")
plot(out[,"real investment"], t= "n", main = "real investment", col = "steelblue")
plot(out[,"real wage"], t= "n", main = "real wage", col = "steelblue")
plot(out[,"employment rate"], t= "n", main = "employment rate", col = "steelblue")


mvSUSY
----

Multivariate Surrogate Synchrony estimates the synchrony within datasets that contain more than two time series. Multivariate SUSY was adapted and extended from SUCO and SUSY with respect to implementing the surrogate controls. 'mvSUSY' works as described in Meier & Tschacher (2021).

[R package website](https://wtschacher.github.io/mvSUSY/)

----

Installation
----

```r
#install.packages("mvSUSY") ## not yet on CRAN

## development version
install.packages("mvSUSY", repos=c("https://wtschacher.github.io/mvSUSY/","https://cloud.r-project.org"))
```

Usage
----

Note that the following example assumes that the source data are in a flat file and it has particular structure (column names in first row, whitespace as field separator). If you do not have such, then use the command in the comment below to mockup random data.

```r
library(mvSUSY)

## read in data from a flat file
data = read.csv(file.choose(), header=TRUE, sep=" ", na.strings=".")

## mockup random data if needed
#data = as.data.frame(replicate(5, sample(10, 500, TRUE)))

## compute mvSUSY using 'lambda_max' method
res = mvsusy(data, segment=10, Hz=10, max_pseudo=100)
res

## plot
plot(res, type="eigenvalue")
plot(res, type="density")
plot(res, type="free scale")
plot(res, type="segment-wise")
plot(res, type="time series")

## compute mvSUSY using 'omega' method
res = mvsusy(data, segment=10, Hz=10, max_pseudo=100, method="omega")
res

plot(res, type="density")
plot(res, type="free scale")
plot(res, type="segment-wise")
plot(res, type="time series")

## export to flat file via data.frame and write.csv
df = as.data.frame(res)
df
```

[`mvsusy` function manual](https://wtschacher.github.io/mvSUSY/reference/mvsusy.html)

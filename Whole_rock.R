library(circular)
library(factoextra)
library(ggpubr)


whole_rock_dataset = read.csv(file.choose())
dim(whole_rock_dataset)
watson.test(whole_rock_dataset$latitude, alpha = 0.01, dist = "vonmises")
watson.test(whole_rock_dataset$longitude, alpha = 0.01, dist = "vonmises")
watson.test(whole_rock_dataset$latitude, alpha = 0.01, dist = "uniform")
watson.test(whole_rock_dataset$longitude, alpha = 0.01, dist = "uniform")


#View(whole_rock_dataset)
wrd = whole_rock_dataset[,-c(1,2,5:57)]

head(wrd)

wrd = wrd[,-c(3:22)]
wrd1 = wrd[,-c(106:117)]
wrd1[is.na(wrd1)] <- 0

sum(is.na(wrd1))

constant_vars <- sapply(wrd1, function(x) max(x) == 0 && min(x) == 0)

# Subset your data frame to keep only non-constant variables
wrd1 <- wrd1[, !constant_vars]

summary(wrd1$ag_ppm)
wrd2 = wrd1

wss = (nrow(wrd2)-1)*sum(apply(wrd2 , 2, var))
for (i in 2:13)
{
  wss[i] =
    sum(kmeans(wrd2 , centers = i)$withinss)
}
plot(1:13 , wss , type = "l" , 
     xlab = "number of clusters" ,
     ylab = "within group of sum of  square")

fit_c = kmeans(wrd2 , 4)

fviz_cluster(fit_c, data = wrd2,
             palette = c("green", "#00AFBB", "#E7B800", "red"), 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw()
)


fit_c$size
#View(wrd1)

f_data = cbind(whole_rock_dataset, fit_c$cluster)
sum(is.na(f_data$rock_origin))

unique(f_data$rock_origin)

library(skmeans)

clus_data = cbind(whole_rock_dataset$latitude, whole_rock_dataset$longitude)

sph_cluster = skmeans(clus_data, 4)


f_data_2 = cbind(whole_rock_dataset , sph_cluster$cluster)

partitions = cbind(whole_rock_dataset$rock_type,whole_rock_dataset$country,whole_rock_dataset$rock_group, whole_rock_dataset$rock_origin, f_data$`fit_c$cluster` , f_data_2$`sph_cluster$cluster`)
head(partitions)
partitions
unique(whole_rock_dataset$rock_type)

library(Directional)
library(movMF)

d = cbind(whole_rock_dataset$latitude, whole_rock_dataset$longitude)
Evmf <- function(K){
  movMF(d, k= K, control = list(nruns = 20))
}

set.seed(123)
Esd = lapply(1:10, Evmf)
Esd
sapply(Esd, BIC)



vmf_density_grid = 
  function(u, ngrid = 100) {
    # Translate to (0,180) and (0,360)
    u[,1] <- u[,1] + 90
    u[,2] <- u[,2] + 180
    res <- vmf.kerncontour(u, thumb = "none", den.ret = T, full = T,
                           ngrid = ngrid)
    
    # Translate back to (-90, 90) and (-180, 180) and create a grid of
    # coordinates
    ret <- expand.grid(Lat = res$lat - 90, Long = res$long - 180)
    ret$Density <- c(res$den)
    ret
  }


wholer.dens = vmf_density_grid(d, ngrid = 300)

library(ggplot2)
library(maps)

world <- map_data("world")
whole_rock_landing <- ggplot() +
  geom_map(data = world, map = world,
           mapping = aes(map_id = region),
           color = "black", fill = "white") +
  geom_point(data = whole_rock_dataset,
             mapping = aes(x = whole_rock_dataset$longitude,
                           y = whole_rock_dataset$latitude),
             color = "red", alpha = .5, size = 1, stroke = 0.1) +
  geom_density_2d(data = whole_rock_dataset,
                  aes(x = whole_rock_dataset$longitude,
                      y = whole_rock_dataset$latitude),
                  color = "#b266ff", alpha = 2) +
  geom_contour(data = wholer.dens, aes(x=Long, y=Lat, z=Density),
               color = "blue") 

plot(whole_rock_landing)



whole_rock_landing <- ggplot() +
  geom_map(data = world, map = world,
           mapping = aes(map_id = region),
           color = "black", fill = "white") +
  geom_point(data = whole_rock_dataset,
             mapping = aes(x = whole_rock_dataset$longitude,
                           y = whole_rock_dataset$latitude),
             color = "red", alpha = .5, size = 1, stroke = 0.1) +
  geom_density_2d(data = whole_rock_dataset,
                  aes(x = whole_rock_dataset$longitude,
                      y = whole_rock_dataset$latitude),
                  color = "#b266ff", alpha = 2) +
  geom_contour(data = wholer.dens, aes(x=Long, y=Lat, z=Density),
               color = "blue") +
coord_map("orthographic", orientation = c(0, 110, 0)) +
  scale_x_continuous(breaks = seq(-180, 180, 20)) +
  scale_y_continuous(breaks = seq(-90, 90, 45)) +
  ggtitle("Orthographic Projection of Spherical Density", "Top / Front View") +
  xlab("") +
  ylab("") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.ontop = TRUE,
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid = element_line(color = "black" ),
        panel.background = element_rect(fill = NA))
whole_rock_landing


library(circular)
library(CircStats)
x= whole_rock_dataset$latitude
mu = mean.circular(whole_rock_dataset$latitude)
kappaa = est.kappa(whole_rock_dataset$latitude)
y <- rvonmises(n=1000, mu, kappaa)
resx <- density(x, bw=25)
resy <- density.circular(y, bw=25)
pp=plot(density.circular(whole_rock_dataset$latitude , bw = 10) ,points.plot=F, pty = 10, lty=1, lwd=2, col="blue",xlim=c(-1.2,1), ylim=c(-1.1, 1.2), main="Comparison of  estimated sample density\n with fitted Von Mises\n for whole rock Latitude data", cex.main=0.25)
lines(resy, points.plot=F, col="red", points.col=2,lwd=2, lty=4, plot.info=pp)

legend("topleft", legend=c("estimated \n kernel density", "von mises density"),
       col=c("blue", "red"), lty=c(1,4),lwd=c(2,2), cex=0.59,
       box.lty=0)


x= whole_rock_dataset$longitude
mu = mean.circular(whole_rock_dataset$longitude)
y <- rvonmises(n=1000, mu, kappaa)
resx <- density(x, bw=25)
resy <- density.circular(y, bw=25)
pp=plot(density.circular(whole_rock_dataset$longitude , bw = 10) ,points.plot=F, pty = 10, lty=1, lwd=2, col="blue",xlim=c(-1.2,1), ylim=c(-1.1, 1.2), main="Comparison of  estimated sample density\n with fitted Von Mises\n for whole rock Longitude data", cex.main=0.25)
lines(resy, points.plot=F, col="red", points.col=2,lwd=2, lty=4, plot.info=pp)

legend("topleft", legend=c("estimated \n kernel density", "von mises density"),
       col=c("blue", "red"), lty=c(1,4),lwd=c(2,2), cex=0.59,
       box.lty=0)

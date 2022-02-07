require(ggplot2)
require(grid)
require(rworldmap)
require(tidyr)

covid <- read.csv(file = file.choose(), header = TRUE, sep = ",")
covid.c <- covid[, c(3, 4, 12)]
colnames(covid.c) <- c("Entity", "Date", "new")
covid.c.r <- subset(covid.c, Entity == "Belgium" |  Entity == "Germany"|  Entity == "Italy"| Entity == "France" |
                      Entity == "Spain"| Entity == "Austria"| Entity == "Greece" | 
                      Entity == "Netherlands" |Entity == "United Kingdom"| Entity == "Ireland"|
                      Entity == "Switzerland" | Entity == "Sweden"|Entity == "Norway"|
                      Entity == "Denmark"| Entity == "Poland" | Entity == "Portugal" | Entity == "Luxembourg"|
                      Entity == "Romania"| Entity == "Slovakia" | Entity == "Bulgaria" | Entity == "Serbia"|
                      Entity == "Hungary" | Entity == "Slovenia" | Entity == "Lithuania" | Entity == "Finland"|
                      Entity == "Latvia" | Entity == "Estonia" | Entity == "Malta" | Entity == "Czechia" | 
                      Entity == "Iceland"| Entity == "Croatia" | Entity == "Bosnia and Herzegovina" | 
                      Entity == "Cyprus" | Entity == "Bulgaria" | Entity == "Albania" )

library(tidyr)
covid.c.r <- spread(covid.c.r, Date, new)
Y = covid.c.r[, 2:424]
Y[Y<0] <- 0 

# fit.q1c <- quan.smsp(Y, alpha = 0.1, interval = c(1e-07, 9e-07))
# fit.q2c <- quan.smsp(Y, alpha = 0.3, interval = c(2e-08, 2e-07))
# fit.q3c <- quan.smsp(Y, alpha = 0.5, interval = c(5e-09, 5e-08))
# fit.q4c <- quan.smsp(Y, alpha = 0.7,  interval = c(8e-10, 5e-09))
# fit.q5c <- quan.smsp(Y, alpha = 0.9,  interval  = c(6e-10, 6e-09))

fit.q1c <- quan.smsp(Y, alpha = 0.1)
fit.q2c <- quan.smsp(Y, alpha = 0.3)
fit.q3c <- quan.smsp(Y, alpha = 0.5)
fit.q4c <- quan.smsp(Y, alpha = 0.7)
fit.q5c <- quan.smsp(Y, alpha = 0.9)

par(mar = c(3.5, 5.7, 3.1, 1.1))
par(mgp = c(3.8, 1, 0))
matplot(t(Y), type = "l", col = "gray", lwd = 3, lty = 1, pch = 20, cex = 1.5,
        ylim = c(0, 2500), cex.lab = 2.5, cex.axis = 2.5, xaxt = "n", xlab = "", ylab = "New cases (per million)") ;
grid()
library(reshape2)
library(ggplot2)
data <- data.frame(t(Y))
data$id <- 1:nrow(data)
#reshape to long format
plot_data <- melt(data, id.var="id" )
#plot
p <- ggplot(plot_data, aes(x=id, y=value, group=variable, colour=variable)) + geom_line(col = "gray", size = 1.3) +labs(x = "", y = "New cases (per million)")
p <- p + scale_x_continuous(breaks = c(60, 150, 240, 330, 420),
                            labels= c("Mar", "June", "Sep", "Dec", "Mar"))
p <- p + theme_bw(base_size = 40)
p




plot(fit.q1c$mu, type = "l", lwd = 3, col = rgb(0, 1, 0), ylim = c(0, 770), xlab = "",
     ylab = "New cases (per million)",  cex.lab = 2.3, cex.axis = 2.3, xaxt = "n", lty = 2)
lines(fit.q2c$mu, type = "l", lwd = 3, col = rgb(0.2470588,  0.7490196, 0.2470588), lty = 3)
lines(fit.q3c$mu, type = "l", lwd = 3, col = rgb(0.4980392, 0.4980392, 0.4980392), lty = 1)
lines(fit.q4c$mu, type = "l", lwd = 3, col = rgb(0.7490196, 0.2470588, 0.7490196), lty = 4)
lines(fit.q5c$mu, type = "l", lwd = 3, col = rgb(1, 0, 1), lty = 5)
grid()
axis(1, at=c(15, 45, 75, 105, 135, 165, 195, 225, 255, 285, 315, 345, 375, 405, 426), 
     labels= c("Jan", "Feb.", "Mar", "Apr", "May", "June", "July", "Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar"), 
     cex = 2.5, cex.axis = 2.5, gap.axis = 0.1)

data <- data.frame(t = 1:ncol(Y), fit1 = fit.q1c$mu, fit2 = fit.q2c$mu, fit3 = fit.q3c$mu, 
                   fit4 = fit.q4c$mu, fit5 = fit.q5c$mu)  
library(ggplot2)
p <- ggplot(data = data, aes(x = t, y = fit1)) + geom_line(size = 1.4, linetype = "dashed",  col = rgb(0, 1, 0)) +  ylim(0, 800) + labs(x = "", 
                                                                                                                 y = "New cases (per million)")
p <- p + scale_x_continuous(breaks = c(60, 150, 240, 330, 420),
                            labels= c("Mar", "June", "Sep", "Dec", "Mar"))
p <- p + geom_line(aes(x = t, y = fit2), size = 1.4, linetype = "dotted", col =  rgb(0.2470588,  0.7490196, 0.2470588))
p <- p + geom_line(aes(x = t, y = fit3), size = 1.4, linetype = "solid", col =  rgb(0.4980392, 0.4980392, 0.4980392))
p <- p +  geom_line(aes(x = t, y = fit4), size = 1.4, linetype = "dotdash", col =   rgb(0.7490196, 0.2470588, 0.7490196))
p <- p +  geom_line(aes(x = t, y = fit5), size = 1.4, linetype = "longdash", col =   rgb(1, 0, 1))
p <- p + theme_bw(base_size = 40)
p

# Deaths
covid.d <- covid[, c(3, 4, 15)]
colnames(covid.d) <- c("Entity", "Date", "new")

covid.d.r <- subset(covid.d, Entity == "Belgium" |  Entity == "Germany"|  Entity == "Italy"| Entity == "France" |
                      Entity == "Spain"| Entity == "Austria"| Entity == "Greece" | 
                      Entity == "Netherlands" |Entity == "United Kingdom"| Entity == "Ireland"|
                      Entity == "Switzerland" | Entity == "Sweden"|Entity == "Norway"|
                      Entity == "Denmark"| Entity == "Poland" | Entity == "Portugal" | Entity == "Luxembourg"|
                      Entity == "Romania"| Entity == "Slovakia" | Entity == "Bulgaria" | Entity == "Serbia"|
                      Entity == "Hungary" | Entity == "Slovenia" | Entity == "Lithuania" | Entity == "Finland"|
                      Entity == "Latvia" | Entity == "Estonia" | Entity == "Malta" | Entity == "Czechia" | 
                      Entity == "Iceland"| Entity == "Croatia" | Entity == "Bosnia and Herzegovina" | 
                      Entity == "Cyprus" | Entity == "Bulgaria" | Entity == "Albania" )

covid.d.r <- spread(covid.d.r, Date, new)
Y = covid.d.r[, 2:424]
Y[Y<0] <- 0 
par(mar = c(3.5, 5.7, 3.1, 1.1))
par(mgp = c(3.8, 1, 0))
matplot(t(Y), type = "l", col = "gray", lwd = 3, lty = 1, pch = 20, cex = 1.5,
        ylim = c(0, 50), cex.lab = 2.5, cex.axis = 2.5, xaxt = "n", xlab = "", ylab = "New cases (per million)") ;
grid()
data <- data.frame(t(Y))
data$id <- 1:nrow(data)
#reshape to long format
plot_data <- melt(data, id.var="id" )
#plot
p <- ggplot(plot_data, aes(x=id, y=value, group=variable, colour=variable)) + geom_line(col = "gray", size = 1.3) +labs(x = "", y = "New deaths (per million)")
p <- p + scale_x_continuous(breaks = c(60, 150, 240, 330, 420),
                            labels= c("Mar", "June", "Sep", "Dec", "Mar"))
p <- p + theme_bw(base_size = 40)
p



fit.q1d <- quan.smsp(Y, alpha = 0.1, interval = c(1e-01, 1.1))
fit.q2d  <- quan.smsp(Y, alpha = 0.3, interval = c(2e-06, 2e-05))
fit.q3d  <- quan.smsp(Y, alpha = 0.5,  interval = c(6e-07, 6e-06))
fit.q4d  <- quan.smsp(Y, alpha = 0.7, interval = c(1e-07, 1e-06))
fit.q5d  <- quan.smsp(Y, alpha = 0.9, interval  = c(6e-07, 6e-08))

# fit.q1d <- quan.smsp(Y, alpha = 0.1)
# fit.q2d  <- quan.smsp(Y, alpha = 0.3)
# fit.q3d  <- quan.smsp(Y, alpha = 0.5)
# fit.q4d  <- quan.smsp(Y, alpha = 0.7)
# fit.q5d  <- quan.smsp(Y, alpha = 0.9)

par(mar = c(3.5, 5.7, 3.1, 1.1))
par(mgp = c(3.8, 1, 0))
matplot(t(Y), type = "l", col = "gray", lwd = 3, lty = 1, pch = 20, cex = 1.5,
        ylim = c(0, 50), cex.lab = 2.5, cex.axis = 2.5, xaxt = "n", xlab = "", ylab = "New deaths (per million)") ;
grid()

plot(fit.q1d$mu, type = "l", lwd = 3, col = rgb(0, 1, 0), ylim = c(0, 17), xlab = "",
     ylab = "New deaths (per million)",  cex.lab = 2.3, cex.axis = 2.3, xaxt = "n", lty = 2)
lines(fit.q2d$mu, type = "l", lwd = 3, col = rgb(0.2470588,  0.7490196, 0.2470588), lty = 3)
lines(fit.q3d$mu, type = "l", lwd = 3, col = rgb(0.4980392, 0.4980392, 0.4980392), lty = 1)
lines(fit.q4d$mu, type = "l", lwd = 3, col = rgb(0.7490196, 0.2470588, 0.7490196), lty = 4)
lines(fit.q5d$mu, type = "l", lwd = 3, col = rgb(1, 0, 1), lty = 5)
grid()
axis(1, at=c(15, 45, 75, 105, 135, 165, 195, 225, 255, 285, 315, 345, 375, 405, 426), 
     labels= c("Jan", "Feb.", "Mar", "Apr", "May", "June", "July", "Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar"), 
     cex = 2.5, cex.axis = 2.5, gap.axis = 0.1)

data <- data.frame(t = 1:ncol(Y), fit1 = fit.q1d$mu, fit2 = fit.q2d$mu, fit3 = fit.q3d$mu, 
                   fit4 = fit.q4d$mu, fit5 = fit.q5d$mu)  
library(ggplot2)
p <- ggplot(data = data, aes(x = t, y = fit1)) + geom_line(size = 1.4, linetype = "dashed",  col = rgb(0, 1, 0)) +  ylim(0, 17.5) + labs(x = "", 
                                                                                                                                        y = "New cases (per million)")
p <- p + scale_x_continuous(breaks = c(60, 150, 240, 330, 420),
                            labels= c("Mar", "June", "Sep", "Dec", "Mar"))
p <- p + geom_line(aes(x = t, y = fit2), size = 1.4, linetype = "dotted", col =  rgb(0.2470588,  0.7490196, 0.2470588))
p <- p + geom_line(aes(x = t, y = fit3), size = 1.4, linetype = "solid", col =  rgb(0.4980392, 0.4980392, 0.4980392))
p <- p +  geom_line(aes(x = t, y = fit4), size = 1.4, linetype = "dotdash", col =   rgb(0.7490196, 0.2470588, 0.7490196))
p <- p +  geom_line(aes(x = t, y = fit5), size = 1.4, linetype = "longdash", col =   rgb(1, 0, 1))
p <- p + theme_bw(base_size = 40) 
p



# Heat maps for the number cases per million

Y = covid.c.r[, 2:424]
Y[Y<0] <- 0 

worldMap <- getMap()
covid.c.r[4,1 ] <- "Bosnia and Herz."
covid.c.r[8,1 ] <- "Czech Rep."
europeanUnion = covid.c.r[, 1]
indEU <- which(worldMap$NAME%in%europeanUnion)
europeCoords <- lapply(indEU, function(i){
  df <- data.frame(worldMap@polygons[[i]]@Polygons[[1]]@coords)
  df$region =as.character(worldMap$NAME[i])
  colnames(df) <- list("long", "lat", "region")
  return(df)
})

library(rnaturalearth)
world_map <- ne_countries(scale = 50, returnclass = 'sf')
europe_map <- world_map %>% filter(name %in% europeanUnion)
cut.offs1 = c(fit.q1c$mu[67], fit.q2c$mu[67], fit.q3c$mu[67], fit.q4c$mu[67], 
              fit.q5c$mu[67]) #30th March 67
cut.offs2 = c( fit.q1c$mu[189], fit.q2c$mu[189], fit.q3c$mu[189], fit.q4c$mu[189], 
               fit.q5c$mu[189] ) # 30th July 189
cut.offs3 = c(fit.q1c$mu[312], fit.q2c$mu[312], fit.q3c$mu[312], fit.q4c$mu[312], 
              fit.q5c$mu[312]) # 30th Nov 312
cut.offs4 = c(fit.q1c$mu[423], fit.q2c$mu[423], fit.q3c$mu[423], fit.q4c$mu[423], 
              fit.q5c$mu[423]) # 21st March 423
covid.c.r[, 1]
value = findInterval(Y[, 67], cut.offs1) + 1
value
europe_map$sovereignt

europe_map$value = c(value[1], value[2], value[3], value[5], value[4], value[33], value[7], value[8], value[13], value[9], value[31], value[10], value[11], value[12],
                     value[34], value[14], value[6], value[15], value[17], value[16], value[18], value[20], value[21],
                     value[19], value[22], value[23], value[24], value[25], value[26], value[27], value[28], value[29], value[30],
                     value[32])

P1 <- ggplot() + geom_sf(data = europe_map, aes(fill = value)) + theme_void() #+ geom_sf(data = circle, fill = 'black', alpha = .2) + 
P1 <- P1 + coord_sf(xlim = c(-08, 28),  ylim = c(32, 62))
P1 <- P1  + scale_fill_gradient(name = "Growth Rate", low = "lightyellow", high = "darkred", na.value = "grey50",
                                guide = FALSE)

P1 <- P1 + theme(#panel.grid.minor = element_line(colour = NA), panel.grid.minor = element_line(colour = NA),
  panel.background = element_rect(colour = NA, fill = "white"), 
  axis.text.x = element_blank(),
  axis.text.y = element_blank(), axis.ticks.x = element_blank(),
  axis.ticks.y = element_blank(), axis.title = element_blank(),
  rect = element_blank(),
  plot.margin = unit( c(-2, -1, -2, -1), "cm" ) )
P1



europeCoords <- do.call("rbind", europeCoords)

cut.offs1 = c(fit.q1c$mu[67], fit.q2c$mu[67], fit.q3c$mu[67], fit.q4c$mu[67], 
              fit.q5c$mu[67]) #30th March 67
cut.offs2 = c( fit.q1c$mu[189], fit.q2c$mu[189], fit.q3c$mu[189], fit.q4c$mu[189], 
               fit.q5c$mu[189] ) # 30th July 189
cut.offs3 = c(fit.q1c$mu[312], fit.q2c$mu[312], fit.q3c$mu[312], fit.q4c$mu[312], 
              fit.q5c$mu[312]) # 30th Nov 312
cut.offs4 = c(fit.q1c$mu[423], fit.q2c$mu[423], fit.q3c$mu[423], fit.q4c$mu[423], 
              fit.q5c$mu[423]) # 21st March 423
# First panel
value = findInterval(Y[, 67], cut.offs1) + 1

europeanUnionTable <- data.frame(country = europeanUnion, value = value)
europeCoords$value <- europeanUnionTable$value[match(europeCoords$region,europeanUnionTable$country)]

P1 <- ggplot() + geom_polygon(data = europeCoords, aes(x = long, y = lat, group = region, fill = value),
                              colour = "black", size = 0.1) +
  coord_map(xlim = c(-08, 28),  ylim = c(32, 62))

P1 <- P1 + scale_fill_gradient(name = "Growth Rate", low = "lightyellow", high = "darkred", na.value = "grey50",
                               guide = FALSE)

P1 <- P1 + theme(#panel.grid.minor = element_line(colour = NA), panel.grid.minor = element_line(colour = NA),
  panel.background = element_rect(colour = NA, fill = "white"), 
  axis.text.x = element_blank(),
  axis.text.y = element_blank(), axis.ticks.x = element_blank(),
  axis.ticks.y = element_blank(), axis.title = element_blank(),
  rect = element_blank(),
  plot.margin = unit( c(-2, -1, -2, -1), "cm" ) )
P1

# Second panel
value = findInterval(Y[, 189], cut.offs2) + 1

value
europe_map$sovereignt

europe_map$value = c(value[1], value[2], value[3], value[5], value[4], value[33], value[7], value[8], value[13], value[9], value[31], value[10], value[11], value[12],
                     value[34], value[14], value[6], value[15], value[17], value[16], value[18], value[20], value[21],
                     value[19], value[22], value[23], value[24], value[25], value[26], value[27], value[28], value[29], value[30],
                     value[32])

P2 <- ggplot() + geom_sf(data = europe_map, aes(fill = value)) + theme_void() #+ geom_sf(data = circle, fill = 'black', alpha = .2) + 
P2 <- P2 + coord_sf(xlim = c(-08, 28),  ylim = c(32, 62))
P2 <- P2  + scale_fill_gradient(name = "Growth Rate", low = "lightyellow", high = "darkred", na.value = "grey50",
                                guide = FALSE)

P2 <- P2 + theme(#panel.grid.minor = element_line(colour = NA), panel.grid.minor = element_line(colour = NA),
  panel.background = element_rect(colour = NA, fill = "white"), 
  axis.text.x = element_blank(),
  axis.text.y = element_blank(), axis.ticks.x = element_blank(),
  axis.ticks.y = element_blank(), axis.title = element_blank(),
  rect = element_blank(),
  plot.margin = unit( c(-2, -1, -2, -1), "cm" ) )
P2



europeanUnionTable <- data.frame(country = europeanUnion, value = value)
europeCoords$value <- europeanUnionTable$value[match(europeCoords$region,europeanUnionTable$country)]

P2 <- ggplot() + geom_polygon(data = europeCoords, aes(x = long, y = lat, group = region, fill = value),
                              colour = "black", size = 0.1) +
  coord_map(xlim = c(-08, 28),  ylim = c(32, 62))

P2 <- P2 + scale_fill_gradient(name = "Growth Rate", low = "lightyellow", high = "darkred", na.value = "grey50",
                               guide = FALSE)

P2 <- P2 + theme(#panel.grid.minor = element_line(colour = NA), panel.grid.minor = element_line(colour = NA),
  panel.background = element_rect(colour = NA, fill = "white"), 
  axis.text.x = element_blank(),
  axis.text.y = element_blank(), axis.ticks.x = element_blank(),
  axis.ticks.y = element_blank(), axis.title = element_blank(),
  rect = element_blank(),
  plot.margin = unit( c(-2, -1, -2, -1), "cm" ) )
P2

# Third panel 

value = findInterval(Y[, 312], cut.offs3) + 1
value
europe_map$sovereignt

europe_map$value = c(value[1], value[2], value[3], value[5], value[4], value[33], value[7], value[8], value[13], value[9], value[31], value[10], value[11], value[12],
                     value[34], value[14], value[6], value[15], value[17], value[16], value[18], value[20], value[21],
                     value[19], value[22], value[23], value[24], value[25], value[26], value[27], value[28], value[29], value[30],
                     value[32])

P3 <- ggplot() + geom_sf(data = europe_map, aes(fill = value)) + theme_void() #+ geom_sf(data = circle, fill = 'black', alpha = .2) + 
P3 <- P3 + coord_sf(xlim = c(-08, 28),  ylim = c(32, 62))
P3 <- P3  + scale_fill_gradient(name = "Growth Rate", low = "lightyellow", high = "darkred", na.value = "grey50",
                                guide = FALSE)

P3 <- P3 + theme(#panel.grid.minor = element_line(colour = NA), panel.grid.minor = element_line(colour = NA),
  panel.background = element_rect(colour = NA, fill = "white"), 
  axis.text.x = element_blank(),
  axis.text.y = element_blank(), axis.ticks.x = element_blank(),
  axis.ticks.y = element_blank(), axis.title = element_blank(),
  rect = element_blank(),
  plot.margin = unit( c(-2, -1, -2, -1), "cm" ) )
P3



europeanUnionTable <- data.frame(country = europeanUnion, value = value)
europeCoords$value <- europeanUnionTable$value[match(europeCoords$region,europeanUnionTable$country)]

P3 <- ggplot() + geom_polygon(data = europeCoords, aes(x = long, y = lat, group = region, fill = value),
                              colour = "black", size = 0.1) +
  coord_map(xlim = c(-08, 28),  ylim = c(32, 62))

P3 <- P3 + scale_fill_gradient(name = "Growth Rate", low = "lightyellow", high = "darkred", na.value = "grey50",
                               guide = FALSE)

P3 <- P3 + theme(#panel.grid.minor = element_line(colour = NA), panel.grid.minor = element_line(colour = NA),
  panel.background = element_rect(colour = NA, fill = "white"), 
  axis.text.x = element_blank(),
  axis.text.y = element_blank(), axis.ticks.x = element_blank(),
  axis.ticks.y = element_blank(), axis.title = element_blank(),
  rect = element_blank(),
  plot.margin = unit( c(-2, -1, -2, -1), "cm" ) )
P3

# Fourth panel 

value = findInterval(Y[, 423], cut.offs4) + 1

value
europe_map$sovereignt

europe_map$value = c(value[1], value[2], value[3], value[5], value[4], value[33], value[7], value[8], value[13], value[9], value[31], value[10], value[11], value[12],
                     value[34], value[14], value[6], value[15], value[17], value[16], value[18], value[20], value[21],
                     value[19], value[22], value[23], value[24], value[25], value[26], value[27], value[28], value[29], value[30],
                     value[32])

P4 <- ggplot() + geom_sf(data = europe_map, aes(fill = value)) + theme_void() #+ geom_sf(data = circle, fill = 'black', alpha = .2) + 
P4 <- P4 + coord_sf(xlim = c(-08, 28),  ylim = c(32, 62))
P4 <- P4  + scale_fill_gradient(name = "Growth Rate", low = "lightyellow", high = "darkred", na.value = "grey50",
                                guide = FALSE)

P4 <- P4 + theme(#panel.grid.minor = element_line(colour = NA), panel.grid.minor = element_line(colour = NA),
  panel.background = element_rect(colour = NA, fill = "white"), 
  axis.text.x = element_blank(),
  axis.text.y = element_blank(), axis.ticks.x = element_blank(),
  axis.ticks.y = element_blank(), axis.title = element_blank(),
  rect = element_blank(),
  plot.margin = unit( c(-2, -1, -2, -1), "cm" ) )
P4



europeanUnionTable <- data.frame(country = europeanUnion, value = value)
europeCoords$value <- europeanUnionTable$value[match(europeCoords$region,europeanUnionTable$country)]

P4 <- ggplot() + geom_polygon(data = europeCoords, aes(x = long, y = lat, group = region, fill = value),
                              colour = "black", size = 0.1) +
  coord_map(xlim = c(-08, 28),  ylim = c(32, 62))

P4 <- P4 + scale_fill_gradient(name = "Growth Rate", low = "lightyellow", high = "darkred", na.value = "grey50",
                               guide = FALSE)

P4 <- P4 + theme(#panel.grid.minor = element_line(colour = NA), panel.grid.minor = element_line(colour = NA),
  panel.background = element_rect(colour = NA, fill = "white"), 
  axis.text.x = element_blank(),
  axis.text.y = element_blank(), axis.ticks.x = element_blank(),
  axis.ticks.y = element_blank(), axis.title = element_blank(),
  rect = element_blank(),
  plot.margin = unit( c(-2, -1, -2, -1), "cm" ) )
P4
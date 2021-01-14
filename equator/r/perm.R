library(foreign)
library(dplyr)
library(magrittr)
library(readstata13)
library(XML)
library(ggplot2)
library(lfe)
library(geosphere)
library(rgdal)
library(viridis)

# PREFERRED COLORS
watercols5 <- viridis.map %>% filter(opt == "C") %>% `[`(,1:3) %>% sapply(function(x) { x[c(1, 64, 128, 192, 256)] }) %>%
  apply(MARGIN = 1, function(x) { paste0("#", round(x*255) %>% as.hexmode %>% paste(collapse = "")) })

# # EXTRACT COUNTRY DATA FROM WORLD BANK
# countries <- readOGR("data/sp/", "ne_50m_admin_0_countries")
# gdp <- xmlToDataFrame(URLencode("http://api.worldbank.org/v2/countries/all/indicator/NY.GDP.MKTP.CD?per_page=20000"))
# pop <- xmlToDataFrame(URLencode("http://api.worldbank.org/v2/countries/all/indicator/SP.POP.TOTL?per_page=20000"))
# gdppc <- gdp %>% filter(date == 2015) %>% rename(gdp = value, isocode = countryiso3code) %>%
#   mutate(isocode = as.character(isocode)) %>% filter(nchar(isocode) == 3) %>%
#   select(isocode, gdp) %>%
#   left_join(pop %>% filter(date == 2015) %>% rename(pop = value, isocode = countryiso3code) %>%
#               mutate(isocode = as.character(isocode)) %>% filter(nchar(isocode) == 3) %>%
#               select(isocode, pop)) %>%
#   mutate(gdp = as.numeric(as.character(gdp)), pop = as.numeric(as.character(pop))) %>%
#   mutate(lgdppc = log(gdp / pop))
# countriescentroid <- countries %>% rgeos::gCentroid(byid = T) %>% slot("coords")
# countries$longitude <- countriescentroid[,"x"]
# countries$latitude <- countriescentroid[,"y"]
# gdppc %<>%
#   left_join(countries@data %>% rename(isocode = ISO_A3) %>% select(isocode, longitude, latitude)) %>%
#   na.omit()
# saveRDS(gdppc, file = "../data/gdppc.rds")

# READ WORLD BANK DATA
data <- readRDS("../data/gdppc.rds")

# READ COUNTRY SHAPEFILE
countries <- readOGR("../data/sp/ne_50m_admin_0_countries.shp")
# GENERATE DATA FRAME FOR PLOTTING COUNTRIES
gcountries <- countries %>% fortify %>%
  left_join(countries@data %>% mutate(id = row.names(countries)) %>%
              rename(isocode = ISO_A3) %>% select(id, isocode)) %>%
  rename(longitude = long, latitude = lat)

# MAIN REGRESSIONS
reg <- lm(lgdppc ~ I(abs(latitude) / 100), data = data)
reg2 <- felm(lgdppc ~ I(abs(latitude) / 100) | 0 | 0 | isocode, data = data)

# GRAB COLUMNS OF DATA FOR SIMULATION
datasim <- data[,c("lgdppc", "latitude", "longitude")]
# SIMULATION WITH RANDOMLY GENERATED EQUATOR
sim.reg <- function(keepdata = F) {
  datasimi <- datasim
  # CONVERT LAT AND LONG TO COORDS ON SPHERE
  datasimi$x <- cos(datasimi$longitude * (pi/180)) * cos(datasimi$latitude * (pi/180))
  datasimi$y <- sin(datasimi$longitude * (pi/180)) * cos(datasimi$latitude * (pi/180))
  datasimi$z <- sin(datasimi$latitude * (pi/180))
  # RANDOMLY DRAW A NEW NORTH POLE
  v <- rnorm(3)
  v <- v / sqrt(sum(v^2))
  # GET ANGLE BETWEEN v AND OLD COORDINATES
  anglei <- acos(v[1] * datasimi$x + v[2] * datasimi$y + v[3] * datasimi$z)
  # CONVERT ANGLE TO LATITUDE
  datasimi$newlatitude <- 90 - (180/pi) * anglei
  # REGRESSION
  regsim <- lm(lgdppc ~ I(abs(newlatitude) / 100), data = datasimi)
  if(!keepdata) {
    return(coef(regsim)[2])
  } else {
    return(list("newpole" = v, "data" = datasimi, "intercept" = coef(regsim)[1], "coef" = coef(regsim)[2]))
  }
}
# SIMULATION WITH NAIVELY PERMUTED LATITUDE
sim.reg2 <- function() {
  datasimi <- datasim
  # RANDOMLY ASSIGN LATITUDE
  datasimi$latitude <- sample(datasimi$latitude)
  # REGRESSION
  regsim <- lm(lgdppc ~ I(abs(latitude) / 100), data = datasimi)
  return(coef(regsim)[2])
}

# RUN SIMULATIONS
set.seed(20190501)
sims <- replicate(100000, sim.reg())
sims2 <- replicate(100000, sim.reg2())
simsplot <- replicate(100, sim.reg(keepdata = T))
# PRINT SPATIAL RANDOMIZATION INFERENCE p value
mean(abs(coef(reg)[2]) < abs(sims))
# PRINT (NAIVE) RANDOMIZATION INFERENCE p value
mean(abs(coef(reg)[2]) < abs(sims2))

# NAIVE VS. SPATIAL
simsdf <- data.frame(name = c(rep("spatial", length(sims)), rep("naive", length(sims2))), "sim" = c(sims, sims2))
gsims <- ggplot(simsdf, aes(x = sim, group = name, linetype = name)) +
  geom_vline(xintercept = coef(reg)[2], col = watercols5[3]) +
  annotate(geom = "label", x = coef(reg)[2], y = 0.4, label = "OLS Estimate") +
  stat_density(data = simsdf %>% subset(name == "spatial"), geom = "line") +
  stat_density(data = simsdf %>% subset(name == "naive"), geom = "line") +
  scale_x_continuous(name = "% change in GDP per degree of absolute latitiude", limits = c(-6, 6)) +
  scale_y_continuous(name = "Density") +
  scale_linetype_discrete(name = "Randomization inference type", labels = c("Naive", "Spatial")) +
  theme(legend.position = "top")
ggsave(filename = "output/permdist.jpg", gsims, width = 6, height = 6)

# CALCULATE DISTANCE MATRIX
datadist <- distm(matrix(c(data$longitude, data$latitude), ncol = 2, byrow = F),
                  matrix(c(data$longitude, data$latitude), ncol = 2, byrow = F),
                  fun = distHaversine)
# CALCULATE RESIDUALS (FOR CONLEY ERRORS)
dataeps <- reg$residuals
datax <- abs(data$latitude) / 100
datax <- datax - mean(datax)

# CONLEY SE
spatialse <- lapply(c(seq(0, 6000, 100)), function(r) {
  r <- r * 1000
  # ROBUST VCOV FORMULA
  ser <- solve(t(datax) %*% datax) %*%
    (t(datax) %*% ((dataeps %*% t(dataeps)) * (datadist <= r)) %*% datax) %*%
    solve(t(datax) %*% datax)
  # SQRT TO GET SE
  ser <- sqrt(ser)
  # EFFECTIVE NUMBER OF CLUSTERS
  nclur <- nrow(data)^2 / sum(datadist <= r)
  # DOF ADJUSTMENT
  sestata <- ser * sqrt((nrow(data) - 1) / (nrow(data) - 2)) * sqrt(nclur / (nclur - 1))
  data.frame(radius = r / 1000, se = ser, sestata = sestata, nclu = nclur)
}) %>% bind_rows() %>%
  # CALCULATE PVALUES
  mutate(pval = 2 * (1 - pnorm(coef(reg)[2] / se)),
         pvalstata = 2 * (1 - pnorm(coef(reg)[2] / sestata)))
# GET P VALUE FROM PERMUTATION INFERENCE
permse <- data.frame(pval = mean(abs(coef(reg)[2]) < abs(sims))) %>%
  # AND "EFFECTIVE" STANDARD ERROR
  mutate(se = coef(reg)[2] / qnorm(1 - pval / 2))

# PLOT PERFORMANCE OF CONLEY STANDARD ERRORS
gsenclu <- ggplot(data = data.frame()) +
  # geom_line(data = spatialse, aes(x = radius, y = se), col = watercols5[3]) +
  geom_line(data = spatialse, aes(x = radius, y = sestata), col = watercols5[3]) +
  annotate(geom = "label", x = 3000, y = spatialse$sestata[which(spatialse$radius == 3000)]-0.2,
           label = "Conley standard error") +
  geom_line(data = spatialse %>% mutate(nclu = (nclu*2/200)+0.4), aes(x = radius, y = nclu)) +
  annotate(geom = "label", x = 3000, y = (spatialse$nclu[which(spatialse$radius == 3000)]*2/200)+0.4+0.2,
           label = "Effective number of clusters") +
  annotate(geom = "segment", x = -Inf, xend = 6000, y = permse$se, yend = permse$se, linetype = 2,
           col = watercols5[3]) +
  annotate(geom = "label", x = 3000, y = permse$se, label = "Randomization inference\neffective standard error") +
  scale_x_continuous(name = "Radius (kilometers)", limits = c(0, 6000)) +
  scale_y_continuous(name = "Conley standard error", limits = c(0.4, 2.4),
                     sec.axis = sec_axis(~(.-0.4)*200/2, name = "Effective number of clusters",
                                         breaks = seq(0, 200, 40)),
                     breaks = seq(0.4, 2.4, 0.4))
ggsave(filename = "output/conleyse.jpg", gsenclu, width = 6, height = 6)

ggsave(filename = "output/sefig.jpg",
       plot = gridExtra::grid.arrange(gsims + ggtitle("Randomization inference") +
                                        theme(plot.title = element_text(face = "bold")),
                                      gsenclu + ggtitle("Comparison with Conley standard errors") +
                                        theme(plot.title = element_text(face = "bold")), ncol = 1),
       width = 6, height = 6)

# VISUALIZING SIMS
getgreatcircle <- function(v) {
  # CREATE A UNIT VECTOR ORTHOGONAL TO NEW POLE
  M2 <- c(v[2], -v[1], 0)
  M2 <- M2 / sqrt(sum(M2^2))
  # GET CROSS PRODUCT OF POLE AND PREVIOUS VECTOR
  M1 <- pracma::cross(v, M2)
  # CREATE ROTATION MATRIX THAT ROTATES (0, 0, 1) TO NEW POLE
  M <- matrix(c(M1, M2, v),
              byrow = T, nrow = 3)
  # CREATE EQUATOR
  circ <- matrix(c(cos(seq(0, 2*pi, 0.001*pi)),
                   sin(seq(0, 2*pi, 0.001*pi)),
                   rep(0, 2001)),
                 byrow = F, ncol = 3)
  # ROTATE TO NEW EQUATOR
  circ <- circ %*% M
  list("gc" = data.frame(x = circ[,1], y = circ[,2], z = circ[,3]) %>%
         mutate(latitude = asin(z)*(180/pi), longitude = atan2(y,x)*(180/pi)),
       "newpole" = data.frame(x = v[1], y = v[2], z = v[3]) %>%
         mutate(latitude = asin(z)*(180/pi), longitude = atan2(y,x)*(180/pi)))
}
for(i in 1:20) {
  # GET COEFFICIENTS FROM FIRST i SIMS
  dfcoefs <- lapply(1:i, function(j) {
    data.frame("intercept" = simsplot[,j]$intercept, "coef" = simsplot[,j]$coef)
  }) %>% bind_rows() %>% mutate(coef = coef / 100)
  # REGRESSIONS FIGURE
  greg <- ggplot(data = data.frame())
  if(i > 1) {
    # ADD TRANSPARENT LINES FOR ALL PREVIOUS SIMS
    greg <- greg + geom_abline(data = dfcoefs[1:(i-1),], aes(intercept = intercept, slope = coef), alpha = 0.3, col = watercols5[3])
  }
  greg <- greg +
    geom_abline(data = dfcoefs[i,], aes(intercept = intercept, slope = coef), col = watercols5[3]) +
    geom_point(data = data %>% mutate(abslat = abs(latitude)), aes(x = abslat, y = lgdppc),
               shape = 1) +
    geom_point(data = simsplot[,i]$data %>% mutate(abslat = abs(newlatitude)), aes(x = abslat, y = lgdppc), col = watercols5[3],
               shape = 1) +
    geom_abline(slope = coef(reg)[2] / 100, intercept = coef(reg)[1]) +
    scale_x_continuous(limits = c(0, 90), name = "Absolute latitude", breaks = seq(0, 80, 20),
                       labels = paste0(seq(0, 80, 20), "°")) +
    scale_y_continuous(name = "log GDP per capita", limits = range(data$lgdppc))
  # MAP FIGURE
  gmap <- ggplot(data = data.frame(), aes(x = longitude, y = latitude)) +
    geom_polygon(data = gcountries %>% left_join(data %>% select(isocode, lgdppc)),
                 aes(fill = lgdppc, group = group), col = "white") +
    annotate(geom = "segment", x = -180, xend = 180, y = 0, yend = 0, col = "black") +
    geom_path(data = getgreatcircle(simsplot[,i]$newpole)$gc %>% arrange(longitude), col = watercols5[3], linetype = 2) +
    # geom_point(data = getgreatcircle(simsplot[,7]$newpole)$newpole, col = "black") +
    # scale_fill_viridis_c(option = "C", na.value = "white") +
    scale_fill_gradient(na.value = "white", low = "grey70", high = "grey10") +
    theme(legend.position = "none") +
    scale_x_continuous(labels = c("100°W", "0°E", "100°E"), breaks = c(-100, 0, 100), limits = c(-180, 180), name = NULL) +
    scale_y_continuous(labels = c("50°S", "0°N", "50°N"), breaks = c(-50, 0, 50), limits = c(-90, 90), name = NULL)
  ggsave(filename = paste0("output/perm", sprintf("%02d", i), ".jpg"),
         plot = gridExtra::grid.arrange(gmap, greg, ncol = 1),
         width = 6, height = 6)
}
# GENERATE GIF FOR SIMS
system("convert -delay 100 -loop 0 output/perm*.jpg output/permgif.gif")

## Create a polygon from the buffer and intesect with a grid to extract wind data
##
## Data has been extracted and figures plotted using 'extractWINDlevel4_csv.R' in 'SatelliteTS'

library(raster)
library(rasterVis)
library(xts)
library(maptools)
library(gridExtra)
library(grid)
library(ncdf4)
library(data.table)
library(rgdal)


setwd("~/R/myFunctions/")

# setup bathymetry map
bathy <- raster("/home/robert/R/x86_64-pc-linux-gnu-library/3.4/bathy-GEBCO/gebco-southernAfrica_1min.nc")
locbb <- matrix(c(10,-40,40,-20),nrow=2,ncol=2)
x <- extent(locbb)
bathy.cont <- crop(bathy,x)
bathy.loc <- crop(bathy,x)      # crop bathy map
bathy.land <- crop(bathy,x)
crs <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
rm(bathy,locbb,x)

# get SA borders
borders <- readShapePoly("SAfrica-borders.shp", proj4string = crs)
borders <- SpatialPolygons(borders@polygons)
minlat <- -40;minlon <- 10;maxlat <- -20;maxlon <- 40
mapThemeBathy <- colorRampPalette(c("blue","white"),alpha=T)

# get station locations
station.loc <- read.table("~/R_projects-CS/ame-temporalbabe/CS-ExtractFeatures/station_locations.Rdata")
station.pts <- SpatialPointsDataFrame(station.loc[,1:2],station.loc,proj4string = crs)

setwd("~/R_projects-CS/ame-temporalbabe/ExtractWindFields/")

# ------------------------- Surface Wind
baseURL <- "/media/robert/Seagate Backup Plus Drive"
type <- "SurfaceWind"
product <- "CCMP"

listFiles <- list.files(paste(baseURL,type,product,"csv/",sep = "/"),full.names = TRUE)
pts <- fread(listFiles[1],header = TRUE, sep = ",")
r.raster <- rasterFromXYZ(cbind(pts$lon,pts$lat,pts$wndSpd),crs=crs)

# if bathy.mask does not exist ---------------------------------------------------
if (!exists("bathy.mask")){
  # extract the coastal quantiles from the raster layer
  buf.dist <- 300000                                  # coastal boudary distance metres
  bathy.land[bathy.land >= 0] = NA                    # set land values to NA
  bathy.loc[bathy.loc > 5] <- NA                      # select 5 m height as coastal limits
  bathy.loc[bathy.loc < -5] <- NA
  bathy.buf <- buffer(bathy.loc, width=buf.dist, doEdge=TRUE)
  bathy.mask.temp <- mask(bathy.buf,bathy.land)       # remove the inland buffer distance
  
  r.raster@crs <- crs
  bathy.mask <- projectRaster(bathy.mask.temp,r.raster)  # resample bathymetry to data resolution
  
  at <- c(6000,0,-200,-500,-1000,-2000,-3000,-4000,-5000,-6000,-7000)
  bathy.mask.plot <- bathy.mask
  bathy.mask.plot[is.na(bathy.mask.plot)] <- 0
  levelplot(bathy.mask.plot, margin=F, col.regions=c("white","yellow"), colorkey=F, 
            xlab="Longitude",ylab="Latitude", contour=T,
            main=paste(product," 300 km Coastal Boundary Mask")) +
    contourplot(bathy.cont,region=T,col.regions=mapThemeBathy, alpha.regions=.3,
                at=at,margin = F, labels = F) +
    layer(sp.lines(borders, col = "black", lwd = 1)) + 
    layer(sp.points(station.pts, col="red", pch=1, cex=1.5, fill=F, lwd=2))
  
  # save the mask image and mask raster layer
  filename1 <- paste0("Output/",product,"_300kmMask.grd")
  writeRaster(bathy.mask,filename1,"raster", overwrite=T)
  
  } else {
    # open existing raster data
    open.file1 <- paste0("Output/",product,"_300kmMask.grd")
    bathy.mask <- raster(open.file1)
    rm(open.file1)
}

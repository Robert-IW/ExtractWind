## To extract Satellite data along a give virtual transect
## Requires loading of saved transects

## Creates a list 'raster.list' for each station and each day containing a list of 4
##  rasterbricks, 1 for each time interval, with 4 layers; 'Uwnd','Vwnd','wndDir'.'wndSpd'

## Creates a list 'wind.along/offshore' similarly with each of the 4 list items containing a
##  matrix with 'cell','Uwnd','Vwnd','wndDir','wndSpd' extracted with the 'trans.along/offshore'
##  transects

## Creates the list 'wind.along.dat/off.dat' similar to 'wind.along/offshore but with the
##  addistion of 'x','y' coords


require(plotrix)
require(ncdf4)
require(rgeos)
require(raster)
require(rasterVis)
require(maptools)
require(ncdf4)
require(grid)

doGraphic = F

# required functions
setwd("~/R/myFunctions/")
source("func_windFromCCMP.R")

# create a border map for plotting
proj <- CRS('+proj=longlat +datum=WGS84')

bathy <- raster("/home/robert/R/x86_64-pc-linux-gnu-library/3.3/bathy-GEBCO/gebco-southernAfrica_1min.nc")
mapThemeBathy <- colorRampPalette(c("blue","white","brown"),alpha=F)
at <- c(5000,2000,1000,500,200,100,50,0,-50,-100,-200,-500,-1000,-2000,-5000,-10000)
my.at <- seq(0,18,length.out = 10)

station.loc <- read.table("~/R_projects-CS/ame-temporalbabe/CS-ExtractFeatures/station_locations.Rdata",stringsAsFactors=F)
station.pts <- SpatialPointsDataFrame(station.loc[,1:2],station.loc,proj4string = proj)

# For each station
for (h in 93:nrow(station.loc)){
  
  stat.name <- station.loc$station[h]
  
  file.load <- paste0("/home/robert/R/myFunctions/Wind_transects/",stat.name,"_transects.Rdata")
  
  load(file.load)
  rm(file.load)
  
  trans.along <- station.list$alongshore
  trans.offshore <- station.list$offshore
  region.wind <- station.list$local
  rm(station.list)
  
  a <- station.pts@data[h,1:2]
  coordinates(a) <- ~lon+lat
  
  # load all the transects
  
  # WIND ----------------------------------------------------------------------
  #sourceURL <- ("/media/robert/KELP-HDD-Portable/SurfaceWind/CCMP/data")
  sourceURL <- ("/media/robert/Seagate Backup Plus Drive/Backup Satellite Data/Wind/CCMP/data/")
  
  #mapTheme <- rasterTheme(region=brewer.pal(8,"Blues"))
  my.at <- seq(0,18,length.out = 10)
  
  # get the list of year directories
  dirList <- list.dirs(sourceURL,recursive = FALSE)
  
  # get the list of day files in each annual directory
  for (i in 1:length(dirList)){
    
    fileList <- list.files(dirList[i])
    year.name <- basename(dirList[i])      # for the saved files
    
    # create filename and check if it exists
    filename1 <- paste0("/media/robert/KELP-HDD-Portable/SurfaceWind/CCMP/extractCSV/",stat.name,"_",
                        year.name,"_alongshore_wind.csv")
    filename2 <- paste0("/media/robert/KELP-HDD-Portable/SurfaceWind/CCMP/extractCSV/",stat.name,"_",
                        year.name,"_offshore_wind.csv")
    filename3 <- paste0("/media/robert/KELP-HDD-Portable/SurfaceWind/CCMP/extractCSV/",stat.name,"_",
                        year.name,"_local_wind.csv")
    
    if (file.exists(filename1) & file.exists(filename2) & file.exists(filename3)){
      print("Files exist")
      next
    }
    
    # for each day get the 4 layers @ 6 hourly intervals
    
    # create an empty data frame
    data.store.along <- data.frame("date"=character(),"lon"=numeric(),"lat"=numeric(),"time"=character(),
                                   "Uwnd"=numeric(),"Vwnd"=numeric(),"wndDir"=numeric(),"wndSpd"=numeric(),
                                   stringsAsFactors = F)
    
    data.store.off <- data.frame("date"=character(),"lon"=numeric(),"lat"=numeric(),"time"=character(),
                                 "Uwnd"=numeric(),"Vwnd"=numeric(),"wndDir"=numeric(),"wndSpd"=numeric(),
                                 stringsAsFactors = F)
    
    data.store.local <- data.frame("date"=character(),"lon"=numeric(),"lat"=numeric(),"time"=character(),
                                   "Uwnd"=numeric(),"Vwnd"=numeric(),"wndDir"=numeric(),"wndSpd"=numeric(),
                                   stringsAsFactors = F)
    
    for (j in 1:length(fileList)){
      
      print(paste("Working on CCMP file ",j," in ",length(fileList),sep = ''))
      
      fileURL <- paste(dirList[i],"/",fileList[j],sep = '')
      
      # get list of raster bricks for each hour interval (brick = 'Uwnd','Vwnd','wndDir','wndSpd')
      raster.list <- windFromCCMP(fileURL)
      
      # ------------- along shore transect
      # extract U wind, V wind, speed and direction along coastal transect
      wind.along <- lapply(raster.list, function(r) extract(r, trans.along, cellnumbers=TRUE))      # a list of matrices
      
      # get the coordinates of the extracted points
      wind.along.loc <- xyFromCell(raster.list[[1]], wind.along[[1]][[1]][,'cell'])   # a matrix
      colnames(wind.along.loc) <- c("lon","lat")
      
      # combine the lon lat and data
      wind.along.dat <- lapply(wind.along, function(r) cbind(wind.along.loc, r[[1]]))
      rm(wind.along, wind.along.loc)
      
      temp <- lapply(seq_along(wind.along.dat), function(r) data.frame("date"=rep(strftime(names(wind.along.dat[r]),format="%Y%m%d"),nrow(wind.along.dat[[r]])),
                                                                       "time"=rep(strftime(names(wind.along.dat[r]),format="%H%M"),nrow(wind.along.dat[[r]])),
                                                                       wind.along.dat[[r]],
                                                                       stringsAsFactors = F))
      wind.along.save <- do.call(rbind, temp)
      wind.along.save <- wind.along.save[,-5]
      rm(temp)
      
      # ------------- offshore transect
      # extract U wind, V wind, speed and direction along offshore transect
      wind.offshore <- lapply(raster.list, function(r) extract(r,trans.offshore,cellnumbers=TRUE))      # a list of matrices
      
      # get the coordinates of the extracted points
      wind.off.loc <- xyFromCell(raster.list[[1]], wind.offshore[[1]][[1]][,'cell'])   # a matrix
      colnames(wind.off.loc) <- c("lon","lat")
      
      # combine the lon lat and data
      wind.off.dat <- lapply(wind.offshore, function(r) cbind(wind.off.loc, r[[1]]))
      rm(wind.offshore, wind.off.loc)
      
      temp <- lapply(seq_along(wind.off.dat), function(r) data.frame("date"=rep(strftime(names(wind.off.dat[r]),format="%Y%m%d"),nrow(wind.off.dat[[r]])),
                                                                     "time"=rep(strftime(names(wind.off.dat[r]),format="%H%M"),nrow(wind.off.dat[[r]])),
                                                                     wind.off.dat[[r]], stringsAsFactors = F))
      wind.off.save <- do.call(rbind, temp)
      wind.off.save <- wind.off.save[,-5]   # remove "cell" column
      rm(temp)
      
      # -------------- local wind region
      # extract U wind, V wind, speed and direction along offshore transect
      wind.local.med <- lapply(raster.list, function(r) extract(r,region.wind, fun=median))      # a list of matrices
      wind.local <- lapply(raster.list, function(r) extract(r, region.wind, cellnumbers=TRUE))
      
      # get the coordinates of the centre of the polygon
      wind.local.loc.med <- gCentroid(region.wind)@coords
      colnames(wind.local.loc.med) <- c("lon","lat")
      
      # combine the lon lat and data
      wind.local.dat.med <- lapply(wind.local.med, function(r) cbind(wind.local.loc.med,r))
      
      temp <- lapply(seq_along(wind.local.dat.med), function(r) data.frame("date"=rep(strftime(names(wind.local.dat.med[r]),format="%Y%m%d"),nrow(wind.local.dat.med[[r]])),
                                                                           "time"=rep(strftime(names(wind.local.dat.med[r]),format="%H%M"),nrow(wind.local.dat.med[[r]])),
                                                                           wind.local.dat.med[[r]], stringsAsFactors = F))
      wind.local.save <- do.call(rbind, temp)
      rm(temp)
      
      # ----------------------------------------------------------------------------------------------------- PRINT  
      # for visualizing the regional wind use all wind in the polygon
      if (doGraphic){
        
        wind.local.loc <- xyFromCell(raster.list[[1]], wind.local[[1]][[1]][,'cell'])
        colnames(wind.local.loc) <- c("lon","lat")
        
        # combine the lon lat and data and convert to data frame
        wind.local.dat <- lapply(wind.local, function(r) cbind(wind.local.loc, r[[1]]))
      }
      
      # -------------- store data for each year
      # bind the data frames for alongshore, offshore and local extractions
      data.store.along <- rbind(data.store.along, wind.along.save)
      data.store.off <- rbind(data.store.off, wind.off.save)
      data.store.local <- rbind(data.store.local, wind.local.save)
      
      rm(wind.along.save,wind.off.save,wind.local.save)
      
      # -------------- plot
      # for the first day of each year plot and save an png of the extracted wind
      
      if (doGraphic & j==1){
        
        print("Creating and saveing sample png")
        # ---------------------- along coast transect
        # combine the lon lat locations with the extracted values for Speed and Direction as list of data frames
        df.along.plot <- lapply(wind.along.dat, function(r) data.frame("lon"=r[,'lon'],"lat"=r[,'lat'],
                                                                       "wndSpd"=r[,'wndSpd'],"wndDir"=r[,'wndDir'],
                                                                       stringsAsFactors = FALSE))
        
        wind.along.plot  <- rasterFromXYZ(df.along.plot[[1]], crs = proj, res = 0.25)
        
        
        # --------------------- offshore transect
        df.off.plot <- lapply(wind.off.dat, function(r) data.frame("lon"=r[,'lon'],"lat"=r[,'lat'],
                                                                   "wndSpd"=r[,'wndSpd'],"wndDir"=r[,'wndDir'],
                                                                   stringsAsFactors = FALSE))
        
        wind.off.plot  <- rasterFromXYZ(df.off.plot[[1]], crs = proj, res = 0.25)
        
        # --------------------- local area
        df.local.plot <- lapply(wind.local.dat, function(r) data.frame("lon"=r[,'lon'],"lat"=r[,'lat'],
                                                                       "wndSpd"=r[,'wndSpd'],"wndDir"=r[,'wndDir'],
                                                                       stringsAsFactors = FALSE))
        
        wind.local.plot  <- rasterFromXYZ(df.local.plot[[1]], crs = proj, res = 0.25)
        
        # ---------------------- merge and plot transect data
        # increase the extent of the images to the maximum that includes the offshore transect
        extent.tmp1 <- extent(wind.off.plot)
        extent.tmp2 <- extent(wind.along.plot)
        
        ext.xmin <- min(c(extent.tmp1@xmin, extent.tmp2@xmin))      # most westerly
        #ext.xmax <- max(c(extent.tmp1@xmax, extent.tmp2@xmax))      # most easterly
        ext.ymin <- min(c(extent.tmp1@ymin, extent.tmp2@ymin))      # most easterly
        #ext.ymax <- max(c(extent.tmp1@ymax, extent.tmp2@ymax))      # most easterly
        
        #extent.tmp <- extent(c(ext.xmin-0.5, ext.xmax+0.5, ext.ymin-0.5, ext.ymax+0.5))
        extent.tmp <- extent(c(ext.xmin-0.7, ext.xmin+4, ext.ymin-0.7, ext.ymin+4))
        rm(extent.tmp1, extent.tmp2, ext.xmin, ext.xmax, ext.ymin, ext.ymax)
        
        terrain.loc <- crop(bathy, extent.tmp)
        
        # merge the rasters
        along.off.plot <- merge(wind.along.plot, wind.off.plot)
        #along.off.plot <- setExtent(along.off.plot, extent.tmp, keepres=F, snap=F)
        
        mapTheme <- rasterTheme(region=colorRampPalette(c("darkturquoise","cornflowerblue","darkmagenta","red","yellow"))(100))
        
        lattice.options(
          layout.heights=list(bottom.padding=list(x=1), top.padding=list(x=1)),
          layout.widths=list(left.padding=list(x=1), right.padding=list(x=2))
        )
        
        filename1 <- paste0("~/R_projects-CS/ame-temporalbabe/ExtractWindFields/Output/",stat.name,"_",
                            year.name,"_wind_sample.png")
        png(filename1, width = 5, height = 5, units = 'in', res = 300)
        
        main.title <- paste0(stat.name," Extracted Wind Vectors\n",names(raster.list[1]))
        scale.min <- min(along.off.plot[[1]]@data@values,na.rm=T)
        
        p <- vectorplot(along.off.plot, par.settings=mapTheme,region=TRUE,scaleSlope=1,
                        isField=TRUE, at=seq(0,18,0.5),
                        narrows=1000,unit="degrees",lwd=3, col="black", alpha.regions=.99,length=.1,
                        colorkey=list(at=my.at),
                        scales=list(tick.number=4, tck=1),
                        #aspX=.05,aspY=.05,
                        #xlim=c(extent.tmp@xmin,extent.tmp@xmax),
                        #ylim=c(extent.tmp@ymin,extent.tmp@ymax),
                        main=list(label=main.title, cex=1.4),
                        ylab=list("Latitude",cex=1.2), xlab=list("Longitude",cex=1.2)) +
          
          contourplot(terrain.loc, region=T,col.regions=mapThemeBathy, alpha.regions=.3,
                      at=at, margin = F, labels = F) +
          
          layer(sp.polygons(region.wind, col="red",lwd=4)) +
          layer(sp.points(a, col="orangered", pch=19, cex=1.5, fill=F, lwd=2)) +
          layer(sp.points(a, col="orangered", pch=1, cex=3, fill=F, lwd=2)) +
          layer(sp.points(a, col="orangered", pch=1, cex=6, fill=F, lwd=2)) +
          
          vectorplot(wind.local.plot, par.settings=mapTheme,col="greenyellow",region=FALSE,
                     scaleSlope=1,isField=T, narrows=8,unit="degrees",lwd=3, length=.1)
        
        p <- update(p, aspect="iso",xlim=c(extent.tmp@xmin,extent.tmp@xmax),
                    ylim=c(extent.tmp@ymin,extent.tmp@ymax))
        print(p)
        trellis.focus("legend", side="right", clipp.off=TRUE, highlight=FALSE)
        grid.text(expression(paste("Wind Speed [",m.s^-1,"]",sep="")),1.4,0.50,rot=90)
        trellis.unfocus()
        dev.off()
        
        rm(wind.along.dat, wind.off.dat, wind.local.dat, filename1)
        
      } else {
        rm(wind.along.dat, wind.off.dat, wind.local.dat, wind.local.dat.med)
      }# if doGraphic
    } # for each day
    #---------------------------------------------------------------------------------------------------- END PRINT
    
    filename1 <- paste0("/media/robert/KELP-HDD-Portable/SurfaceWind/CCMP/extractCSV/",stat.name,"_",
                        year.name,"_alongshore_wind.csv")
    filename2 <- paste0("/media/robert/KELP-HDD-Portable/SurfaceWind/CCMP/extractCSV/",stat.name,"_",
                        year.name,"_offshore_wind.csv")
    filename3 <- paste0("/media/robert/KELP-HDD-Portable/SurfaceWind/CCMP/extractCSV/",stat.name,"_",
                        year.name,"_local_wind.csv")
    
    print("Writing CSV files")
    write.csv(data.store.along, file = filename1, row.names = FALSE)
    write.csv(data.store.off, file = filename2, row.names = FALSE)
    write.csv(data.store.local, file = filename3, row.names = FALSE)
    
    rm(data.store.along,data.store.off,data.store.local)
    rm(filename1,filename2,filename3)
    gc()
  } # for each year
} # for each station

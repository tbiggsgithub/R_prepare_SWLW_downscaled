# R_prepare_SWLW_downscaled

```r
#  Load, subset and prepare Rn images for SEBAL study area
# For other years, modify the year on line "patt" and "indir.albedo"

#  MCD43A described in http://www-modis.bu.edu/brdf/userguide/param.html#sds
#  24-hour Rn taken from Bisht and Bras 2010, Eq 12
#  Rn = 2Rn.instant / (pi*sin(pi*((top-trise)/(tset-trise))))
#  trise = tsunrise+1hr   tset = tsunset-1hr 

if (exists(as.character(substitute(sw1030stack)))) {remove(sw1030stack)}
if (exists(as.character(substitute(lw1030stack)))) {remove(lw1030stack)}
if (exists(as.character(substitute(sw24stack)))) {remove(sw24stack)}
if (exists(as.character(substitute(lw24stack)))) {remove(lw24stack)}
if (exists(as.character(substitute(Rn1030stack)))) {remove(sw1030stack)}
if (exists(as.character(substitute(Rn1030stack)))) {remove(sw1030stack)}

library(raster)
library(rgeos)
library(maptools)
library(ncdf)
library(rgdal)

start.end = c(2011181,2011305) #  Year and julian day to downscale

#  Load shapefile with the ground reference tower locations
indir.shp = "G:/mydocuments/SDSU/research/globalET/fluxtowers_USA/fluxtowerALL"
fname.shp = "cropfluxUSlatlon"
shpin = readShapePoints(paste(indir.shp,fname.shp,sep="/"))
shp = shpin[shpin$Location=="California",]

indir.alb = "G:/large_datasets/USA/california/modis_fluxtower_sites/MCD43A3/projected/"
indir.LST = "G:/large_datasets/USA/california/modis_fluxtower_sites/MOD11A1/projected/"
outdir.1030 = "G:/large_datasets/USA/california/modis_fluxtower_sites/SEBAL/Rn/1030am/"
outdir.24 = "G:/large_datasets/USA/california/modis_fluxtower_sites/SEBAL/Rn/24hmean/"

indir.sw = "G:/large_datasets/USA/california/modis_fluxtower_sites/GMAO/"
indir.lw = indir.sw 
indir.lwup = "G:/large_datasets/USA/california/modis_fluxtower_sites/GMAO/LWup/"

# Load albedo list
setwd(indir.alb)
file.list.alb.in = list.files(pattern="WSA")
file.list.alb = file.list.alb.in[grep("shortwave.tif",file.list.alb.in)]
yyyyjjj.albedo = as.numeric(substr(file.list.alb,10,16))
file.list.alb.sub = file.list.alb[(yyyyjjj.albedo>=start.end[1]) & (yyyyjjj.albedo<start.end[2])]
yyyyjjj.alb.sub = as.numeric(substr(file.list.alb.sub,10,16))

# Load example albedo file, serves as a template for resampling MERRA data
x=file.list.alb.sub[1]
balc = raster(x)
balc[balc==32767]=NA

#  Load SW, LW list
patt="nc"
flist.swlw <- list.files(indir.sw,pattern = patt)
dates.swlw = strptime(substr(flist.swlw,37,45),format="%Y%m%d")
yyyyjjj.swlw = as.numeric(format(dates.swlw,"%Y%j"))
index.swlw = which((yyyyjjj.swlw>=start.end[1])&(yyyyjjj.swlw<=start.end[2]))
file.list.swlw.sub = flist.swlw[index.swlw]  # File names in yyyymmdd
yyyyjjj.swlw.sub=yyyyjjj.swlw[index.swlw]

flist.lwup = list.files(indir.lwup,pattern="nc")
dates.lwup = strptime(substr(flist.lwup,37,44),format="%Y%m%d")
yyyyjjj.lwup = as.numeric(format(dates.lwup,"%Y%j"))
index.lwup = which((yyyyjjj.lwup>=start.end[1])&(yyyyjjj.lwup<=start.end[2]))
file.list.lwup.sub = flist.lwup[index.lwup]  # File names in yyyymmdd
yyyyjjj.lwup.sub=yyyyjjj.lwup[index.lwup]

#  Check out the variables in the SW and LW file (MERRA)
setwd(indir.sw)
fid = open.ncdf(file.list.swlw.sub[1])
fid

setwd(indir.lwup)
foo.lwup = raster(file.list.lwup.sub[1])

#  LWGAB = longwave absorbed = LWin
#  LWGEM = longwave emitted = LWout
#  LWGNT = LWGAB - LWGEM

#  MERRA data are in 24 layers per day, layer1 = 0:30 GMT
#  0:00 GMT  =  4pm PST, 5pm PDT
#  MODIS Terra overpass is 1030am PST, which is 630pm GMT, which is layer 19

#  Load SWin and LWin, resample to create stacks of the subsetted data
cnt=0
N = length(file.list.swlw.sub) # For whole list
#N = 5  #  Partial list
for (x in 1:N) {
	setwd(indir.sw)
	bsw.1030am = raster(file.list.swlw.sub[x], band=19, varname = "SWGDN")  # Shortwave at ground down
	# Metadata:  http://globalchange.nasa.gov/KeywordSearch/Metadata.do?Portal=daacs&KeywordPath=[Keyword%3D%27LWGAB%27]&EntryId=GES_DISC_MATUNXRAD_V5.2.0&MetadataView=Text&MetadataType=0&lbnode=mdlb3

	bsw.24 = stack(file.list.swlw.sub[x],varname="SWGDN")
	mean.sw.at.towers = as.data.frame(extract(bsw.24,shp,FUN=mean))
	dayindex = which(mean.sw.at.towers[1,]>0)
	bsw.daymean = round(mean(bsw.24[[dayindex]]),1)
	
	cnt=cnt+1
	print(paste("sw ", cnt," of ", N, sep=""))
  	flush.console()
	projection(bsw.daymean)=projection(balc)
	bsw.day.resam = round(resample(bsw.daymean,balc,method="bilinear"),1)
	if (cnt==1){
		swdaystack = stack(bsw.day.resam)
		}
	if (cnt>1){
		swdaystack = addLayer(swdaystack,bsw.day.resam)
		}
}

writeRaster(swdaystack,paste(outdir.24,"swinWm2_day_",as.character(start.end[1]),"_",as.character(start.end[2]),sep=""),overwrite=TRUE)

rm(swdaystack)
cnt=0
for (x in 1:N) {
	setwd(indir.sw)
	bsw.24 = stack(file.list.swlw.sub[x],varname="SWGDN")
	mean.sw.at.towers = as.data.frame(extract(bsw.24,shp,FUN=mean))
	dayindex = which(mean.sw.at.towers[1,]>0)
	setwd(indir.lwup)
	blw.in.24 = stack(file.list.lwup.sub[x],varname="lwgab")
	blw.in.daymean = round(mean(blw.in.24[[dayindex]]),1)
	blw.out.24 = stack(file.list.lwup.sub[x],varname="lwgem")
	blw.out.daymean = round(mean(blw.out.24[[dayindex]]),1)
	cnt=cnt+1
	print(paste("lw ", cnt," of ", N, sep=""))
  	flush.console()
	projection(blw.in.daymean)=projection(balc)
	projection(blw.out.daymean)=projection(balc)
	blw.in.day.resam = round(resample(blw.in.daymean,balc,method="bilinear"),1)
	blw.out.day.resam = round(resample(blw.out.daymean,balc,method="bilinear"),1)
if (cnt==1){
	lwindaystack = stack(blw.in.day.resam)
	lwoutdaystack = stack(blw.out.day.resam)
	}
if (cnt>1){
	lwindaystack = addLayer(lwindaystack,blw.in.day.resam)
	lwoutdaystack = addLayer(lwoutdaystack,blw.out.day.resam)
		}
}

writeRaster(lwindaystack,paste(outdir.24,"lwinWm2_day_",as.character(start.end[1]),"_",as.character(start.end[2]),sep=""),overwrite=TRUE)
writeRaster(lwoutdaystack,paste(outdir.24,"lwoutWm2_day_",as.character(start.end[1]),"_",as.character(start.end[2]),sep=""),overwrite=TRUE)

## Linear Unmixing of Hyperspectral Point Data

### PRE-PROCESSING RAW DATA ###

#Import data where each row is a data point/coordinate and each column is a wavelength
data<-read.csv("Data.csv")
Wavelength<-read.csv("Wavelength.csv")

#Trim to desired wavelengths (columns) ie. 400-740nm
data1<-data[,c(1:843)]
Wavelength<-Wavelength[,c(1:843)]

# Savitzy Golay smoothing and filtering
library(prospectr)
Spectra<-list()
Spectra$spec<-data1
sg <- savitzkyGolay(Spectra$spec, p = 3, w = 11, m = 0)

#Normalise for vector length

# Define function to find vector length
veclen=function(vec) {
  sqrt(sum(vec^2))
}

# Find vector length for spectrum of each pixel
sgt<-t(sg) #Transpose the dataset so each column is a spectra, each row is the wavelength
sgt<-as.data.frame(sgt)
sgt$Wavelength<-Wavelength

# Define function that return normalized vector
vecnorm=function(vector) {
  vector/veclen(vector)
}

#find vector lengths
veclengths <- sgt %>%
  group_by(Wavelength) %>%
  as.data.frame(summarise(across(sgt[,1:2553], veclen)))

# Normalize by dividing each reflectance value by the vector’s length
df <- veclengths %>%
  group_by(Wavelength) %>%
  summarise(across(veclengths[,1:2553], vecnorm)) #Where there is 2553 data points in this case

dfFinal<-t(df) #transpose
dfFinal<-as.data.frame(dfFinal)

### ENDMEMBER DETERMINATION ###

#Rows are endmembers, columns are wavelengths)
EndmemberLibFull<-read.csv("EndmemberLibFull.csv")
WavelengthEnd<-read.csv("WavelengthEnd")

#Scale endmembers before clustering 
Endscaled<-scale(EndmemberLibFull,center=TRUE,scale=TRUE)
Endscaled<-as.data.frame(Endscaled)

#Principal Component Analysis
library(ggfortify)
pca_res <- prcomp(Endscaled, scale. = TRUE)
Endscaled$Colour<-Class #Can add a list of substrate class groups, ie. Sand, Coral, Rock, Algae) This is to label colours on the plot

#Plot PCA
autoplot(pca_res, data=Endscaled, colour= 'Colour', label=TRUE, label.size = 5)

##Create spectral library
Speclib<-EndmemberLibFull[c(48, 56, 67, 80),] #Pick rows of endmembers wanted (48:Rock, 56:Coral, 67:Algae, 80:Sand)

#Smooth spectral library spectra
SpectraEnd<-list()
SpectraEnd$spec<-Speclib
sgEnd <- savitzkyGolay(SpectraEnd$spec, p = 3, w = 11, m = 0)

#Normalise for vector length
sgtEnd<-t(sgEnd) #Transpose the dataset so each column is a spectra, each row is the wavelength
sgtEnd<-as.data.frame(sgtEnd)
sgtEnd$Wavelength<-WavelengthEnd

#find vector lengths for endmember spectra
veclengthsEnd <- sgtEnd %>%
  group_by(Wavelength) %>%
  as.data.frame(summarise(across(sgtEnd[,1:4], veclen))) #4 endmembers

# Normalize by dividing each reflectance value by the vector’s length
dfEnd <- veclengthsEnd %>%
  group_by(Wavelength) %>%
  summarise(across(veclengthsEnd[,1:4], vecnorm)) #Where there is 4 endmembers

EndmemberdfFinal<-t(dfEnd) #transpose
EndmemberdfFinal<-as.data.frame(EndmemberdfFinal)

#Check for linear dependence

library(plm)
detect.lindep(dfEnd) ##Should return "No linear dependent column(s) detected.", Otherwise, try finding new endmembers that are not collinear

#Resample endmember spectra to drone data wavelengths
library(spectrolab)
EndmemberdfFinal<-as.matrix(EndmemberdfFinal)
EndLibFinal<-resample(EndmemberdfFinal, Wavelength) #Resample through spline smoothing

### SPECTRAL UNMIXING ###
library(RStoolbox)
EndLib<-speclib(EndLibFinal, Wavelength) #Creation of a speclib (format)
dfFinal<-as.matrix(dfFinal) #Must convert to matrix before converting to speclib format
DroneLib<-speclib(dfFinal, Wavelength)

#Unmixing algorithm
Unmix_res<-unmix(DroneLib, EndLib,  returnHCR = "auto", scale = TRUE)

#Extract fractional contributions
fracs<-Unmix_res$fractions
fracst<-t(fracs)
colnames(fracst)<-c("Rock", "Coral", "Algae", "Sand")
fracst<-as.data.frame(fracst)

#Import coordinates from raw data
latlon<-read.csv("coordinates.csv")
fracst$lon<-latlon$lon[c(1:843),] #Trim to the number of datapoints kept in the initial trimming step
fracst$lat<-latlon$lat[c(1:843),]

#Plot fractional contributions of Rock
library(ggplot2)
ggplot10<-ggplot(data = fracst) + geom_point(aes(x = lon, y = lat, colour = fracst$Rock),data = fracst) 
ggplot10+scale_color_gradientn(colours = rainbow(5))

#Plot fractional contributions of Coral
ggplot10<-ggplot(data = fracst) + geom_point(aes(x = lon, y = lat, colour = fracst$Coral),data = fracst) 
ggplot10+scale_color_gradientn(colours = rainbow(5))

#Plot fractional contributions of Algae
ggplot10<-ggplot(data = fracst) + geom_point(aes(x = lon, y = lat, colour = fracst$Algae),data = fracst) 
ggplot10+scale_color_gradientn(colours = rainbow(5))

#Plot fractional contributions of Sand
ggplot10<-ggplot(data = fracst) + geom_point(aes(x = lon, y = lat, colour = fracst$Sand),data = fracst) 
ggplot10+scale_color_gradientn(colours = rainbow(5))

##Export figures on RStudio

#Extract root mean square error (RMSE) values for each coordinate in the unmixing 
error<-fracs$error


### ACCURACY ASSESSMENT ###

#Import results from RGB Classification on Google Earth Engine
RGBClass<-read.csv("RGBClass") #Columns are substrate classes, rows are the randomly generated points, lon and lat are included


fracstACC<-fracst[c(),] #Within c(), input all rows of the randomly generated coordinates (in this study 50 were chosen)

Acc<-merge(fracstACC, RGBClass, by=c("lon","lat"))
colnames(Acc)<-c("Hyp_Rock", "Hyp_Coral", "Hyp_Algae", "Hyp_Sand", "RGB_Rock", "RGB_Coral", "RGB_Algae", "RGB_Sand") #rename columns to suit endmember names corresponding to the test

#Test correlation for rock

Rock_ACCASS<-cor.test(Acc$RGB_Rock, Acc$Hyp_Rock, method="spearman")
Rock_ACCASS #prints results of Spearman's correlation test

#Plot correlation 
ggplot(Acc, aes(x = RGB_Rock, y =Hyp_Rock)) +
  geom_point() +geom_smooth(method='lm',formula= y~x) + labs(x="Rock Cover from RGB Classification (%)", "y" = "Rock Cover from Linear Unmixing (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#Test correlation for Coral

LC_ACCASS<-cor.test(Acc$RGB_Coral, Acc$Hyp_Coral, method="spearman")
LC_ACCASS #prints results of Spearman's correlation test

#Plot correlation 
ggplot(Acc, aes(x = RGB_Coral, y =Hyp_Coral)) +
  geom_point() +geom_smooth(method='lm',formula= y~x) + labs(x="Live Coral Cover from RGB Classification (%)", "y" = "Live Coral Cover from Linear Unmixing (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#Test correlation for algae

Algae_ACCASS<-cor.test(Acc$RGB_Algae, Acc$Hyp_Algae, method="spearman")
Algae_ACCASS #prints results of Spearman's correlation test

#Plot correlation 
ggplot(Acc, aes(x = RGB_Algae, y =Hyp_Algae)) +
  geom_point() +geom_smooth(method='lm',formula= y~x) + labs(x="Algae Cover from RGB Classification (%)", "y" = "Algae Cover from Linear Unmixing (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#Test correlation for sand

Sand_ACCASS<-cor.test(Acc$RGB_Sand, Acc$Hyp_Sand, method="spearman")
Sand_ACCASS #prints results of Spearman's correlation test

#Plot correlation 
ggplot(Acc, aes(x = RGB_Sand, y =Hyp_Sand)) +
  geom_point() +geom_smooth(method='lm',formula= y~x) + labs(x="Sand Cover from RGB Classification (%)", "y" = "Sand Cover from Linear Unmixing (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


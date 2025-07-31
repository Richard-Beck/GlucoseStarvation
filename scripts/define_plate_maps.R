# platemaps are lookup tables informed by the plate map pngs in the corresponding source folders

m <- list()

x <- array(NA, dim = c(6, 10, 2),
           dimnames = list(
             LETTERS[2:7],              
             paste0("C", 2:11),          
             c("ploidy", "glucose")       
           ))
x[,1:5,"ploidy"] <- "parental"
x[,6:10,"ploidy"] <- "3N"

for(i in 1:5){x[1:4,c(i,i+5),"glucose"] <- c(0,0.1,0.25,0.5,1)[i] }
x["F",,"glucose"] <- 5
x["G",,"glucose"] <- 25

m[["MDA-MB-231_B00-IncucyteRawDataLiveDead-varyGlucose-250213"]] <- x

#MCF10A_A00-IncucyteRawDataLiveDead-varyGlucose-241015
x <- array(NA, dim = c(4, 6, 1),
           dimnames = list(
             LETTERS[2:5],              
             paste0("C", 2:7),          
             c("glucose")       
           ))

for(i in 1:7) x[LETTERS[2+(i-1)%%4],(1:3)+3*floor((i-1)/4),"glucose"] <- c(0,0.1,0.25,0.5,1,5,25)[i]


m[["MCF10A_A00-IncucyteRawDataLiveDead-varyGlucose-241015"]] <- x

## these two assumed to have same platemap:
#MCF10A_A00b-IncucyteRawDataLiveDead-varyGlucose-241015
#MCF10A_A00c-IncucyteRawDataLiveDead-varyGlucose-241015
x <- array(NA, dim = c(4, 7, 1),
           dimnames = list(
             LETTERS[2:5],              
             paste0("C", 1:7),          
             c("glucose")       
           ))
for(i in 1:7) x[,i,"glucose"] <- c(0,0.1,0.25,0.5,1,5,25)[i]

m[["MCF10A_A00b-IncucyteRawDataLiveDead-varyGlucose-241015"]] <- x
m[["MCF10A_A00c-IncucyteRawDataLiveDead-varyGlucose-241015"]] <- x



#SNU668_A00-IncucyteRawDataLiveDead-varyGlucose-250324
x <- array(NA, dim = c(6, 10, 2),
           dimnames = list(
             LETTERS[2:7],              
             paste0("C", 2:11),          
             c("ploidy", "glucose")       
           ))
x[,1:5,"ploidy"] <- "low"
x[,6:10,"ploidy"] <- "high"

for(i in 1:5){x[1:4,c(i,i+5),"glucose"] <- c(0,0.1,0.25,0.5,1)[i] }
x["F",,"glucose"] <- 5
x["G",,"glucose"] <- 25
m[["SNU668_A00-IncucyteRawDataLiveDead-varyGlucose-250324"]] <- x

#SUM-159_M00a-IncucyteRawDataLiveDead-varyGlucose and
#SUM-159_M00b-IncucyteRawDataLiveDead-varyGlucose and
x <- array(NA, dim = c(6, 10, 2),
           dimnames = list(
             LETTERS[2:7],              
             paste0("C", 2:11),          
             c("ploidy", "glucose")       
           ))
x[,1:5,"ploidy"] <- "2N"
x[,6:10,"ploidy"] <- "4N"

for(i in 1:5){x[1:4,c(i,i+5),"glucose"] <- c(0,0.1,0.25,0.5,1)[i] }
x["F",,"glucose"] <- 5
x["G",,"glucose"] <- 25
m[["SUM-159_M00a-IncucyteRawDataLiveDead-varyGlucose"]] <- x
m[["SUM-159_M00b-IncucyteRawDataLiveDead-varyGlucose"]] <- x

# "SUM-159_I00-IncucyteRawDataLiveDead-varyGlucose" 
x <- array(NA, dim = c(4, 10,2),
           dimnames = list(
             LETTERS[2:5],              
             paste0("C", 2:11),          
             c("ploidy", "glucose")       
           ))
x[,1:5,"ploidy"] <- "2N"
x[,6:10,"ploidy"] <- "4N"
for(i in 1:5){x[1:4,c(i,i+5),"glucose"] <- c(0,0.1,0.25,0.5,1)[i] }
m[["SUM-159_I00-IncucyteRawDataLiveDead-varyGlucose"]] <- x
#  SUM-159_C00-IncucyteRawDataLiveDead-varyGlucose
x <- array(NA, dim = c(6, 10,2),
           dimnames = list(
             LETTERS[2:7],              
             paste0("C", 2:11),          
             c("ploidy", "glucose")       
           ))
x[,1:5,"ploidy"] <- "2N"
x[,6:10,"ploidy"] <- "4N"
for(i in 1:5){x[1:6,c(i,i+5),"glucose"] <- c(0.1,0.5,1,5,25)[i] }
m[["SUM-159_C00-IncucyteRawDataLiveDead-varyGlucose"]] <- x


get_meta <- function(filename){
  
  fparts <- strsplit(filename,"_") |> unlist()
  pmap <- m[[paste(fparts[1:2],collapse="_")]]
  
  rc <- fparts[length(fparts)-2]
  row <- substr(rc,1,1)
  col <- paste0("C",substr(rc,2,nchar(rc)))
  
  glucose <- pmap[row,col,"glucose"]
  
  ploidy <- NaN
  if(!"ploidy"%in%unlist(dimnames(pmap))) {
    ploidy <- fparts[3]
  } else ploidy <- pmap[row,col,"ploidy"]
  
  tt <- tail(fparts,1) ## expect time in format e.g. "05d00h00m"
  days <- as.numeric(gsub("d.*", "", tt))
  hours <- as.numeric(gsub(".*d|h.*", "", tt))
  mins <- as.numeric(gsub(".*h|m", "", tt))
  
  hours <- days * 24 + hours + mins / 60
  
  if(fparts[1] == "SUM-159"){
    if(fparts[2]=="M00b-IncucyteRawDataLiveDead-varyGlucose"){
      fparts[1] <- paste0(fparts[1],"-chem")
    } else fparts[1] <- paste0(fparts[1],"-fuse") 
  }
  if(fparts[1] == "MCF10A"){
    if(grepl("Ctrl",ploidy)){
      fparts[1] <- paste0(fparts[1],"-ctrl")
      ploidy <- gsub("-Ctrl","",ploidy)
    }
    if(grepl("HRas",ploidy)){
      fparts[1] <- paste0(fparts[1],"-hras")
      ploidy <- gsub("-HRas","",ploidy)
    }
  }
  
  data.frame(cellLine = fparts[1],experiment=fparts[2],plateID=rc,ploidy,glucose,hours)
  
}


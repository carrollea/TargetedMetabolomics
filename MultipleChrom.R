#!/usr/bin/env Rscript --vanilla

#define packages to install
packages <- c('ggplot2', 'dplyr', 'BiocManager')
biopackages <- c('xcms', 'SummarizedExperiment', 'msPurity', 'Spectra')
#install all packages that are not already installed
install.packages(setdiff(packages, rownames(installed.packages())))
BiocManager::install(setdiff(biopackages, rownames(BiocManager::install())))

#The --vanilla on the end, tells Rscript to run without saving or restoring anything in the process. This just keeps things nice a clean.
suppressMessages(library(xcms)) #will make chromatograms 
library(SummarizedExperiment)
suppressMessages(library(msPurity))
suppressMessages(library(Spectra)) #will make spectra
library(dplyr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
#The ‘trailingOnly = TRUE’ option means that the first element of ‘args’ is the first argument, instead of the name of the command itself.

if (length(args) < 4){
  stop("You should have three arguments.\n")
}

#args[1] = directory the mzML files are stored in 
#args[2] = data frame of target ion masses
#args[3] = retention time start
#args[4] = retention time end
#args[5] = data frame describing how the files are related to the sets in the experiment 


invitro<-list.files(path = args[1], pattern = "mzML", full.names = TRUE)

mzlist1<-read.csv(args[2], header = TRUE)
rtr<-c(as.numeric(args[3]), as.numeric(args[4]))
#pd<-read.csv(args[5], header = TRUE)
pd<-data.frame(sample_name = sub(basename(invitro), pattern =".mzML", replacement = "", fixed = TRUE),sample_group = sub(basename(invitro), pattern =".mzML", replacement = "", fixed = TRUE), stringsAsFactors = FALSE)


label_fun <- function(x) {
  ints <- unlist(intensity(x))
  mzs <- format(unlist(mz(x)), digits = 7)
  mzs[ints > 100 & ints < 1000] <- ""
  mzs
}


#set up variables for graph
x<-c(0, 1000)
#y<-c(0, 2000000)
ts = 50 #font size
lw = 2 #width of plotted line
aw = 1  #width of axis line
cl<-"black"
yfs = 20


peak_param = CentWaveParam(peakwidth = c(5, 15), prefilter = c(5, 1000))

Exp_data<-readMSData(files = invitro, pdata = new("NAnnotatedDataFrame", pd), mode = "onDisk")

for (c in 2:length(mzlist1)) {
  for (m in 1:nrow(mzlist1)) {
    mzp<-c(mzlist1[m,c])
    mzp5<-((5*mzp)+(mzp*1000000))/1000000
    mzm5<-((-5*mzp)+(mzp*1000000))/1000000
    mzr<-c(mzm5, mzp5)
    graph_chrom<-chromatogram(Exp_data, aggregationFun = "sum", rt = rtr, mz =mzr) #make chromatograms of all the samples
    xchrs <- findChromPeaks(graph_chrom, param = peak_param) #find peaks in chromatograms
    peakDF<-as.data.frame(chromPeaks(xchrs, rt = rtr)) #or can use pander(chromPeaks(xchrs, rt = rtr)). Make peak information into dataframe
    peakDF$rtminute<-peakDF$rt/60 #find the retention time of peaks in minutes
    peakDF$intmax<-peakDF$maxo + 1000000
    peakDF$ID<-Exp_data$sample_name[peakDF$column]  #put the file name of where the peak is in the dataframe
    if (nrow(peakDF)>0) {
      dirID1<-paste(mzlist1[m,1], mzlist1[m,c], "mz", sep="_") #Make name for new directory of the set that you're analyzing
      dir.create(dirID1) #create a directory with name from above
      peakfile<-file.path(dirID1, paste("Peaks",  ".csv", sep=""))
      write.csv(peakDF, peakfile) #make a file of the data frame
      graph_chrom_clean<-clean(graph_chrom, na.rm = TRUE)
      colid<-unique(peakDF$column) #make a list of the ids associated with each file so we know which files have peaks 
      for (i in 1:length(colid)) {
        peaksub<-subset.data.frame(peakDF, subset = peakDF$column == colid[i]) #only take the peaks from the one file that we want to work with
        peaksub$rtminute<-as.numeric(format(round(peaksub$rtminute, 2), nsmall = 2))
        Test<-as.data.frame(intensity(xchrs[,colid[i]])) #make data frame of the intensity and retention time of the file we are working with. colid ensures we are only getting EIC for files that have peaks
        Test$Rtime<-rtime(xchrs[,colid[i]])
        colnames(Test)<-c("Intensity", "Rtime")
        Final<-Test[order(Test$Rtime),]
        Final$RTimeMin<-Final$Rtime/60 #Get retention time in minutes 
        maxint<-max(Final$Intensity)
        Final[is.na(Final)]=0
        fileeic<-file.path(dirID1, paste(Exp_data$sample_name[colid[i]],"_EIC", ".csv", sep="")) #make file of EIC info
        write.csv(Final, file = fileeic)
        fileid<-file.path(dirID1, paste(Exp_data$sample_name[colid[i]], ".pdf", sep=""))
        temp_plot<-ggplot(Final, aes(y=Intensity, x=RTimeMin)) +
          #  geom_point() +
          geom_line(linewidth = 1) +
          # scale_y_continuous(limits = c(0, 50000000)) +
          #  stat_smooth(aes(y=Intensity, x=RTimeMin), formula = y ~ s(x, k = 75), method = "gam", se = FALSE, col = cl) +
          #  coord_cartesian(xlim = x, ylim = c(0,40000000)) +
          #geom_text(data = peaksub, label = peaksub$rtminute, x = peaksub$rtminute, y = peaksub$maxo + 5000000, size = 8) +  #Add labels to the graph. can turn off if wanted or use geom_label for a label with a background
          labs(x = "Retention Time (Min)") +
          theme(panel.background = element_rect(fill='transparent'), plot.background = element_rect(fill='transparent', color=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill='transparent'), legend.box.background = element_rect(fill='transparent'), text = element_text(size = 26), axis.line = element_line(color="black", size = 1)) #this makes the background transparent 
        ggsave(temp_plot, file=fileid, width = 20, height = 10, units = "cm") #save the ggplot
      }
      for (i in 1:length(colid)) {
        Sp<-Spectra(invitro[colid[i]]) #get spectra from one file
        dirID2<-paste(Exp_data$sample_name[colid[i]],"_Spectra", sep="") #Make name for new directory of the set that you're analyzing
   #     dir.create(file.path(dirID1, dirID2)) #create a directory with name from above
        peaksub<-subset.data.frame(peakDF, subset = peakDF$column == colid[i]) #get info for peaks from the peakDF file for only our set of interest
        for (a in 1:nrow(peaksub)) {
          TargetRT<-c(peaksub[a,2], peaksub[a,3]) #selects retention times min and max around the peak that you want to grab the ms2 spectra of
          SpFil<-filterRt(Sp, rt=TargetRT)#filter the spectra you want based on TargetRT range
          SpFilmz<-filterPrecursorMzValues(Sp, mz=mzp, ppm=5) #further filter SpFil spectra for spectra with specific parent ion
          if (length(SpFilmz)>0) {
            for (y in 1:length(SpFilmz)) {
              dir.create(file.path(dirID1, dirID2)) #create a directory with name from above
              intlist<-unlist(intensity(SpFilmz[y])) #for the filtered spectra get only the intensity
              mzlist<-unlist(mz(SpFilmz[y])) #for the filtered spectra get only the m/z values
              SpFilmzData<-data.frame(MZ=mzlist, Intensity=intlist) #make dataframe of intensity and m/z values
              SpecrawID<-file.path(dirID1, dirID2, paste(Exp_data$sample_name[colid[i]], rtime(SpFilmz[y]), ".csv", sep="_")) #make a .csv file of the intensity and m/z dataframe
              write.csv(SpFilmzData, SpecrawID)
              SpecID<-file.path(dirID1, dirID2, paste(Exp_data$sample_name[colid[i]], rtime(SpFilmz[y]), ".pdf", sep="_")) #make a name for the pdf of the filtered spectra
              pdf(SpecID, width = 8, height = 6) #make a pdf of the filtered spectra
              plotSpectra(SpFilmz[y], ppm = 5, labels = label_fun, labelPos = 2, labelOffset = 0.4, labelSrt = -30, labelCex = 0.8)
              dev.off()
            }
          }
        }
      }
    }
  }
}

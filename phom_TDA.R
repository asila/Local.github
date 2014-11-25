#install.packages("phom")
#library(lle)
library(phom)
library(lsa)
ir<-read.csv(file.choose())
#ir<-ir[,-c(2:10)]
#k<-calc_k(ir[,-1],m=2)
#lle.m<-lle(nir.cos,m=2,k=16)
mir<-ir[,-1]

mir.cos<-round(cosine(t(mir)),3)
length(mir.cos)
#Homology computation
mir.ph<-pHom(mir.cos,6,max_filtration_value=0.6,metric="distance_matrix",mode="lw",maxmin_samples=nrow(mir.cos))


plotPersistenceDiagram(mir.ph,max_dim=3,max_f=0.5,title="MIR Persistence Diagram")


#############################################################3
##  Newbies Guide to Homology - The Pursuit of Topological Data Analysis
## 
##  Run example code from PHOM, and then create a shape that is
##  not micky mouse, and create Betti Barcodes on them too
##  
##  Ryan Anderson - March 2014 - Dream to Lear Data Viz Community
##  Enjoy!  (and thanks to authors of package and example group 1 )
##############################################################
##
##
install.packages("phom")
library(phom)

#   EXAMPLE GROUP 1 - Source: http://cran.at.r-project.org/web/packages/phom/phom.pdf examples sourced
#   Title Persistent Homology in R > Version 1.0.3 Date 2011-11-02 Author Andrew Tausz
#   And see also http://comptop.stanford.edu/ 
### This generates a simple set of points, circular, so we ca see Betti Barcodes evolve
t <- 2 * pi * runif(30)
x <- cos(t); y <- sin(t)
M <- t(as.matrix(rbind(x, y)))
plot(M)
max_dim <- 1
max_f <-  .6
intervals <- pHom(M, max_dim, max_f)
plotPersistenceDiagram(intervals, max_dim, max_f,
                       title="Random Points on S^1 with Euclidean Norm")
#BETTI 0 - every point gets 'represented' until the bubbles join them all up
plotBarcodeDiagram(intervals, 0, max_f, title="") # with small Epsilon, every point gets represented, as epsilon 'bubbles' grow - our N points get lumped into one 'item' (from many to one)
#BETTI 1 - only when the last 'join' happens, does it 'become' a circle.  If you increase max_F to 2 or so (which crashes my R) the cicle becomes a blob and the band in Betti 1 ends just before 2
plotBarcodeDiagram(intervals, 1, max_f, title="") # the circle becomes a coherent loop (joined up) at right point, and then 'holds together'

#### OK - let's tweak values
max_dim <- 1
max_f <-  2 ## this will make R work hard  
# increased this from 0.6 to see if we could get Betti 1 to collapse as we filled in center of circle
# above 1 it takes a while to run - 2 looks like end of the line around 1.8 (the circle got 'filled in)
intervals <- pHom(M, max_dim, max_f)
plotBarcodeDiagram(intervals, 0, max_f, title="") # with small Epsilon, every point gets represented, as epsilon 'bubbles' grow - our N points get lumped into one 'item' (from many to one)
plotBarcodeDiagram(intervals, 1, max_f, title="") # the circle becomes a coherent loop (joined up) at right point, and then 'holds together'

plotPersistenceDiagram(intervals, max_dim, max_f,
                       title="Random Points on S^1 with Euclidean Norm")
#######################
## Example GROUP 2: let's get fancy with our points and create a shape that is 
## not a mouse
t <- 2 * pi * runif(100)
x <- cos(t); y <- sin(t)
HEAD <- t(as.matrix(rbind(x, y)))
HEAD <- HEAD + 3
plot(HEAD)
HEAD
right_ear <- (HEAD / 1.5)+2.2
plot(right_ear)
left_ear <- right_ear
plot(left_ear)
for (i in 1:length(right_ear)) 
    { 
    left_ear[i,1] = left_ear[i,1]-2.4
    }
plot(left_ear)
not_mickey_mouse <- rbind(HEAD,right_ear,left_ear)
plot(not_mickey_mouse, main = "Not a mouse")
max_dim <- 1
max_f <-  .6
intervals2 <- pHom(not_mickey_mouse, max_dim, max_f)
plotPersistenceDiagram(intervals2, max_dim, max_f,
                       title="Not a Mouse Persistence Diagram - Euclidean Norm")
#BETTI 0 - every point gets 'represented' until the bubbles join them all up
plotBarcodeDiagram(intervals2, 0, max_f, title="Betti 0 Barcode - Not a Mouse") # with small Epsilon, every point gets represented, as epsilon 'bubbles' grow - our N points get lumped into one 'item' (from many to one)
#BETTI 1 - only when the last 'join' happens, does it 'become' a circle.  If you increase max_F to 2 or so (which crashes my R) the cicle becomes a blob and the band in Betti 1 ends just before 2
plotBarcodeDiagram(intervals2, 1, max_f, title="Betti 1 Barcode - Not a Mouse") # the circle becomes a coherent loop (joined up) at right point, and then 'holds together'

##### More reading: 
#  http://en.wikipedia.org/wiki/Betti_number - nice 
#  "for higher dimensions "In many biological settings, Betti numbers are used 
#  to understand the properties of genes located in breast cancer, by creating a curve 
# of Betti numbers. For gene expression and copy number data sets, a graph is created 
# out of one patients' log{2} ratios (y-axis), calculated in a DNA microarray, and the
# location of these log2 ratios in a specific chromosome (x-axis). Using a window of 
# size 1,2,3,....,or n-dimensions, a point cloud can be constructed...."
# The filtration is denoted by  \epsilon , which is a very small number, usually
# ranging from .0000001 to .1.  \epsilon  is now considered the RADIUS OF A CIRCLE
# centered at each point in the point cloud. When two circles OVERLAP, this forms a
# connection between the two points, as this process is continued, more simplices will
# show up with more Betti numbers as well. As  \epsilon  increases, there are more 
# "births" and "deaths" in the data, meaning that as the filtration changes, certain
#Betti numbers will decrease and others will increase.  SOURCE: Wikipedia Betti_number
# See also:  http://www.youtube.com/watch?v=dDQgA8Qbaho - helps visualize Barcodes

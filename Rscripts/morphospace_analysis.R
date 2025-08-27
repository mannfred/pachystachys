library(geomorph)
library(here)
library(tidyverse)

#import data
paths <- list.files(path=here("data/raw_data/"), pattern="*.TPS", full.names=TRUE, recursive=TRUE) 
lms <- readmulti.tps(paths, specID="imageID", readcurves=TRUE)




# -----------------------------------------------------------
# estimate bill arclengths (from non-procrustes-aligned data)
# don't need to scale because geomorph::readland.tps already scaled landmarks

# index the distances between landmarks placed on the dorsal side of the bill
# `dorsal_distances` is a matrix with rows as species and columns as landmark distances
lmks <- matrix(c(1:10, 2:11), ncol=2)
colnames(lmks) <- c('start', 'end')
dorsal_distances <- geomorph::interlmkdist(A=lms, lmks=lmks)


# calculate arclength for each specimen by summing the dorsal landmark distances
arclengths <- 
  dorsal_distances |> 
  as.data.frame() |> 
  tibble::rownames_to_column(var = "species") |> 
  dplyr::rowwise() |> 
  dplyr::mutate(arclength = sum(c_across(V1:V10))) |> 
  dplyr::mutate(species = str_remove(species, ".PNG")) |>
  dplyr::select(species, bill_arclength)


# -------------------------------------------------------------------
# organize landmark data for curvr
# subset to dorsal landmarks (LMs 1-11, LM11 is the apex of the bill)
# Momocs::Ldk() converts to Coo object
dorsal <- Momocs::Ldk(lms[c(1:11),,])

# align the longest axis of each shape along the x-axis
dorsal_x <- Momocs::coo_alignxax(dorsal)

# simplify nested coo objects
dorsal_xsimple <- unlist(dorsal_x, recursive=F) 

# x axis boundaries
dorsal_baselines <- 
  dorsal_xsimple %>% 
  lapply(., function(b) c(unlist(b[,1])[1], unlist(b[,1])[11]))


# fit interpolating splines and compute total curvature 
dorsal_k_spline <-
  mapply(curvr::curvature_spline, dorsal_xsimple, dorsal_baselines, 'smooth') |> 
  tibble::as_tibble() |>  
  dplyr::slice(1) |> 
  unlist() |> 
  tibble::as_tibble() |>  
  dplyr::mutate(Ktot = abs(value*(180/pi))) |> 
  dplyr::mutate(species = sub(" ", "_", dimnames(lms)[[3]])) |> 
  dplyr::mutate(species = str_remove(species, ".PNG")) 
  

# visualize fitted splines
plot(lms[,,6])
lines(smooth.spline(lms[,,6]), lwd=2, col='red')


# ---------------------------------------------
# stage dataset and plot

plotdata <- dplyr::left_join(arclengths, dorsal_k_spline)
plot(plotdata$bill_arclength, plotdata$Ktot)
text(plotdata$bill_arclength, plotdata$Ktot, labels=plotdata$species)

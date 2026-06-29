###############################################################################################################################
##
##  (02) This script explores the issue of depth stratified anisotropy for the Chatham Rise case
##
###############################################################################################################################

######## Load required libraries
library( tinyVAST )
library( sf )
library( fmesher )
library( ggplot2 )
library( tidyverse )
library( dplyr )
library( profvis )
library( sp )
library( rnaturalearth )
library( geosphere )
library( GISTools )
library( raster ) 
library( ggspatial )
library( Kendall )
library( pdp )

######## Load utility functions
source( make.filename( "utilities.R", DIR$Functions ) )
source( make.filename( "plot utility functions.R", DIR$Functions ) )

######## Define other functions
#### Define the "as.matrix.sgbp" function 
as.matrix.sgbp <- function( x, ... ) {

	nc <- attr( x, "ncol" )
	get_vec <- function( x, n ) { v = rep( FALSE, n ); v[x] = TRUE; v }
	do.call( rbind, lapply( x, get_vec, n = nc ) )

}

######## Set up the folder where the maps produced within this script will be saved
MapsDir <- paste0( DIR$Maps, "Chatham Rise case" )
if ( !dir.exists( MapsDir ) ) { dir.create( MapsDir ) } 

######## Load and process shapefiles
contour_200m <- read_sf( dsn = make.filename( "NZ_200_m.shp", DIR$Input ) )
crs_contour_200m <- st_crs( contour_200m )
contour_200m <- st_transform( contour_200m, crs = "+init=epsg:4326" )
contour_200m <- sf::st_shift_longitude( contour_200m )
contour_500m <- read_sf( dsn = make.filename( "NZ_500_m.shp", DIR$Input ) )
contour_500m <- st_transform( contour_500m, crs = "+init=epsg:4326" )
contour_500m <- sf::st_shift_longitude( contour_500m )
depth_contours <- read_sf( dsn = make.filename( "NZ_region_bathymetry_contours.shp", DIR$Input ) )
st_crs( depth_contours ) <- crs_contour_200m
contour_800m <- depth_contours %>% filter( ELEV == -800 )
contour_800m <- st_transform( contour_800m, crs = "+init=epsg:4326" )
contour_800m <- sf::st_shift_longitude( contour_800m )

######## Load mapping parameters 
load( file = make.filename( "Mapping_parameters_Chatham_Rise_hake.RData", DIR$Input ) )

####### Define or modify a few parameters/settings
location_1 <- "Shallow Chatham Rise"
location_2 <- "Middle-depth Chatham Rise"
location_3 <- "Deep Chatham Rise"
location_of_interest_1 <- c( 183.2, -43.8 )
location_of_interest_2 <- c( 181, -44 )
location_of_interest_3 <- c( 179, -45 )
Mapping_parameters[["legend_x"]] <- c( 0.84, 0.89 ) 
Mapping_parameters[["legend_y"]] <- c( 0.05, 0.25 ) 
Mapping_parameters[["xlim_map"]] <- c( 172.5, 187.0 )
map_sf <- Mapping_parameters[["land_data"]]
Mapping_parameters[["xaxis"]] <- seq( from = 170, to = 190, by = 5 )
Mapping_parameters[["xaxis_label"]] <- as.character( Mapping_parameters[["xaxis"]] ) 
Mapping_parameters[["yaxis"]] <- seq( from = -46, to = -42, by = 1 )
Mapping_parameters[["yaxis_label"]] <- as.character( Mapping_parameters[["yaxis"]] ) 
cutoff <- 0.1
barrier_stiffness <- 1e-16
cellsize <- 0.1 * c( 1, 1 )

######## Load and process data 
load( file = make.filename( "Chatham_Rise_hake_groomed_survey_data.RData", DIR$Input ) )
Data <- Survey_data
Data$Year <- Data$fyear 
Data <- Data[,c( 'Year', 'Lon', 'Lat', 'Catch_rate' )]
Data$var <- "density"

######## Define settings for tinyVAST models
formula <- Catch_rate ~ 1
space_term <- "
"
spacetime_term <- "
" 
time_term <- "
" 

######## Generate a mesh for the tinyVAST models using the SPDE method and covariate-based regional anisotropy 
mesh <- fm_mesh_2d( Data[,c( "Lon", "Lat" )], cutoff = cutoff )
nn <- RANN::nn2( data = Data[,c( 'Lon', 'Lat' )], query = mesh$loc[,1 : 2], k = 1 )

Data_vect <- terra::vect( Data, geom = c( "Lon", "Lat" ), crs = "EPSG:4326" )
Data_vect_projected <- terra::project( Data_vect, "EPSG:9191" )
inputfile <- "nzbathymetry_2016_Albers1km_SCC.nc"
file_name <- make.filename( inputfile, DIR$Input ) 
Env_data_raster <- terra::rast( file_name ) 
crs( Env_data_raster ) <- "EPSG:9191"
Bathymetry_extract <- terra::extract( Env_data_raster, Data_vect_projected )$bathy
Depth <- -Bathymetry_extract
Depth <- ifelse( Depth < 0, 0, Depth )
summary( Depth )
hist( Depth )
pause( 3 ); graphics.off()

Data$Depth <- Depth

mesh_cov <- add_mesh_covariates(
  	mesh = mesh,
  	data = Data,
  	covariates = c(),
  	coords = c( "Lon", "Lat" )
)
mesh_cov$vertex_covariates <- data.frame(
  	x = ( mesh$loc[,1] ),
  	y = ( mesh$loc[,2] )
)
mesh_cov$vertex_covariates <- cbind(
  	mesh_cov$vertex_covariates,
	Depth = Data[nn$nn.idx, 'Depth']
)

plot( mesh_cov )
pause( 3 ); graphics.off()

######## Prepare the tinyVAST models using SPDE method and covariate-based regional anisotropy 
fit_norun <- tinyVAST( 
	data = droplevels( Data ), 
	formula = formula, 
	space_term = space_term,
  	spacetime_term = spacetime_term,
  	time_term = time_term,
  	spatial_domain = mesh_cov,
  	family = tweedie( "log" ),
	space_columns = c( "Lon", "Lat" ),
  	time_column = "Year",
  	times = min( Data$Year ) : max( Data$Year ),
  	variable_column = "var",
  	development = list(
		vertex_formula = ~ Depth
  	),
	control = tinyVASTcontrol( 
		verbose = TRUE, 
		opt_loops = 0,
		newton_loops = 0, 
		trace = 1,
    		use_anisotropy = FALSE, 
		getsd = FALSE,
    		calculate = FALSE
	)
)

######## Fit the tinyVAST model using the SPDE method and covariate-based regional anisotropy under Scenario 1
tmb_par <- fit_norun$tmb_inputs$tmb_par
tmb_par$ln_H_input[3] <- 0
tmb_par$log_kappa <- 0

Regional_Anisotropy_Scenario1 <- tinyVAST( 
	data = droplevels( Data ), 
	formula = formula, 
	space_term = space_term,
  	spacetime_term = spacetime_term,
  	time_term = time_term,
  	spatial_domain = mesh_cov,
  	family = tweedie( "log" ),
	space_columns = c( "Lon", "Lat" ),
  	time_column = "Year",
  	times = min( Data$Year ) : max( Data$Year ),
  	variable_column = "var",
  	development = list(
		vertex_formula = ~ Depth
  	),
	control = tinyVASTcontrol( 
		verbose = TRUE, 
		opt_loops = 0,
		newton_loops = 0, 
		trace = 1,
    		use_anisotropy = FALSE, 
		getsd = FALSE,
    		calculate = FALSE, 
		tmb_par = tmb_par
	)
)

######## Approximate spatial correlation for the tinyVAST model using the SPDE method and covariate-based regional anisotropy 
######## for the first location of interest under Scenario 1
mesh_sf <- fm_as_sfc( mesh )
domain_sf <- st_union( mesh_sf )
grid_sf <- st_make_grid( domain_sf, cellsize = cellsize )
grid_sf <- st_intersection( grid_sf, domain_sf )
pred <- st_coordinates( st_centroid( grid_sf ) )
fit <- Regional_Anisotropy_Scenario1
r1 <- spatial_cor(
  	fit$rep$Q,
  	mesh = mesh,
  	coord = location_of_interest_1,
  	pred = pred
)
plot_sf <- st_sf( st_geometry( grid_sf ), r1 = r1 )
r_Regional_Anisotropy_Scenario1_1 <- plot_sf

plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
pause( 3 ); graphics.off()

######## Approximate spatial correlation for the tinyVAST model using the SPDE method and covariate-based regional anisotropy 
######## for the second location of interest under Scenario 1
mesh_sf <- fm_as_sfc( mesh )
domain_sf <- st_union( mesh_sf )
grid_sf <- st_make_grid( domain_sf, cellsize = cellsize )
grid_sf <- st_intersection( grid_sf, domain_sf )
pred <- st_coordinates( st_centroid( grid_sf ) )
fit <- Regional_Anisotropy_Scenario1
r1 <- spatial_cor(
  	fit$rep$Q,
  	mesh = mesh,
  	coord = location_of_interest_2,
  	pred = pred
)
plot_sf <- st_sf( st_geometry( grid_sf ), r1 = r1 )
r_Regional_Anisotropy_Scenario1_2 <- plot_sf

plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
pause( 3 ); graphics.off()

######## Approximate spatial correlation for the tinyVAST model using the SPDE method and covariate-based regional anisotropy 
######## for the third location of interest under Scenario 1
mesh_sf <- fm_as_sfc( mesh )
domain_sf <- st_union( mesh_sf )
grid_sf <- st_make_grid( domain_sf, cellsize = cellsize )
grid_sf <- st_intersection( grid_sf, domain_sf )
pred <- st_coordinates( st_centroid( grid_sf ) )
fit <- Regional_Anisotropy_Scenario1
r1 <- spatial_cor(
  	fit$rep$Q,
  	mesh = mesh,
  	coord = location_of_interest_3,
  	pred = pred
)
plot_sf <- st_sf( st_geometry( grid_sf ), r1 = r1 )
r_Regional_Anisotropy_Scenario1_3 <- plot_sf

plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
pause( 3 ); graphics.off()

######## Fit the tinyVAST model using the SPDE method and covariate-based regional anisotropy under Scenario 2
tmb_par <- fit_norun$tmb_inputs$tmb_par
tmb_par$ln_H_input[3] <- 5
tmb_par$log_kappa <- 7

Regional_Anisotropy_Scenario2 <- tinyVAST( 
	data = droplevels( Data ), 
	formula = formula, 
	space_term = space_term,
  	spacetime_term = spacetime_term,
  	time_term = time_term,
  	spatial_domain = mesh_cov,
  	family = tweedie( "log" ),
	space_columns = c( "Lon", "Lat" ),
  	time_column = "Year",
  	times = min( Data$Year ) : max( Data$Year ),
  	variable_column = "var",
  	development = list(
		vertex_formula = ~ Depth
  	),
	control = tinyVASTcontrol( 
		verbose = TRUE, 
		opt_loops = 0,
		newton_loops = 0, 
		trace = 1,
    		use_anisotropy = FALSE, 
		getsd = FALSE,
    		calculate = FALSE, 
		tmb_par = tmb_par
	)
)

######## Approximate spatial correlation for the tinyVAST model using the SPDE method and covariate-based regional anisotropy 
######## for the first location of interest under Scenario 2
mesh_sf <- fm_as_sfc( mesh )
domain_sf <- st_union( mesh_sf )
grid_sf <- st_make_grid( domain_sf, cellsize = cellsize )
grid_sf <- st_intersection( grid_sf, domain_sf )
pred <- st_coordinates( st_centroid( grid_sf ) )
fit <- Regional_Anisotropy_Scenario2
r1 <- spatial_cor(
  	fit$rep$Q,
  	mesh = mesh,
  	coord = location_of_interest_1,
  	pred = pred
)
plot_sf <- st_sf( st_geometry( grid_sf ), r1 = r1 )
r_Regional_Anisotropy_Scenario2_1 <- plot_sf

plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
pause( 3 ); graphics.off()

######## Approximate spatial correlation for the tinyVAST model using the SPDE method and covariate-based regional anisotropy 
######## for the second location of interest under Scenario 2
mesh_sf <- fm_as_sfc( mesh )
domain_sf <- st_union( mesh_sf )
grid_sf <- st_make_grid( domain_sf, cellsize = cellsize )
grid_sf <- st_intersection( grid_sf, domain_sf )
pred <- st_coordinates( st_centroid( grid_sf ) )
fit <- Regional_Anisotropy_Scenario2
r1 <- spatial_cor(
  	fit$rep$Q,
  	mesh = mesh,
  	coord = location_of_interest_2,
  	pred = pred
)
plot_sf <- st_sf( st_geometry( grid_sf ), r1 = r1 )
r_Regional_Anisotropy_Scenario2_2 <- plot_sf

plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
pause( 3 ); graphics.off()

######## Approximate spatial correlation for the tinyVAST model using the SPDE method and covariate-based regional anisotropy 
######## for the third location of interest under Scenario 2
mesh_sf <- fm_as_sfc( mesh )
domain_sf <- st_union( mesh_sf )
grid_sf <- st_make_grid( domain_sf, cellsize = cellsize )
grid_sf <- st_intersection( grid_sf, domain_sf )
pred <- st_coordinates( st_centroid( grid_sf ) )
fit <- Regional_Anisotropy_Scenario2
r1 <- spatial_cor(
  	fit$rep$Q,
  	mesh = mesh,
  	coord = location_of_interest_3,
  	pred = pred
)
plot_sf <- st_sf( st_geometry( grid_sf ), r1 = r1 )
r_Regional_Anisotropy_Scenario2_3 <- plot_sf

plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
pause( 3 ); graphics.off()

######## Fit the tinyVAST model using the SPDE method and covariate-based regional anisotropy under Scenario 3
tmb_par <- fit_norun$tmb_inputs$tmb_par
tmb_par$ln_H_input[3] <- 25
tmb_par$log_kappa <- 7

Regional_Anisotropy_Scenario3 <- tinyVAST( 
	data = droplevels( Data ), 
	formula = formula, 
	space_term = space_term,
  	spacetime_term = spacetime_term,
  	time_term = time_term,
  	spatial_domain = mesh_cov,
  	family = tweedie( "log" ),
	space_columns = c( "Lon", "Lat" ),
  	time_column = "Year",
  	times = min( Data$Year ) : max( Data$Year ),
  	variable_column = "var",
  	development = list(
		vertex_formula = ~ Depth
  	),
	control = tinyVASTcontrol( 
		verbose = TRUE, 
		opt_loops = 0,
		newton_loops = 0, 
		trace = 1,
    		use_anisotropy = FALSE, 
		getsd = FALSE,
    		calculate = FALSE, 
		tmb_par = tmb_par
	)
)

######## Approximate spatial correlation for the tinyVAST model using the SPDE method and covariate-based regional anisotropy 
######## for the first location of interest under Scenario 3
mesh_sf <- fm_as_sfc( mesh )
domain_sf <- st_union( mesh_sf )
grid_sf <- st_make_grid( domain_sf, cellsize = cellsize )
grid_sf <- st_intersection( grid_sf, domain_sf )
pred <- st_coordinates( st_centroid( grid_sf ) )
fit <- Regional_Anisotropy_Scenario3
r1 <- spatial_cor(
  	fit$rep$Q,
  	mesh = mesh,
  	coord = location_of_interest_1,
  	pred = pred
)
plot_sf <- st_sf( st_geometry( grid_sf ), r1 = r1 )
r_Regional_Anisotropy_Scenario3_1 <- plot_sf

plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
pause( 3 ); graphics.off()

######## Approximate spatial correlation for the tinyVAST model using the SPDE method and covariate-based regional anisotropy 
######## for the second location of interest under Scenario 3
mesh_sf <- fm_as_sfc( mesh )
domain_sf <- st_union( mesh_sf )
grid_sf <- st_make_grid( domain_sf, cellsize = cellsize )
grid_sf <- st_intersection( grid_sf, domain_sf )
pred <- st_coordinates( st_centroid( grid_sf ) )
fit <- Regional_Anisotropy_Scenario3
r1 <- spatial_cor(
  	fit$rep$Q,
  	mesh = mesh,
  	coord = location_of_interest_2,
  	pred = pred
)
plot_sf <- st_sf( st_geometry( grid_sf ), r1 = r1 )
r_Regional_Anisotropy_Scenario3_2 <- plot_sf

plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
pause( 3 ); graphics.off()

######## Approximate spatial correlation for the tinyVAST model using the SPDE method and covariate-based regional anisotropy 
######## for the third location of interest under Scenario 3
mesh_sf <- fm_as_sfc( mesh )
domain_sf <- st_union( mesh_sf )
grid_sf <- st_make_grid( domain_sf, cellsize = cellsize )
grid_sf <- st_intersection( grid_sf, domain_sf )
pred <- st_coordinates( st_centroid( grid_sf ) )
fit <- Regional_Anisotropy_Scenario3
r1 <- spatial_cor(
  	fit$rep$Q,
  	mesh = mesh,
  	coord = location_of_interest_3,
  	pred = pred
)
plot_sf <- st_sf( st_geometry( grid_sf ), r1 = r1 )
r_Regional_Anisotropy_Scenario3_3 <- plot_sf

plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
pause( 3 ); graphics.off()

######## Fit the tinyVAST model using the SPDE method and covariate-based regional anisotropy under Scenario 4
tmb_par <- fit_norun$tmb_inputs$tmb_par
tmb_par$ln_H_input[3] <- 50
tmb_par$log_kappa <- 7

Regional_Anisotropy_Scenario4 <- tinyVAST( 
	data = droplevels( Data ), 
	formula = formula, 
	space_term = space_term,
  	spacetime_term = spacetime_term,
  	time_term = time_term,
  	spatial_domain = mesh_cov,
  	family = tweedie( "log" ),
	space_columns = c( "Lon", "Lat" ),
  	time_column = "Year",
  	times = min( Data$Year ) : max( Data$Year ),
  	variable_column = "var",
  	development = list(
		vertex_formula = ~ Depth
  	),
	control = tinyVASTcontrol( 
		verbose = TRUE, 
		opt_loops = 0,
		newton_loops = 0, 
		trace = 1,
    		use_anisotropy = FALSE, 
		getsd = FALSE,
    		calculate = FALSE, 
		tmb_par = tmb_par
	)
)

######## Approximate spatial correlation for the tinyVAST model using the SPDE method and covariate-based regional anisotropy 
######## for the first location of interest under Scenario 4
mesh_sf <- fm_as_sfc( mesh )
domain_sf <- st_union( mesh_sf )
grid_sf <- st_make_grid( domain_sf, cellsize = cellsize )
grid_sf <- st_intersection( grid_sf, domain_sf )
pred <- st_coordinates( st_centroid( grid_sf ) )
fit <- Regional_Anisotropy_Scenario4
r1 <- spatial_cor(
  	fit$rep$Q,
  	mesh = mesh,
  	coord = location_of_interest_1,
  	pred = pred
)
plot_sf <- st_sf( st_geometry( grid_sf ), r1 = r1 )
r_Regional_Anisotropy_Scenario4_1 <- plot_sf

plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
pause( 3 ); graphics.off()

######## Approximate spatial correlation for the tinyVAST model using the SPDE method and covariate-based regional anisotropy 
######## for the second location of interest under Scenario 4
mesh_sf <- fm_as_sfc( mesh )
domain_sf <- st_union( mesh_sf )
grid_sf <- st_make_grid( domain_sf, cellsize = cellsize )
grid_sf <- st_intersection( grid_sf, domain_sf )
pred <- st_coordinates( st_centroid( grid_sf ) )
fit <- Regional_Anisotropy_Scenario4
r1 <- spatial_cor(
  	fit$rep$Q,
  	mesh = mesh,
  	coord = location_of_interest_2,
  	pred = pred
)
plot_sf <- st_sf( st_geometry( grid_sf ), r1 = r1 )
r_Regional_Anisotropy_Scenario4_2 <- plot_sf

plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
pause( 3 ); graphics.off()

######## Approximate spatial correlation for the tinyVAST model using the SPDE method and covariate-based regional anisotropy 
######## for the third location of interest under Scenario 4
mesh_sf <- fm_as_sfc( mesh )
domain_sf <- st_union( mesh_sf )
grid_sf <- st_make_grid( domain_sf, cellsize = cellsize )
grid_sf <- st_intersection( grid_sf, domain_sf )
pred <- st_coordinates( st_centroid( grid_sf ) )
fit <- Regional_Anisotropy_Scenario4
r1 <- spatial_cor(
  	fit$rep$Q,
  	mesh = mesh,
  	coord = location_of_interest_3,
  	pred = pred
)
plot_sf <- st_sf( st_geometry( grid_sf ), r1 = r1 )
r_Regional_Anisotropy_Scenario4_3 <- plot_sf

plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
pause( 3 ); graphics.off()

######## Produce a map of the approximated spatial correlation for the tinyVAST model using the SPDE method 
######## and covariate-based regional anisotropy under Scenario 1 for the first location of interest and save it
scenario <- "Regional_Anisotropy_Scenario1"
plot_sf <- r_Regional_Anisotropy_Scenario1_1
zlim <- range( c( r_Regional_Anisotropy_Scenario1_1$r1, r_Regional_Anisotropy_Scenario1_2$r1, 
	r_Regional_Anisotropy_Scenario1_3$r1, 
	r_Regional_Anisotropy_Scenario2_1$r1, r_Regional_Anisotropy_Scenario2_2$r1, 
	r_Regional_Anisotropy_Scenario2_3$r1, 
	r_Regional_Anisotropy_Scenario3_1$r1, r_Regional_Anisotropy_Scenario3_2$r1, 
	r_Regional_Anisotropy_Scenario3_3$r1, 
	r_Regional_Anisotropy_Scenario4_1$r1, r_Regional_Anisotropy_Scenario4_2$r1, 
	r_Regional_Anisotropy_Scenario4_3$r1 ), na.rm = TRUE ) 

#### Produce a version of the map with a legend and save it
new.graph( 0.5 )
par( oma = c( 3, 3, 0, 0 ), cex = 1.1 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, main = "", axes = FALSE, 
	xlim = st_bbox(plot_sf)[c( 1, 3 )], ylim = st_bbox(plot_sf)[c( 2, 4 )] ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
plot( contour_200m, col = "white", lwd = 3.5, lty = 2, add = TRUE )
plot( contour_500m, col = "gray48", lwd = 3.5, lty = 2, add = TRUE )
box()
axis( side = 1, at = Mapping_parameters[["xaxis"]], labels = Mapping_parameters[["xaxis_label"]] )
axis( side = 2, at = Mapping_parameters[["yaxis"]], labels = Mapping_parameters[["yaxis_label"]] )
minor.tick( nx = 2, ny = 1 )
mtext( side = 1, outer = TRUE, "Longitude (°E)", cex = 1.25, line = par()$oma[1]/2, padj = 0 )
mtext( side = 2, outer = TRUE, "Latitude (°N)", cex = 1.25, line = par()$oma[2]/2, padj = 0, las = 0 )
xl = ( 1 - Mapping_parameters[["legend_x"]][1] ) * par( 'usr' )[1] + 
	( Mapping_parameters[["legend_x"]][1] ) * par( 'usr' )[2]
xr = ( 1 - Mapping_parameters[["legend_x"]][2] ) * par( 'usr' )[1] + 
	( Mapping_parameters[["legend_x"]][2] ) * par( 'usr' )[2]
yb = ( 1 - Mapping_parameters[["legend_y"]][1] ) * par( 'usr' )[3] + 
	( Mapping_parameters[["legend_y"]][1] ) * par( 'usr' )[4]
yt = ( 1 - Mapping_parameters[["legend_y"]][2] ) * par( 'usr' )[3] + 
	( Mapping_parameters[["legend_y"]][2] ) * par( 'usr' )[4]
align = c( "lt","rb" )[2]
gradient = c( "x", "y" )[2]
plotrix::color.legend( xl = xl, yb = yb, xr = xr, yt = yt, 
	legend = round( seq( zlim[1], zlim[2], length = 5 ), digits = 2 ),
	rect.col = sf.colors( n = 10, cutoff.tails = c( 0.35, 0.2 ), alpha = 1, categorical = FALSE ),
	cex = 1.2, align = align, gradient = gradient, col = "black" )
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_1, "_", 
	scenario, "_with a legend" ), path = MapsDir, type = "png" )
pause( 2 ); graphics.off()

#### Produce a version of the map without a legend and save it 
new.graph( 0.5 )
par( oma = c( 3, 3, 0, 0 ), cex = 1.1 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, main = "", axes = FALSE, 
	xlim = st_bbox(plot_sf)[c( 1, 3 )], ylim = st_bbox(plot_sf)[c( 2, 4 )] ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
plot( contour_200m, col = "white", lwd = 3.5, lty = 2, add = TRUE )
plot( contour_500m, col = "gray48", lwd = 3.5, lty = 2, add = TRUE )
box()
axis( side = 1, at = Mapping_parameters[["xaxis"]], labels = Mapping_parameters[["xaxis_label"]] )
axis( side = 2, at = Mapping_parameters[["yaxis"]], labels = Mapping_parameters[["yaxis_label"]] )
minor.tick( nx = 2, ny = 1 )
mtext( side = 1, outer = TRUE, "Longitude (°E)", cex = 1.25, line = par()$oma[1]/2, padj = 0 )
mtext( side = 2, outer = TRUE, "Latitude (°N)", cex = 1.25, line = par()$oma[2]/2, padj = 0, las = 0 )
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_1, "_", 
	scenario, "_without a legend" ), path = MapsDir, type = "png" )
pause( 2 ); graphics.off()

######## Produce a map of the approximated spatial correlation for the tinyVAST model using the SPDE method 
######## and covariate-based regional anisotropy under Scenario 1 for the second location of interest and save it
scenario <- "Regional_Anisotropy_Scenario1"
plot_sf <- r_Regional_Anisotropy_Scenario1_2
zlim <- range( c( r_Regional_Anisotropy_Scenario1_1$r1, r_Regional_Anisotropy_Scenario1_2$r1, 
	r_Regional_Anisotropy_Scenario1_3$r1, 
	r_Regional_Anisotropy_Scenario2_1$r1, r_Regional_Anisotropy_Scenario2_2$r1, 
	r_Regional_Anisotropy_Scenario2_3$r1, 
	r_Regional_Anisotropy_Scenario3_1$r1, r_Regional_Anisotropy_Scenario3_2$r1, 
	r_Regional_Anisotropy_Scenario3_3$r1, 
	r_Regional_Anisotropy_Scenario4_1$r1, r_Regional_Anisotropy_Scenario4_2$r1, 
	r_Regional_Anisotropy_Scenario4_3$r1 ), na.rm = TRUE ) 

#### Produce a version of the map with a legend and save it
new.graph( 0.5 )
par( oma = c( 3, 3, 0, 0 ), cex = 1.1 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, main = "", axes = FALSE, 
	xlim = st_bbox(plot_sf)[c( 1, 3 )], ylim = st_bbox(plot_sf)[c( 2, 4 )] ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
plot( contour_200m, col = "white", lwd = 3.5, lty = 2, add = TRUE )
plot( contour_500m, col = "gray48", lwd = 3.5, lty = 2, add = TRUE )
box()
axis( side = 1, at = Mapping_parameters[["xaxis"]], labels = Mapping_parameters[["xaxis_label"]] )
axis( side = 2, at = Mapping_parameters[["yaxis"]], labels = Mapping_parameters[["yaxis_label"]] )
minor.tick( nx = 2, ny = 1 )
mtext( side = 1, outer = TRUE, "Longitude (°E)", cex = 1.25, line = par()$oma[1]/2, padj = 0 )
mtext( side = 2, outer = TRUE, "Latitude (°N)", cex = 1.25, line = par()$oma[2]/2, padj = 0, las = 0 )
xl = ( 1 - Mapping_parameters[["legend_x"]][1] ) * par( 'usr' )[1] + 
	( Mapping_parameters[["legend_x"]][1] ) * par( 'usr' )[2]
xr = ( 1 - Mapping_parameters[["legend_x"]][2] ) * par( 'usr' )[1] + 
	( Mapping_parameters[["legend_x"]][2] ) * par( 'usr' )[2]
yb = ( 1 - Mapping_parameters[["legend_y"]][1] ) * par( 'usr' )[3] + 
	( Mapping_parameters[["legend_y"]][1] ) * par( 'usr' )[4]
yt = ( 1 - Mapping_parameters[["legend_y"]][2] ) * par( 'usr' )[3] + 
	( Mapping_parameters[["legend_y"]][2] ) * par( 'usr' )[4]
align = c( "lt","rb" )[2]
gradient = c( "x", "y" )[2]
plotrix::color.legend( xl = xl, yb = yb, xr = xr, yt = yt, 
	legend = round( seq( zlim[1], zlim[2], length = 5 ), digits = 2 ),
	rect.col = sf.colors( n = 10, cutoff.tails = c( 0.35, 0.2 ), alpha = 1, categorical = FALSE ),
	cex = 1.2, align = align, gradient = gradient, col = "black" )
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_2, "_", 
	scenario, "_with a legend" ), path = MapsDir, type = "png" )
pause( 2 ); graphics.off()

#### Produce a version of the map without a legend and save it 
new.graph( 0.5 )
par( oma = c( 3, 3, 0, 0 ), cex = 1.1 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, main = "", axes = FALSE, 
	xlim = st_bbox(plot_sf)[c( 1, 3 )], ylim = st_bbox(plot_sf)[c( 2, 4 )] ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
plot( contour_200m, col = "white", lwd = 3.5, lty = 2, add = TRUE )
plot( contour_500m, col = "gray48", lwd = 3.5, lty = 2, add = TRUE )
box()
axis( side = 1, at = Mapping_parameters[["xaxis"]], labels = Mapping_parameters[["xaxis_label"]] )
axis( side = 2, at = Mapping_parameters[["yaxis"]], labels = Mapping_parameters[["yaxis_label"]] )
minor.tick( nx = 2, ny = 1 )
mtext( side = 1, outer = TRUE, "Longitude (°E)", cex = 1.25, line = par()$oma[1]/2, padj = 0 )
mtext( side = 2, outer = TRUE, "Latitude (°N)", cex = 1.25, line = par()$oma[2]/2, padj = 0, las = 0 )
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_2, "_", 
	scenario, "_without a legend" ), path = MapsDir, type = "png" )
pause( 2 ); graphics.off()

######## Produce a map of the approximated spatial correlation for the tinyVAST model using the SPDE method 
######## and covariate-based regional anisotropy under Scenario 1 for the third location of interest and save it
scenario <- "Regional_Anisotropy_Scenario1"
plot_sf <- r_Regional_Anisotropy_Scenario1_3
zlim <- range( c( r_Regional_Anisotropy_Scenario1_1$r1, r_Regional_Anisotropy_Scenario1_2$r1, 
	r_Regional_Anisotropy_Scenario1_3$r1, 
	r_Regional_Anisotropy_Scenario2_1$r1, r_Regional_Anisotropy_Scenario2_2$r1, 
	r_Regional_Anisotropy_Scenario2_3$r1, 
	r_Regional_Anisotropy_Scenario3_1$r1, r_Regional_Anisotropy_Scenario3_2$r1, 
	r_Regional_Anisotropy_Scenario3_3$r1, 
	r_Regional_Anisotropy_Scenario4_1$r1, r_Regional_Anisotropy_Scenario4_2$r1, 
	r_Regional_Anisotropy_Scenario4_3$r1 ), na.rm = TRUE ) 

#### Produce a version of the map with a legend and save it
new.graph( 0.5 )
par( oma = c( 3, 3, 0, 0 ), cex = 1.1 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, main = "", axes = FALSE, 
	xlim = st_bbox(plot_sf)[c( 1, 3 )], ylim = st_bbox(plot_sf)[c( 2, 4 )] ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
plot( contour_200m, col = "white", lwd = 3.5, lty = 2, add = TRUE )
plot( contour_500m, col = "gray48", lwd = 3.5, lty = 2, add = TRUE )
box()
axis( side = 1, at = Mapping_parameters[["xaxis"]], labels = Mapping_parameters[["xaxis_label"]] )
axis( side = 2, at = Mapping_parameters[["yaxis"]], labels = Mapping_parameters[["yaxis_label"]] )
minor.tick( nx = 2, ny = 1 )
mtext( side = 1, outer = TRUE, "Longitude (°E)", cex = 1.25, line = par()$oma[1]/2, padj = 0 )
mtext( side = 2, outer = TRUE, "Latitude (°N)", cex = 1.25, line = par()$oma[2]/2, padj = 0, las = 0 )
xl = ( 1 - Mapping_parameters[["legend_x"]][1] ) * par( 'usr' )[1] + 
	( Mapping_parameters[["legend_x"]][1] ) * par( 'usr' )[2]
xr = ( 1 - Mapping_parameters[["legend_x"]][2] ) * par( 'usr' )[1] + 
	( Mapping_parameters[["legend_x"]][2] ) * par( 'usr' )[2]
yb = ( 1 - Mapping_parameters[["legend_y"]][1] ) * par( 'usr' )[3] + 
	( Mapping_parameters[["legend_y"]][1] ) * par( 'usr' )[4]
yt = ( 1 - Mapping_parameters[["legend_y"]][2] ) * par( 'usr' )[3] + 
	( Mapping_parameters[["legend_y"]][2] ) * par( 'usr' )[4]
align = c( "lt","rb" )[2]
gradient = c( "x", "y" )[2]
plotrix::color.legend( xl = xl, yb = yb, xr = xr, yt = yt, 
	legend = round( seq( zlim[1], zlim[2], length = 5 ), digits = 2 ),
	rect.col = sf.colors( n = 10, cutoff.tails = c( 0.35, 0.2 ), alpha = 1, categorical = FALSE ),
	cex = 1.2, align = align, gradient = gradient, col = "black" )
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_3, "_", 
	scenario, "_with a legend" ), path = MapsDir, type = "png" )
pause( 2 ); graphics.off()

#### Produce a version of the map without a legend and save it 
new.graph( 0.5 )
par( oma = c( 3, 3, 0, 0 ), cex = 1.1 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, main = "", axes = FALSE, 
	xlim = st_bbox(plot_sf)[c( 1, 3 )], ylim = st_bbox(plot_sf)[c( 2, 4 )] ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
plot( contour_200m, col = "white", lwd = 3.5, lty = 2, add = TRUE )
plot( contour_500m, col = "gray48", lwd = 3.5, lty = 2, add = TRUE )
box()
axis( side = 1, at = Mapping_parameters[["xaxis"]], labels = Mapping_parameters[["xaxis_label"]] )
axis( side = 2, at = Mapping_parameters[["yaxis"]], labels = Mapping_parameters[["yaxis_label"]] )
minor.tick( nx = 2, ny = 1 )
mtext( side = 1, outer = TRUE, "Longitude (°E)", cex = 1.25, line = par()$oma[1]/2, padj = 0 )
mtext( side = 2, outer = TRUE, "Latitude (°N)", cex = 1.25, line = par()$oma[2]/2, padj = 0, las = 0 )
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_3, "_", 
	scenario, "_without a legend" ), path = MapsDir, type = "png" )
pause( 2 ); graphics.off()

######## Produce a map of the approximated spatial correlation for the tinyVAST model using the SPDE method 
######## and covariate-based regional anisotropy under Scenario 2 for the first location of interest and save it
scenario <- "Regional_Anisotropy_Scenario2"
plot_sf <- r_Regional_Anisotropy_Scenario2_1
zlim <- range( c( r_Regional_Anisotropy_Scenario1_1$r1, r_Regional_Anisotropy_Scenario1_2$r1, 
	r_Regional_Anisotropy_Scenario1_3$r1, 
	r_Regional_Anisotropy_Scenario2_1$r1, r_Regional_Anisotropy_Scenario2_2$r1, 
	r_Regional_Anisotropy_Scenario2_3$r1, 
	r_Regional_Anisotropy_Scenario3_1$r1, r_Regional_Anisotropy_Scenario3_2$r1, 
	r_Regional_Anisotropy_Scenario3_3$r1, 
	r_Regional_Anisotropy_Scenario4_1$r1, r_Regional_Anisotropy_Scenario4_2$r1, 
	r_Regional_Anisotropy_Scenario4_3$r1 ), na.rm = TRUE ) 

#### Produce a version of the map with a legend and save it
new.graph( 0.5 )
par( oma = c( 3, 3, 0, 0 ), cex = 1.1 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, main = "", axes = FALSE, 
	xlim = st_bbox(plot_sf)[c( 1, 3 )], ylim = st_bbox(plot_sf)[c( 2, 4 )] ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
plot( contour_200m, col = "white", lwd = 3.5, lty = 2, add = TRUE )
plot( contour_500m, col = "gray48", lwd = 3.5, lty = 2, add = TRUE )
box()
axis( side = 1, at = Mapping_parameters[["xaxis"]], labels = Mapping_parameters[["xaxis_label"]] )
axis( side = 2, at = Mapping_parameters[["yaxis"]], labels = Mapping_parameters[["yaxis_label"]] )
minor.tick( nx = 2, ny = 1 )
mtext( side = 1, outer = TRUE, "Longitude (°E)", cex = 1.25, line = par()$oma[1]/2, padj = 0 )
mtext( side = 2, outer = TRUE, "Latitude (°N)", cex = 1.25, line = par()$oma[2]/2, padj = 0, las = 0 )
xl = ( 1 - Mapping_parameters[["legend_x"]][1] ) * par( 'usr' )[1] + 
	( Mapping_parameters[["legend_x"]][1] ) * par( 'usr' )[2]
xr = ( 1 - Mapping_parameters[["legend_x"]][2] ) * par( 'usr' )[1] + 
	( Mapping_parameters[["legend_x"]][2] ) * par( 'usr' )[2]
yb = ( 1 - Mapping_parameters[["legend_y"]][1] ) * par( 'usr' )[3] + 
	( Mapping_parameters[["legend_y"]][1] ) * par( 'usr' )[4]
yt = ( 1 - Mapping_parameters[["legend_y"]][2] ) * par( 'usr' )[3] + 
	( Mapping_parameters[["legend_y"]][2] ) * par( 'usr' )[4]
align = c( "lt","rb" )[2]
gradient = c( "x", "y" )[2]
plotrix::color.legend( xl = xl, yb = yb, xr = xr, yt = yt, 
	legend = round( seq( zlim[1], zlim[2], length = 5 ), digits = 2 ),
	rect.col = sf.colors( n = 10, cutoff.tails = c( 0.35, 0.2 ), alpha = 1, categorical = FALSE ),
	cex = 1.2, align = align, gradient = gradient, col = "black" )
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_1, "_", 
	scenario, "_with a legend" ), path = MapsDir, type = "png" )
pause( 2 ); graphics.off()

#### Produce a version of the map without a legend and save it 
new.graph( 0.5 )
par( oma = c( 3, 3, 0, 0 ), cex = 1.1 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, main = "", axes = FALSE, 
	xlim = st_bbox(plot_sf)[c( 1, 3 )], ylim = st_bbox(plot_sf)[c( 2, 4 )] ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
plot( contour_200m, col = "white", lwd = 3.5, lty = 2, add = TRUE )
plot( contour_500m, col = "gray48", lwd = 3.5, lty = 2, add = TRUE )
box()
axis( side = 1, at = Mapping_parameters[["xaxis"]], labels = Mapping_parameters[["xaxis_label"]] )
axis( side = 2, at = Mapping_parameters[["yaxis"]], labels = Mapping_parameters[["yaxis_label"]] )
minor.tick( nx = 2, ny = 1 )
mtext( side = 1, outer = TRUE, "Longitude (°E)", cex = 1.25, line = par()$oma[1]/2, padj = 0 )
mtext( side = 2, outer = TRUE, "Latitude (°N)", cex = 1.25, line = par()$oma[2]/2, padj = 0, las = 0 )
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_1, "_", 
	scenario, "_without a legend" ), path = MapsDir, type = "png" )
pause( 2 ); graphics.off()

######## Produce a map of the approximated spatial correlation for the tinyVAST model using the SPDE method 
######## and covariate-based regional anisotropy under Scenario 2 for the second location of interest and save it
scenario <- "Regional_Anisotropy_Scenario2"
plot_sf <- r_Regional_Anisotropy_Scenario2_2
zlim <- range( c( r_Regional_Anisotropy_Scenario1_1$r1, r_Regional_Anisotropy_Scenario1_2$r1, 
	r_Regional_Anisotropy_Scenario1_3$r1, 
	r_Regional_Anisotropy_Scenario2_1$r1, r_Regional_Anisotropy_Scenario2_2$r1, 
	r_Regional_Anisotropy_Scenario2_3$r1, 
	r_Regional_Anisotropy_Scenario3_1$r1, r_Regional_Anisotropy_Scenario3_2$r1, 
	r_Regional_Anisotropy_Scenario3_3$r1, 
	r_Regional_Anisotropy_Scenario4_1$r1, r_Regional_Anisotropy_Scenario4_2$r1, 
	r_Regional_Anisotropy_Scenario4_3$r1 ), na.rm = TRUE ) 

#### Produce a version of the map with a legend and save it
new.graph( 0.5 )
par( oma = c( 3, 3, 0, 0 ), cex = 1.1 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, main = "", axes = FALSE, 
	xlim = st_bbox(plot_sf)[c( 1, 3 )], ylim = st_bbox(plot_sf)[c( 2, 4 )] ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
plot( contour_200m, col = "white", lwd = 3.5, lty = 2, add = TRUE )
plot( contour_500m, col = "gray48", lwd = 3.5, lty = 2, add = TRUE )
box()
axis( side = 1, at = Mapping_parameters[["xaxis"]], labels = Mapping_parameters[["xaxis_label"]] )
axis( side = 2, at = Mapping_parameters[["yaxis"]], labels = Mapping_parameters[["yaxis_label"]] )
minor.tick( nx = 2, ny = 1 )
mtext( side = 1, outer = TRUE, "Longitude (°E)", cex = 1.25, line = par()$oma[1]/2, padj = 0 )
mtext( side = 2, outer = TRUE, "Latitude (°N)", cex = 1.25, line = par()$oma[2]/2, padj = 0, las = 0 )
xl = ( 1 - Mapping_parameters[["legend_x"]][1] ) * par( 'usr' )[1] + 
	( Mapping_parameters[["legend_x"]][1] ) * par( 'usr' )[2]
xr = ( 1 - Mapping_parameters[["legend_x"]][2] ) * par( 'usr' )[1] + 
	( Mapping_parameters[["legend_x"]][2] ) * par( 'usr' )[2]
yb = ( 1 - Mapping_parameters[["legend_y"]][1] ) * par( 'usr' )[3] + 
	( Mapping_parameters[["legend_y"]][1] ) * par( 'usr' )[4]
yt = ( 1 - Mapping_parameters[["legend_y"]][2] ) * par( 'usr' )[3] + 
	( Mapping_parameters[["legend_y"]][2] ) * par( 'usr' )[4]
align = c( "lt","rb" )[2]
gradient = c( "x", "y" )[2]
plotrix::color.legend( xl = xl, yb = yb, xr = xr, yt = yt, 
	legend = round( seq( zlim[1], zlim[2], length = 5 ), digits = 2 ),
	rect.col = sf.colors( n = 10, cutoff.tails = c( 0.35, 0.2 ), alpha = 1, categorical = FALSE ),
	cex = 1.2, align = align, gradient = gradient, col = "black" )
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_2, "_", 
	scenario, "_with a legend" ), path = MapsDir, type = "png" )
pause( 2 ); graphics.off()

#### Produce a version of the map without a legend and save it 
new.graph( 0.5 )
par( oma = c( 3, 3, 0, 0 ), cex = 1.1 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, main = "", axes = FALSE, 
	xlim = st_bbox(plot_sf)[c( 1, 3 )], ylim = st_bbox(plot_sf)[c( 2, 4 )] ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
plot( contour_200m, col = "white", lwd = 3.5, lty = 2, add = TRUE )
plot( contour_500m, col = "gray48", lwd = 3.5, lty = 2, add = TRUE )
box()
axis( side = 1, at = Mapping_parameters[["xaxis"]], labels = Mapping_parameters[["xaxis_label"]] )
axis( side = 2, at = Mapping_parameters[["yaxis"]], labels = Mapping_parameters[["yaxis_label"]] )
minor.tick( nx = 2, ny = 1 )
mtext( side = 1, outer = TRUE, "Longitude (°E)", cex = 1.25, line = par()$oma[1]/2, padj = 0 )
mtext( side = 2, outer = TRUE, "Latitude (°N)", cex = 1.25, line = par()$oma[2]/2, padj = 0, las = 0 )
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_2, "_", 
	scenario, "_without a legend" ), path = MapsDir, type = "png" )
pause( 2 ); graphics.off()

######## Produce a map of the approximated spatial correlation for the tinyVAST model using the SPDE method 
######## and covariate-based regional anisotropy under Scenario 2 for the third location of interest and save it
scenario <- "Regional_Anisotropy_Scenario2"
plot_sf <- r_Regional_Anisotropy_Scenario2_3
zlim <- range( c( r_Regional_Anisotropy_Scenario1_1$r1, r_Regional_Anisotropy_Scenario1_2$r1, 
	r_Regional_Anisotropy_Scenario1_3$r1, 
	r_Regional_Anisotropy_Scenario2_1$r1, r_Regional_Anisotropy_Scenario2_2$r1, 
	r_Regional_Anisotropy_Scenario2_3$r1, 
	r_Regional_Anisotropy_Scenario3_1$r1, r_Regional_Anisotropy_Scenario3_2$r1, 
	r_Regional_Anisotropy_Scenario3_3$r1, 
	r_Regional_Anisotropy_Scenario4_1$r1, r_Regional_Anisotropy_Scenario4_2$r1, 
	r_Regional_Anisotropy_Scenario4_3$r1 ), na.rm = TRUE ) 

#### Produce a version of the map with a legend and save it
new.graph( 0.5 )
par( oma = c( 3, 3, 0, 0 ), cex = 1.1 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, main = "", axes = FALSE, 
	xlim = st_bbox(plot_sf)[c( 1, 3 )], ylim = st_bbox(plot_sf)[c( 2, 4 )] ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
plot( contour_200m, col = "white", lwd = 3.5, lty = 2, add = TRUE )
plot( contour_500m, col = "gray48", lwd = 3.5, lty = 2, add = TRUE )
box()
axis( side = 1, at = Mapping_parameters[["xaxis"]], labels = Mapping_parameters[["xaxis_label"]] )
axis( side = 2, at = Mapping_parameters[["yaxis"]], labels = Mapping_parameters[["yaxis_label"]] )
minor.tick( nx = 2, ny = 1 )
mtext( side = 1, outer = TRUE, "Longitude (°E)", cex = 1.25, line = par()$oma[1]/2, padj = 0 )
mtext( side = 2, outer = TRUE, "Latitude (°N)", cex = 1.25, line = par()$oma[2]/2, padj = 0, las = 0 )
xl = ( 1 - Mapping_parameters[["legend_x"]][1] ) * par( 'usr' )[1] + 
	( Mapping_parameters[["legend_x"]][1] ) * par( 'usr' )[2]
xr = ( 1 - Mapping_parameters[["legend_x"]][2] ) * par( 'usr' )[1] + 
	( Mapping_parameters[["legend_x"]][2] ) * par( 'usr' )[2]
yb = ( 1 - Mapping_parameters[["legend_y"]][1] ) * par( 'usr' )[3] + 
	( Mapping_parameters[["legend_y"]][1] ) * par( 'usr' )[4]
yt = ( 1 - Mapping_parameters[["legend_y"]][2] ) * par( 'usr' )[3] + 
	( Mapping_parameters[["legend_y"]][2] ) * par( 'usr' )[4]
align = c( "lt","rb" )[2]
gradient = c( "x", "y" )[2]
plotrix::color.legend( xl = xl, yb = yb, xr = xr, yt = yt, 
	legend = round( seq( zlim[1], zlim[2], length = 5 ), digits = 2 ),
	rect.col = sf.colors( n = 10, cutoff.tails = c( 0.35, 0.2 ), alpha = 1, categorical = FALSE ),
	cex = 1.2, align = align, gradient = gradient, col = "black" )
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_3, "_", 
	scenario, "_with a legend" ), path = MapsDir, type = "png" )
pause( 2 ); graphics.off()

#### Produce a version of the map without a legend and save it 
new.graph( 0.5 )
par( oma = c( 3, 3, 0, 0 ), cex = 1.1 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, main = "", axes = FALSE, 
	xlim = st_bbox(plot_sf)[c( 1, 3 )], ylim = st_bbox(plot_sf)[c( 2, 4 )] ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
plot( contour_200m, col = "white", lwd = 3.5, lty = 2, add = TRUE )
plot( contour_500m, col = "gray48", lwd = 3.5, lty = 2, add = TRUE )
box()
axis( side = 1, at = Mapping_parameters[["xaxis"]], labels = Mapping_parameters[["xaxis_label"]] )
axis( side = 2, at = Mapping_parameters[["yaxis"]], labels = Mapping_parameters[["yaxis_label"]] )
minor.tick( nx = 2, ny = 1 )
mtext( side = 1, outer = TRUE, "Longitude (°E)", cex = 1.25, line = par()$oma[1]/2, padj = 0 )
mtext( side = 2, outer = TRUE, "Latitude (°N)", cex = 1.25, line = par()$oma[2]/2, padj = 0, las = 0 )
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_3, "_", 
	scenario, "_without a legend" ), path = MapsDir, type = "png" )
pause( 2 ); graphics.off()

######## Produce a map of the approximated spatial correlation for the tinyVAST model using the SPDE method 
######## and covariate-based regional anisotropy under Scenario 3 for the first location of interest and save it
scenario <- "Regional_Anisotropy_Scenario3"
plot_sf <- r_Regional_Anisotropy_Scenario3_1
zlim <- range( c( r_Regional_Anisotropy_Scenario1_1$r1, r_Regional_Anisotropy_Scenario1_2$r1, 
	r_Regional_Anisotropy_Scenario1_3$r1, 
	r_Regional_Anisotropy_Scenario2_1$r1, r_Regional_Anisotropy_Scenario2_2$r1, 
	r_Regional_Anisotropy_Scenario2_3$r1, 
	r_Regional_Anisotropy_Scenario3_1$r1, r_Regional_Anisotropy_Scenario3_2$r1, 
	r_Regional_Anisotropy_Scenario3_3$r1, 
	r_Regional_Anisotropy_Scenario4_1$r1, r_Regional_Anisotropy_Scenario4_2$r1, 
	r_Regional_Anisotropy_Scenario4_3$r1 ), na.rm = TRUE ) 

#### Produce a version of the map with a legend and save it
new.graph( 0.5 )
par( oma = c( 3, 3, 0, 0 ), cex = 1.1 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, main = "", axes = FALSE, 
	xlim = st_bbox(plot_sf)[c( 1, 3 )], ylim = st_bbox(plot_sf)[c( 2, 4 )] ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
plot( contour_200m, col = "white", lwd = 3.5, lty = 2, add = TRUE )
plot( contour_500m, col = "gray48", lwd = 3.5, lty = 2, add = TRUE )
box()
axis( side = 1, at = Mapping_parameters[["xaxis"]], labels = Mapping_parameters[["xaxis_label"]] )
axis( side = 2, at = Mapping_parameters[["yaxis"]], labels = Mapping_parameters[["yaxis_label"]] )
minor.tick( nx = 2, ny = 1 )
mtext( side = 1, outer = TRUE, "Longitude (°E)", cex = 1.25, line = par()$oma[1]/2, padj = 0 )
mtext( side = 2, outer = TRUE, "Latitude (°N)", cex = 1.25, line = par()$oma[2]/2, padj = 0, las = 0 )
xl = ( 1 - Mapping_parameters[["legend_x"]][1] ) * par( 'usr' )[1] + 
	( Mapping_parameters[["legend_x"]][1] ) * par( 'usr' )[2]
xr = ( 1 - Mapping_parameters[["legend_x"]][2] ) * par( 'usr' )[1] + 
	( Mapping_parameters[["legend_x"]][2] ) * par( 'usr' )[2]
yb = ( 1 - Mapping_parameters[["legend_y"]][1] ) * par( 'usr' )[3] + 
	( Mapping_parameters[["legend_y"]][1] ) * par( 'usr' )[4]
yt = ( 1 - Mapping_parameters[["legend_y"]][2] ) * par( 'usr' )[3] + 
	( Mapping_parameters[["legend_y"]][2] ) * par( 'usr' )[4]
align = c( "lt","rb" )[2]
gradient = c( "x", "y" )[2]
plotrix::color.legend( xl = xl, yb = yb, xr = xr, yt = yt, 
	legend = round( seq( zlim[1], zlim[2], length = 5 ), digits = 2 ),
	rect.col = sf.colors( n = 10, cutoff.tails = c( 0.35, 0.2 ), alpha = 1, categorical = FALSE ),
	cex = 1.2, align = align, gradient = gradient, col = "black" )
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_1, "_", 
	scenario, "_with a legend" ), path = MapsDir, type = "png" )
pause( 2 ); graphics.off()

#### Produce a version of the map without a legend and save it 
new.graph( 0.5 )
par( oma = c( 3, 3, 0, 0 ), cex = 1.1 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, main = "", axes = FALSE, 
	xlim = st_bbox(plot_sf)[c( 1, 3 )], ylim = st_bbox(plot_sf)[c( 2, 4 )] ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
plot( contour_200m, col = "white", lwd = 3.5, lty = 2, add = TRUE )
plot( contour_500m, col = "gray48", lwd = 3.5, lty = 2, add = TRUE )
box()
axis( side = 1, at = Mapping_parameters[["xaxis"]], labels = Mapping_parameters[["xaxis_label"]] )
axis( side = 2, at = Mapping_parameters[["yaxis"]], labels = Mapping_parameters[["yaxis_label"]] )
minor.tick( nx = 2, ny = 1 )
mtext( side = 1, outer = TRUE, "Longitude (°E)", cex = 1.25, line = par()$oma[1]/2, padj = 0 )
mtext( side = 2, outer = TRUE, "Latitude (°N)", cex = 1.25, line = par()$oma[2]/2, padj = 0, las = 0 )
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_1, "_", 
	scenario, "_without a legend" ), path = MapsDir, type = "png" )
pause( 2 ); graphics.off()

######## Produce a map of the approximated spatial correlation for the tinyVAST model using the SPDE method 
######## and covariate-based regional anisotropy under Scenario 3 for the second location of interest and save it
scenario <- "Regional_Anisotropy_Scenario3"
plot_sf <- r_Regional_Anisotropy_Scenario3_2
zlim <- range( c( r_Regional_Anisotropy_Scenario1_1$r1, r_Regional_Anisotropy_Scenario1_2$r1, 
	r_Regional_Anisotropy_Scenario1_3$r1, 
	r_Regional_Anisotropy_Scenario2_1$r1, r_Regional_Anisotropy_Scenario2_2$r1, 
	r_Regional_Anisotropy_Scenario2_3$r1, 
	r_Regional_Anisotropy_Scenario3_1$r1, r_Regional_Anisotropy_Scenario3_2$r1, 
	r_Regional_Anisotropy_Scenario3_3$r1, 
	r_Regional_Anisotropy_Scenario4_1$r1, r_Regional_Anisotropy_Scenario4_2$r1, 
	r_Regional_Anisotropy_Scenario4_3$r1 ), na.rm = TRUE ) 

#### Produce a version of the map with a legend and save it
new.graph( 0.5 )
par( oma = c( 3, 3, 0, 0 ), cex = 1.1 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, main = "", axes = FALSE, 
	xlim = st_bbox(plot_sf)[c( 1, 3 )], ylim = st_bbox(plot_sf)[c( 2, 4 )] ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
plot( contour_200m, col = "white", lwd = 3.5, lty = 2, add = TRUE )
plot( contour_500m, col = "gray48", lwd = 3.5, lty = 2, add = TRUE )
box()
axis( side = 1, at = Mapping_parameters[["xaxis"]], labels = Mapping_parameters[["xaxis_label"]] )
axis( side = 2, at = Mapping_parameters[["yaxis"]], labels = Mapping_parameters[["yaxis_label"]] )
minor.tick( nx = 2, ny = 1 )
mtext( side = 1, outer = TRUE, "Longitude (°E)", cex = 1.25, line = par()$oma[1]/2, padj = 0 )
mtext( side = 2, outer = TRUE, "Latitude (°N)", cex = 1.25, line = par()$oma[2]/2, padj = 0, las = 0 )
xl = ( 1 - Mapping_parameters[["legend_x"]][1] ) * par( 'usr' )[1] + 
	( Mapping_parameters[["legend_x"]][1] ) * par( 'usr' )[2]
xr = ( 1 - Mapping_parameters[["legend_x"]][2] ) * par( 'usr' )[1] + 
	( Mapping_parameters[["legend_x"]][2] ) * par( 'usr' )[2]
yb = ( 1 - Mapping_parameters[["legend_y"]][1] ) * par( 'usr' )[3] + 
	( Mapping_parameters[["legend_y"]][1] ) * par( 'usr' )[4]
yt = ( 1 - Mapping_parameters[["legend_y"]][2] ) * par( 'usr' )[3] + 
	( Mapping_parameters[["legend_y"]][2] ) * par( 'usr' )[4]
align = c( "lt","rb" )[2]
gradient = c( "x", "y" )[2]
plotrix::color.legend( xl = xl, yb = yb, xr = xr, yt = yt, 
	legend = round( seq( zlim[1], zlim[2], length = 5 ), digits = 2 ),
	rect.col = sf.colors( n = 10, cutoff.tails = c( 0.35, 0.2 ), alpha = 1, categorical = FALSE ),
	cex = 1.2, align = align, gradient = gradient, col = "black" )
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_2, "_", 
	scenario, "_with a legend" ), path = MapsDir, type = "png" )
pause( 2 ); graphics.off()

#### Produce a version of the map without a legend and save it 
new.graph( 0.5 )
par( oma = c( 3, 3, 0, 0 ), cex = 1.1 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, main = "", axes = FALSE, 
	xlim = st_bbox(plot_sf)[c( 1, 3 )], ylim = st_bbox(plot_sf)[c( 2, 4 )] ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
plot( contour_200m, col = "white", lwd = 3.5, lty = 2, add = TRUE )
plot( contour_500m, col = "gray48", lwd = 3.5, lty = 2, add = TRUE )
box()
axis( side = 1, at = Mapping_parameters[["xaxis"]], labels = Mapping_parameters[["xaxis_label"]] )
axis( side = 2, at = Mapping_parameters[["yaxis"]], labels = Mapping_parameters[["yaxis_label"]] )
minor.tick( nx = 2, ny = 1 )
mtext( side = 1, outer = TRUE, "Longitude (°E)", cex = 1.25, line = par()$oma[1]/2, padj = 0 )
mtext( side = 2, outer = TRUE, "Latitude (°N)", cex = 1.25, line = par()$oma[2]/2, padj = 0, las = 0 )
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_2, "_", 
	scenario, "_without a legend" ), path = MapsDir, type = "png" )
pause( 2 ); graphics.off()

######## Produce a map of the approximated spatial correlation for the tinyVAST model using the SPDE method 
######## and covariate-based regional anisotropy under Scenario 3 for the third location of interest and save it
scenario <- "Regional_Anisotropy_Scenario3"
plot_sf <- r_Regional_Anisotropy_Scenario3_3
zlim <- range( c( r_Regional_Anisotropy_Scenario1_1$r1, r_Regional_Anisotropy_Scenario1_2$r1, 
	r_Regional_Anisotropy_Scenario1_3$r1, 
	r_Regional_Anisotropy_Scenario2_1$r1, r_Regional_Anisotropy_Scenario2_2$r1, 
	r_Regional_Anisotropy_Scenario2_3$r1, 
	r_Regional_Anisotropy_Scenario3_1$r1, r_Regional_Anisotropy_Scenario3_2$r1, 
	r_Regional_Anisotropy_Scenario3_3$r1, 
	r_Regional_Anisotropy_Scenario4_1$r1, r_Regional_Anisotropy_Scenario4_2$r1, 
	r_Regional_Anisotropy_Scenario4_3$r1 ), na.rm = TRUE ) 

#### Produce a version of the map with a legend and save it
new.graph( 0.5 )
par( oma = c( 3, 3, 0, 0 ), cex = 1.1 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, main = "", axes = FALSE, 
	xlim = st_bbox(plot_sf)[c( 1, 3 )], ylim = st_bbox(plot_sf)[c( 2, 4 )] ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
plot( contour_200m, col = "white", lwd = 3.5, lty = 2, add = TRUE )
plot( contour_500m, col = "gray48", lwd = 3.5, lty = 2, add = TRUE )
box()
axis( side = 1, at = Mapping_parameters[["xaxis"]], labels = Mapping_parameters[["xaxis_label"]] )
axis( side = 2, at = Mapping_parameters[["yaxis"]], labels = Mapping_parameters[["yaxis_label"]] )
minor.tick( nx = 2, ny = 1 )
mtext( side = 1, outer = TRUE, "Longitude (°E)", cex = 1.25, line = par()$oma[1]/2, padj = 0 )
mtext( side = 2, outer = TRUE, "Latitude (°N)", cex = 1.25, line = par()$oma[2]/2, padj = 0, las = 0 )
xl = ( 1 - Mapping_parameters[["legend_x"]][1] ) * par( 'usr' )[1] + 
	( Mapping_parameters[["legend_x"]][1] ) * par( 'usr' )[2]
xr = ( 1 - Mapping_parameters[["legend_x"]][2] ) * par( 'usr' )[1] + 
	( Mapping_parameters[["legend_x"]][2] ) * par( 'usr' )[2]
yb = ( 1 - Mapping_parameters[["legend_y"]][1] ) * par( 'usr' )[3] + 
	( Mapping_parameters[["legend_y"]][1] ) * par( 'usr' )[4]
yt = ( 1 - Mapping_parameters[["legend_y"]][2] ) * par( 'usr' )[3] + 
	( Mapping_parameters[["legend_y"]][2] ) * par( 'usr' )[4]
align = c( "lt","rb" )[2]
gradient = c( "x", "y" )[2]
plotrix::color.legend( xl = xl, yb = yb, xr = xr, yt = yt, 
	legend = round( seq( zlim[1], zlim[2], length = 5 ), digits = 2 ),
	rect.col = sf.colors( n = 10, cutoff.tails = c( 0.35, 0.2 ), alpha = 1, categorical = FALSE ),
	cex = 1.2, align = align, gradient = gradient, col = "black" )
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_3, "_", 
	scenario, "_with a legend" ), path = MapsDir, type = "png" )
pause( 2 ); graphics.off()

#### Produce a version of the map without a legend and save it 
new.graph( 0.5 )
par( oma = c( 3, 3, 0, 0 ), cex = 1.1 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, main = "", axes = FALSE, 
	xlim = st_bbox(plot_sf)[c( 1, 3 )], ylim = st_bbox(plot_sf)[c( 2, 4 )] ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
plot( contour_200m, col = "white", lwd = 3.5, lty = 2, add = TRUE )
plot( contour_500m, col = "gray48", lwd = 3.5, lty = 2, add = TRUE )
box()
axis( side = 1, at = Mapping_parameters[["xaxis"]], labels = Mapping_parameters[["xaxis_label"]] )
axis( side = 2, at = Mapping_parameters[["yaxis"]], labels = Mapping_parameters[["yaxis_label"]] )
minor.tick( nx = 2, ny = 1 )
mtext( side = 1, outer = TRUE, "Longitude (°E)", cex = 1.25, line = par()$oma[1]/2, padj = 0 )
mtext( side = 2, outer = TRUE, "Latitude (°N)", cex = 1.25, line = par()$oma[2]/2, padj = 0, las = 0 )
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_3, "_", 
	scenario, "_without a legend" ), path = MapsDir, type = "png" )
pause( 2 ); graphics.off()

######## Produce a map of the approximated spatial correlation for the tinyVAST model using the SPDE method 
######## and covariate-based regional anisotropy under Scenario 4 for the first location of interest and save it
scenario <- "Regional_Anisotropy_Scenario4"
plot_sf <- r_Regional_Anisotropy_Scenario4_1
zlim <- range( c( r_Regional_Anisotropy_Scenario1_1$r1, r_Regional_Anisotropy_Scenario1_2$r1, 
	r_Regional_Anisotropy_Scenario1_3$r1, 
	r_Regional_Anisotropy_Scenario2_1$r1, r_Regional_Anisotropy_Scenario2_2$r1, 
	r_Regional_Anisotropy_Scenario2_3$r1, 
	r_Regional_Anisotropy_Scenario3_1$r1, r_Regional_Anisotropy_Scenario3_2$r1, 
	r_Regional_Anisotropy_Scenario3_3$r1, 
	r_Regional_Anisotropy_Scenario4_1$r1, r_Regional_Anisotropy_Scenario4_2$r1, 
	r_Regional_Anisotropy_Scenario4_3$r1 ), na.rm = TRUE ) 

#### Produce a version of the map with a legend and save it
new.graph( 0.5 )
par( oma = c( 3, 3, 0, 0 ), cex = 1.1 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, main = "", axes = FALSE, 
	xlim = st_bbox(plot_sf)[c( 1, 3 )], ylim = st_bbox(plot_sf)[c( 2, 4 )] ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
plot( contour_200m, col = "white", lwd = 3.5, lty = 2, add = TRUE )
plot( contour_500m, col = "gray48", lwd = 3.5, lty = 2, add = TRUE )
box()
axis( side = 1, at = Mapping_parameters[["xaxis"]], labels = Mapping_parameters[["xaxis_label"]] )
axis( side = 2, at = Mapping_parameters[["yaxis"]], labels = Mapping_parameters[["yaxis_label"]] )
minor.tick( nx = 2, ny = 1 )
mtext( side = 1, outer = TRUE, "Longitude (°E)", cex = 1.25, line = par()$oma[1]/2, padj = 0 )
mtext( side = 2, outer = TRUE, "Latitude (°N)", cex = 1.25, line = par()$oma[2]/2, padj = 0, las = 0 )
xl = ( 1 - Mapping_parameters[["legend_x"]][1] ) * par( 'usr' )[1] + 
	( Mapping_parameters[["legend_x"]][1] ) * par( 'usr' )[2]
xr = ( 1 - Mapping_parameters[["legend_x"]][2] ) * par( 'usr' )[1] + 
	( Mapping_parameters[["legend_x"]][2] ) * par( 'usr' )[2]
yb = ( 1 - Mapping_parameters[["legend_y"]][1] ) * par( 'usr' )[3] + 
	( Mapping_parameters[["legend_y"]][1] ) * par( 'usr' )[4]
yt = ( 1 - Mapping_parameters[["legend_y"]][2] ) * par( 'usr' )[3] + 
	( Mapping_parameters[["legend_y"]][2] ) * par( 'usr' )[4]
align = c( "lt","rb" )[2]
gradient = c( "x", "y" )[2]
plotrix::color.legend( xl = xl, yb = yb, xr = xr, yt = yt, 
	legend = round( seq( zlim[1], zlim[2], length = 5 ), digits = 2 ),
	rect.col = sf.colors( n = 10, cutoff.tails = c( 0.35, 0.2 ), alpha = 1, categorical = FALSE ),
	cex = 1.2, align = align, gradient = gradient, col = "black" )
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_1, "_", 
	scenario, "_with a legend" ), path = MapsDir, type = "png" )
pause( 2 ); graphics.off()

#### Produce a version of the map without a legend and save it 
new.graph( 0.5 )
par( oma = c( 3, 3, 0, 0 ), cex = 1.1 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, main = "", axes = FALSE, 
	xlim = st_bbox(plot_sf)[c( 1, 3 )], ylim = st_bbox(plot_sf)[c( 2, 4 )] ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
plot( contour_200m, col = "white", lwd = 3.5, lty = 2, add = TRUE )
plot( contour_500m, col = "gray48", lwd = 3.5, lty = 2, add = TRUE )
box()
axis( side = 1, at = Mapping_parameters[["xaxis"]], labels = Mapping_parameters[["xaxis_label"]] )
axis( side = 2, at = Mapping_parameters[["yaxis"]], labels = Mapping_parameters[["yaxis_label"]] )
minor.tick( nx = 2, ny = 1 )
mtext( side = 1, outer = TRUE, "Longitude (°E)", cex = 1.25, line = par()$oma[1]/2, padj = 0 )
mtext( side = 2, outer = TRUE, "Latitude (°N)", cex = 1.25, line = par()$oma[2]/2, padj = 0, las = 0 )
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_1, "_", 
	scenario, "_without a legend" ), path = MapsDir, type = "png" )
pause( 2 ); graphics.off()

######## Produce a map of the approximated spatial correlation for the tinyVAST model using the SPDE method 
######## and covariate-based regional anisotropy under Scenario 4 for the second location of interest and save it
scenario <- "Regional_Anisotropy_Scenario4"
plot_sf <- r_Regional_Anisotropy_Scenario4_2
zlim <- range( c( r_Regional_Anisotropy_Scenario1_1$r1, r_Regional_Anisotropy_Scenario1_2$r1, 
	r_Regional_Anisotropy_Scenario1_3$r1, 
	r_Regional_Anisotropy_Scenario2_1$r1, r_Regional_Anisotropy_Scenario2_2$r1, 
	r_Regional_Anisotropy_Scenario2_3$r1, 
	r_Regional_Anisotropy_Scenario3_1$r1, r_Regional_Anisotropy_Scenario3_2$r1, 
	r_Regional_Anisotropy_Scenario3_3$r1, 
	r_Regional_Anisotropy_Scenario4_1$r1, r_Regional_Anisotropy_Scenario4_2$r1, 
	r_Regional_Anisotropy_Scenario4_3$r1 ), na.rm = TRUE ) 

#### Produce a version of the map with a legend and save it
new.graph( 0.5 )
par( oma = c( 3, 3, 0, 0 ), cex = 1.1 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, main = "", axes = FALSE, 
	xlim = st_bbox(plot_sf)[c( 1, 3 )], ylim = st_bbox(plot_sf)[c( 2, 4 )] ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
plot( contour_200m, col = "white", lwd = 3.5, lty = 2, add = TRUE )
plot( contour_500m, col = "gray48", lwd = 3.5, lty = 2, add = TRUE )
box()
axis( side = 1, at = Mapping_parameters[["xaxis"]], labels = Mapping_parameters[["xaxis_label"]] )
axis( side = 2, at = Mapping_parameters[["yaxis"]], labels = Mapping_parameters[["yaxis_label"]] )
minor.tick( nx = 2, ny = 1 )
mtext( side = 1, outer = TRUE, "Longitude (°E)", cex = 1.25, line = par()$oma[1]/2, padj = 0 )
mtext( side = 2, outer = TRUE, "Latitude (°N)", cex = 1.25, line = par()$oma[2]/2, padj = 0, las = 0 )
xl = ( 1 - Mapping_parameters[["legend_x"]][1] ) * par( 'usr' )[1] + 
	( Mapping_parameters[["legend_x"]][1] ) * par( 'usr' )[2]
xr = ( 1 - Mapping_parameters[["legend_x"]][2] ) * par( 'usr' )[1] + 
	( Mapping_parameters[["legend_x"]][2] ) * par( 'usr' )[2]
yb = ( 1 - Mapping_parameters[["legend_y"]][1] ) * par( 'usr' )[3] + 
	( Mapping_parameters[["legend_y"]][1] ) * par( 'usr' )[4]
yt = ( 1 - Mapping_parameters[["legend_y"]][2] ) * par( 'usr' )[3] + 
	( Mapping_parameters[["legend_y"]][2] ) * par( 'usr' )[4]
align = c( "lt","rb" )[2]
gradient = c( "x", "y" )[2]
plotrix::color.legend( xl = xl, yb = yb, xr = xr, yt = yt, 
	legend = round( seq( zlim[1], zlim[2], length = 5 ), digits = 2 ),
	rect.col = sf.colors( n = 10, cutoff.tails = c( 0.35, 0.2 ), alpha = 1, categorical = FALSE ),
	cex = 1.2, align = align, gradient = gradient, col = "black" )
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_2, "_", 
	scenario, "_with a legend" ), path = MapsDir, type = "png" )
pause( 2 ); graphics.off()

#### Produce a version of the map without a legend and save it 
new.graph( 0.5 )
par( oma = c( 3, 3, 0, 0 ), cex = 1.1 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, main = "", axes = FALSE, 
	xlim = st_bbox(plot_sf)[c( 1, 3 )], ylim = st_bbox(plot_sf)[c( 2, 4 )] ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
plot( contour_200m, col = "white", lwd = 3.5, lty = 2, add = TRUE )
plot( contour_500m, col = "gray48", lwd = 3.5, lty = 2, add = TRUE )
box()
axis( side = 1, at = Mapping_parameters[["xaxis"]], labels = Mapping_parameters[["xaxis_label"]] )
axis( side = 2, at = Mapping_parameters[["yaxis"]], labels = Mapping_parameters[["yaxis_label"]] )
minor.tick( nx = 2, ny = 1 )
mtext( side = 1, outer = TRUE, "Longitude (°E)", cex = 1.25, line = par()$oma[1]/2, padj = 0 )
mtext( side = 2, outer = TRUE, "Latitude (°N)", cex = 1.25, line = par()$oma[2]/2, padj = 0, las = 0 )
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_2, "_", 
	scenario, "_without a legend" ), path = MapsDir, type = "png" )
pause( 2 ); graphics.off()

######## Produce a map of the approximated spatial correlation for the tinyVAST model using the SPDE method 
######## and covariate-based regional anisotropy under Scenario 4 for the third location of interest and save it
scenario <- "Regional_Anisotropy_Scenario4"
plot_sf <- r_Regional_Anisotropy_Scenario4_3
zlim <- range( c( r_Regional_Anisotropy_Scenario1_1$r1, r_Regional_Anisotropy_Scenario1_2$r1, 
	r_Regional_Anisotropy_Scenario1_3$r1, 
	r_Regional_Anisotropy_Scenario2_1$r1, r_Regional_Anisotropy_Scenario2_2$r1, 
	r_Regional_Anisotropy_Scenario2_3$r1, 
	r_Regional_Anisotropy_Scenario3_1$r1, r_Regional_Anisotropy_Scenario3_2$r1, 
	r_Regional_Anisotropy_Scenario3_3$r1, 
	r_Regional_Anisotropy_Scenario4_1$r1, r_Regional_Anisotropy_Scenario4_2$r1, 
	r_Regional_Anisotropy_Scenario4_3$r1 ), na.rm = TRUE ) 

#### Produce a version of the map with a legend and save it
new.graph( 0.5 )
par( oma = c( 3, 3, 0, 0 ), cex = 1.1 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, main = "", axes = FALSE, 
	xlim = st_bbox(plot_sf)[c( 1, 3 )], ylim = st_bbox(plot_sf)[c( 2, 4 )] ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
plot( contour_200m, col = "white", lwd = 3.5, lty = 2, add = TRUE )
plot( contour_500m, col = "gray48", lwd = 3.5, lty = 2, add = TRUE )
box()
axis( side = 1, at = Mapping_parameters[["xaxis"]], labels = Mapping_parameters[["xaxis_label"]] )
axis( side = 2, at = Mapping_parameters[["yaxis"]], labels = Mapping_parameters[["yaxis_label"]] )
minor.tick( nx = 2, ny = 1 )
mtext( side = 1, outer = TRUE, "Longitude (°E)", cex = 1.25, line = par()$oma[1]/2, padj = 0 )
mtext( side = 2, outer = TRUE, "Latitude (°N)", cex = 1.25, line = par()$oma[2]/2, padj = 0, las = 0 )
xl = ( 1 - Mapping_parameters[["legend_x"]][1] ) * par( 'usr' )[1] + 
	( Mapping_parameters[["legend_x"]][1] ) * par( 'usr' )[2]
xr = ( 1 - Mapping_parameters[["legend_x"]][2] ) * par( 'usr' )[1] + 
	( Mapping_parameters[["legend_x"]][2] ) * par( 'usr' )[2]
yb = ( 1 - Mapping_parameters[["legend_y"]][1] ) * par( 'usr' )[3] + 
	( Mapping_parameters[["legend_y"]][1] ) * par( 'usr' )[4]
yt = ( 1 - Mapping_parameters[["legend_y"]][2] ) * par( 'usr' )[3] + 
	( Mapping_parameters[["legend_y"]][2] ) * par( 'usr' )[4]
align = c( "lt","rb" )[2]
gradient = c( "x", "y" )[2]
plotrix::color.legend( xl = xl, yb = yb, xr = xr, yt = yt, 
	legend = round( seq( zlim[1], zlim[2], length = 5 ), digits = 2 ),
	rect.col = sf.colors( n = 10, cutoff.tails = c( 0.35, 0.2 ), alpha = 1, categorical = FALSE ),
	cex = 1.2, align = align, gradient = gradient, col = "black" )
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_3, "_", 
	scenario, "_with a legend" ), path = MapsDir, type = "png" )
pause( 2 ); graphics.off()

#### Produce a version of the map without a legend and save it 
new.graph( 0.5 )
par( oma = c( 3, 3, 0, 0 ), cex = 1.1 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, main = "", axes = FALSE, 
	xlim = st_bbox(plot_sf)[c( 1, 3 )], ylim = st_bbox(plot_sf)[c( 2, 4 )] ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
plot( contour_200m, col = "white", lwd = 3.5, lty = 2, add = TRUE )
plot( contour_500m, col = "gray48", lwd = 3.5, lty = 2, add = TRUE )
box()
axis( side = 1, at = Mapping_parameters[["xaxis"]], labels = Mapping_parameters[["xaxis_label"]] )
axis( side = 2, at = Mapping_parameters[["yaxis"]], labels = Mapping_parameters[["yaxis_label"]] )
minor.tick( nx = 2, ny = 1 )
mtext( side = 1, outer = TRUE, "Longitude (°E)", cex = 1.25, line = par()$oma[1]/2, padj = 0 )
mtext( side = 2, outer = TRUE, "Latitude (°N)", cex = 1.25, line = par()$oma[2]/2, padj = 0, las = 0 )
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_3, "_", 
	scenario, "_without a legend" ), path = MapsDir, type = "png" )
pause( 2 ); graphics.off()


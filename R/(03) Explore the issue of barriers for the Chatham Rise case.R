###############################################################################################################################
##
##  (03) This script explores the issue of barriers for the Chatham Rise case
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
MapsDir <- paste0( DIR$Maps, "Chatham Rise case_2" )
if ( !dir.exists( MapsDir ) ) { dir.create( MapsDir ) } 

######## Load and process shapefiles
contour_200m <- read_sf( dsn = make.filename( "NZ_200_m.shp", DIR$Input ) )
contour_200m <- st_transform( contour_200m, crs = "+init=epsg:4326" )
contour_200m <- sf::st_shift_longitude( contour_200m )
contour_500m <- read_sf( dsn = make.filename( "NZ_500_m.shp", DIR$Input ) )
contour_500m <- st_transform( contour_500m, crs = "+init=epsg:4326" )
contour_500m <- sf::st_shift_longitude( contour_500m )

######## Load mapping parameters 
load( file = make.filename( "Mapping_parameters_Chatham_Rise_hake.RData", DIR$Input ) )

####### Define or modify a few parameters/settings
location_1 <- "Shallow Chatham Rise"
location_2 <- "Middle-depth Chatham Rise"
location_3 <- "Deep Chatham Rise"
location_of_interest_1 <- c( 183.2, -43.8 )
location_of_interest_2 <- c( 181, -44 )
location_of_interest_3 <- c( 179, -45 )
Mapping_parameters[["legend_x"]] <- c( 0.89, 0.94 ) 
Mapping_parameters[["legend_y"]] <- c( 0.02, 0.22 ) 
Mapping_parameters[["xlim_map"]] <- c( 172.5, 187.0 )
map_sf <- Mapping_parameters[["land_data"]]
cutoff <- 0.1
barrier_stiffness <- 0.01
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

######## Define the mesh for the tinyVAST model using the SPDE method and assuming isotropy
mesh <- fm_mesh_2d( Data[,c( 'Lon', 'Lat' )], cutoff = cutoff )
plot( mesh )
pause( 3 ); graphics.off()

######## Fit the tinyVAST model using the SPDE method and assuming isotropy
SPDE_Isotropy <- tinyVAST( 
	data = droplevels( Data ), 
	formula = formula, 
	space_term = space_term,
  	spacetime_term = spacetime_term,
  	time_term = time_term,
  	spatial_domain = mesh,
  	family = tweedie( "log" ),
	space_columns = c( "Lon", "Lat" ),
  	time_column = "Year",
  	times = min( Data$Year ) : max( Data$Year ),
  	variable_column = "var",
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

######## Approximate spatial correlation for the tinyVAST model using the SPDE method and assuming isotropy
######## for the first location of interest 
mesh_sf <- fm_as_sfc( mesh )
domain_sf <- st_union( mesh_sf )
grid_sf <- st_make_grid( domain_sf, cellsize = cellsize )
grid_sf <- st_intersection( grid_sf, domain_sf )
pred <- st_coordinates( st_centroid( grid_sf ) )
fit <- SPDE_Isotropy
r1 <- spatial_cor(
  	fit$rep$Q,
  	mesh = mesh,
  	coord = location_of_interest_1,
  	pred = pred
)
plot_sf <- st_sf( st_geometry( grid_sf ), r1 = r1 )
r_SPDE_Isotropy_1 <- plot_sf

plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
pause( 3 ); graphics.off()

######## Approximate spatial correlation for the tinyVAST model using the SPDE method and assuming isotropy
######## for the second location of interest 
mesh_sf <- fm_as_sfc( mesh )
domain_sf <- st_union( mesh_sf )
grid_sf <- st_make_grid( domain_sf, cellsize = cellsize )
grid_sf <- st_intersection( grid_sf, domain_sf )
pred <- st_coordinates( st_centroid( grid_sf ) )
fit <- SPDE_Isotropy
r1 <- spatial_cor(
  	fit$rep$Q,
  	mesh = mesh,
  	coord = location_of_interest_2,
  	pred = pred
)
plot_sf <- st_sf( st_geometry( grid_sf ), r1 = r1 )
r_SPDE_Isotropy_2 <- plot_sf

plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
pause( 3 ); graphics.off()

######## Approximate spatial correlation for the tinyVAST model using the SPDE method and assuming isotropy
######## for the third location of interest 
mesh_sf <- fm_as_sfc( mesh )
domain_sf <- st_union( mesh_sf )
grid_sf <- st_make_grid( domain_sf, cellsize = cellsize )
grid_sf <- st_intersection( grid_sf, domain_sf )
pred <- st_coordinates( st_centroid( grid_sf ) )
fit <- SPDE_Isotropy
r1 <- spatial_cor(
  	fit$rep$Q,
  	mesh = mesh,
  	coord = location_of_interest_3,
  	pred = pred
)
plot_sf <- st_sf( st_geometry( grid_sf ), r1 = r1 )
r_SPDE_Isotropy_3 <- plot_sf

plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
pause( 3 ); graphics.off()

######## Generate a mesh for the tinyVAST model using the SPDE-Barrier method under Scenario 1
mesh <- fm_mesh_2d( Data[,c( "Lon", "Lat" )], cutoff = cutoff )

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

Data$Region <- ifelse( Depth > 800, "DEEP", "MIDDEPTH" )
Data$Region <- ifelse( Depth <= 500, "SHALLOW", Data$Region )
Table( Data$Region )
####      DEEP  MIDDEPTH   SHALLOW TOTAL.OBS 
####       394      1115      1991      3500 

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
nn <- RANN::nn2( data = Data[,c( 'Lon', 'Lat' )], query = mesh_cov$triangle_covariates[,1 : 2], k = 1 )
Region <- Data[nn$nn.idx, 'Region']
mesh_cov$triangle_covariates$DEPTH_STRATUM <- ifelse( Region == "SHALLOW", 0.5, 0 )
mesh_cov$triangle_covariates$DEPTH_STRATUM <- ifelse( Region == "MIDDEPTH", 0.25, mesh_cov$triangle_covariates$DEPTH_STRATUM )

plot( mesh_cov )
pause( 3 ); graphics.off()

######## Prepare the tinyVAST models using SPDE-Barrier method under Scenario 1
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
		triangle_formula = ~ offset( DEPTH_STRATUM )
  	),
	control = tinyVASTcontrol( 
		verbose = TRUE, 
		opt_loops = 0,
		newton_loops = 0, 
		trace = 1,
    		use_anisotropy = FALSE, 
		getsd = FALSE,
    		calculate = FALSE, 
		barrier_stiffness = barrier_stiffness, 
	)
)
tmb_par <- fit_norun$tmb_inputs$tmb_par

######## Fit the tinyVAST model using SPDE-Barrier method under Scenario 1
SPDE_Barrier_Scenario1 <- tinyVAST( 
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
		triangle_formula = ~ offset( DEPTH_STRATUM )
  	),
	control = tinyVASTcontrol( 
		verbose = TRUE, 
		opt_loops = 0,
		newton_loops = 0, 
		trace = 1,
    		use_anisotropy = FALSE, 
		getsd = FALSE,
    		calculate = FALSE, 
		barrier_stiffness = barrier_stiffness, 
		tmb_par = tmb_par 
	)
)

######## Approximate spatial correlation for the tinyVAST model using the SPDE-Barrier method under Scenario 1
######## for the first location of interest 
mesh_sf <- fm_as_sfc( mesh )
domain_sf <- st_union( mesh_sf )
grid_sf <- st_make_grid( domain_sf, cellsize = cellsize )
grid_sf <- st_intersection( grid_sf, domain_sf )
pred <- st_coordinates( st_centroid( grid_sf ) )
fit <- SPDE_Barrier_Scenario1
r1 <- spatial_cor(
  	fit$rep$Q,
  	mesh = mesh,
  	coord = location_of_interest_1,
  	pred = pred
)
plot_sf <- st_sf( st_geometry( grid_sf ), r1 = r1 )
r_SPDE_Barrier_Scenario1_1 <- plot_sf

plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
pause( 3 ); graphics.off()

######## Approximate spatial correlation for the tinyVAST model using the SPDE-Barrier method under Scenario 1
######## for the second location of interest 
mesh_sf <- fm_as_sfc( mesh )
domain_sf <- st_union( mesh_sf )
grid_sf <- st_make_grid( domain_sf, cellsize = cellsize )
grid_sf <- st_intersection( grid_sf, domain_sf )
pred <- st_coordinates( st_centroid( grid_sf ) )
fit <- SPDE_Barrier_Scenario1
r1 <- spatial_cor(
  	fit$rep$Q,
  	mesh = mesh,
  	coord = location_of_interest_2,
  	pred = pred
)
plot_sf <- st_sf( st_geometry( grid_sf ), r1 = r1 )
r_SPDE_Barrier_Scenario1_2 <- plot_sf

plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
pause( 3 ); graphics.off()

######## Approximate spatial correlation for the tinyVAST model using the SPDE-Barrier method under Scenario 1
######## for the third location of interest 
mesh_sf <- fm_as_sfc( mesh )
domain_sf <- st_union( mesh_sf )
grid_sf <- st_make_grid( domain_sf, cellsize = cellsize )
grid_sf <- st_intersection( grid_sf, domain_sf )
pred <- st_coordinates( st_centroid( grid_sf ) )
fit <- SPDE_Barrier_Scenario1
r1 <- spatial_cor(
  	fit$rep$Q,
  	mesh = mesh,
  	coord = location_of_interest_3,
  	pred = pred
)
plot_sf <- st_sf( st_geometry( grid_sf ), r1 = r1 )
r_SPDE_Barrier_Scenario1_3 <- plot_sf

plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
pause( 3 ); graphics.off()

######## Generate a mesh for the tinyVAST model using the SPDE-Barrier method under Scenario 2
mesh <- fm_mesh_2d( Data[,c( "Lon", "Lat" )], cutoff = cutoff )

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

Data$Region <- ifelse( Depth > 800, "DEEP", "MIDDEPTH" )
Data$Region <- ifelse( Depth <= 500, "SHALLOW", Data$Region )
Table( Data$Region )
####      DEEP  MIDDEPTH   SHALLOW TOTAL.OBS 
####       394      1115      1991      3500 

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
nn <- RANN::nn2( data = Data[,c( 'Lon', 'Lat' )], query = mesh_cov$triangle_covariates[,1 : 2], k = 1 )
Region <- Data[nn$nn.idx, 'Region']
mesh_cov$triangle_covariates$DEPTH_STRATUM <- ifelse( Region == "DEEP", 0.5, 0 )
mesh_cov$triangle_covariates$DEPTH_STRATUM <- ifelse( Region == "MIDDEPTH", 0.25, mesh_cov$triangle_covariates$DEPTH_STRATUM )

plot( mesh_cov )
pause( 3 ); graphics.off()

######## Prepare the tinyVAST models using SPDE-Barrier method under Scenario 2
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
		triangle_formula = ~ offset( DEPTH_STRATUM )
  	),
	control = tinyVASTcontrol( 
		verbose = TRUE, 
		opt_loops = 0,
		newton_loops = 0, 
		trace = 1,
    		use_anisotropy = FALSE, 
		getsd = FALSE,
    		calculate = FALSE, 
		barrier_stiffness = barrier_stiffness, 
	)
)
tmb_par <- fit_norun$tmb_inputs$tmb_par

######## Fit the tinyVAST model using SPDE-Barrier method under Scenario 2
SPDE_Barrier_Scenario2 <- tinyVAST( 
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
		triangle_formula = ~ offset( DEPTH_STRATUM )
  	),
	control = tinyVASTcontrol( 
		verbose = TRUE, 
		opt_loops = 0,
		newton_loops = 0, 
		trace = 1,
    		use_anisotropy = FALSE, 
		getsd = FALSE,
    		calculate = FALSE, 
		barrier_stiffness = barrier_stiffness, 
		tmb_par = tmb_par 
	)
)

######## Approximate spatial correlation for the tinyVAST model using the SPDE-Barrier method under Scenario 2
######## for the first location of interest 
mesh_sf <- fm_as_sfc( mesh )
domain_sf <- st_union( mesh_sf )
grid_sf <- st_make_grid( domain_sf, cellsize = cellsize )
grid_sf <- st_intersection( grid_sf, domain_sf )
pred <- st_coordinates( st_centroid( grid_sf ) )
fit <- SPDE_Barrier_Scenario2
r1 <- spatial_cor(
  	fit$rep$Q,
  	mesh = mesh,
  	coord = location_of_interest_1,
  	pred = pred
)
plot_sf <- st_sf( st_geometry( grid_sf ), r1 = r1 )
r_SPDE_Barrier_Scenario2_1 <- plot_sf

plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
pause( 3 ); graphics.off()

######## Approximate spatial correlation for the tinyVAST model using the SPDE-Barrier method under Scenario 2
######## for the second location of interest 
mesh_sf <- fm_as_sfc( mesh )
domain_sf <- st_union( mesh_sf )
grid_sf <- st_make_grid( domain_sf, cellsize = cellsize )
grid_sf <- st_intersection( grid_sf, domain_sf )
pred <- st_coordinates( st_centroid( grid_sf ) )
fit <- SPDE_Barrier_Scenario2
r1 <- spatial_cor(
  	fit$rep$Q,
  	mesh = mesh,
  	coord = location_of_interest_2,
  	pred = pred
)
plot_sf <- st_sf( st_geometry( grid_sf ), r1 = r1 )
r_SPDE_Barrier_Scenario2_2 <- plot_sf

plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
pause( 3 ); graphics.off()

######## Approximate spatial correlation for the tinyVAST model using the SPDE-Barrier method under Scenario 2
######## for the third location of interest 
mesh_sf <- fm_as_sfc( mesh )
domain_sf <- st_union( mesh_sf )
grid_sf <- st_make_grid( domain_sf, cellsize = cellsize )
grid_sf <- st_intersection( grid_sf, domain_sf )
pred <- st_coordinates( st_centroid( grid_sf ) )
fit <- SPDE_Barrier_Scenario2
r1 <- spatial_cor(
  	fit$rep$Q,
  	mesh = mesh,
  	coord = location_of_interest_3,
  	pred = pred
)
plot_sf <- st_sf( st_geometry( grid_sf ), r1 = r1 )
r_SPDE_Barrier_Scenario2_3 <- plot_sf

plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
pause( 3 ); graphics.off()

######## Produce a map of the approximated spatial correlation for the tinyVAST model using the SPDE method 
######## and assuming isotropy for the first location of interest and save it
scenario <- "SPDE_Isotropy"
plot_sf <- r_SPDE_Isotropy_1
zlim <- range( c( r_SPDE_Isotropy_1$r1, r_SPDE_Isotropy_2$r1, r_SPDE_Isotropy_3$r1, 
	r_SPDE_Barrier_Scenario1_1$r1, r_SPDE_Barrier_Scenario1_2$r1, r_SPDE_Barrier_Scenario1_3$r1, 
	r_SPDE_Barrier_Scenario2_1$r1, r_SPDE_Barrier_Scenario2_2$r1, r_SPDE_Barrier_Scenario2_3$r1 ), na.rm = TRUE ) 

#### Produce a version of the map with a legend and save it 
new.graph( 0.5 )
par( mar = c( 4.5, 4.5, 1, 1 ) )
plot( Mapping_parameters[["xlim_map"]], Mapping_parameters[["ylim_map"]], 
	xlim = Mapping_parameters[["xlim_map"]], ylim = Mapping_parameters[["ylim_map"]], 
	type = "n", ylab = "Latitude (°N)", xlab = "Longitude (°E)", cex.lab = 1.25, font.lab = 2 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, add = TRUE ) 
plot( map_sf, col = makeTransparent( "gray88", alpha = 170 ), border = NA, add = TRUE )
plot( contour_200m, col = "white", lwd = 3.5, lty = 2, add = TRUE )
plot( contour_500m, col = "gray48", lwd = 3.5, lty = 2, add = TRUE )
box()
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
	cex = 1., align = align, gradient = gradient, col = "black" )
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_1, "_", 
	scenario, "_with a legend" ), path = MapsDir, type = "png" )
pause( 2 ); graphics.off()

#### Produce a version of the map without a legend and save it 
new.graph( 0.5 )
par( mar = c( 4.5, 4.5, 1, 1 ) )
plot( Mapping_parameters[["xlim_map"]], Mapping_parameters[["ylim_map"]], 
	xlim = Mapping_parameters[["xlim_map"]], ylim = Mapping_parameters[["ylim_map"]], 
	type = "n", ylab = "Latitude (°N)", xlab = "Longitude (°E)", cex.lab = 1.25, font.lab = 2 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, add = TRUE ) 
plot( map_sf, col = makeTransparent( "gray88", alpha = 170 ), border = NA, add = TRUE )
plot( contour_200m, col = "white", lwd = 3.5, lty = 2, add = TRUE )
plot( contour_500m, col = "gray48", lwd = 3.5, lty = 2, add = TRUE )
box()
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_1, "_", 
	scenario, "_without a legend" ), path = MapsDir, type = "png" )
pause( 2 ); graphics.off()

######## Produce a map of the approximated spatial correlation for the tinyVAST model using the SPDE method 
######## and assuming isotropy for the second location of interest and save it
scenario <- "SPDE_Isotropy"
plot_sf <- r_SPDE_Isotropy_2
zlim <- range( c( r_SPDE_Isotropy_1$r1, r_SPDE_Isotropy_2$r1, r_SPDE_Isotropy_3$r1, 
	r_SPDE_Barrier_Scenario1_1$r1, r_SPDE_Barrier_Scenario1_2$r1, r_SPDE_Barrier_Scenario1_3$r1, 
	r_SPDE_Barrier_Scenario2_1$r1, r_SPDE_Barrier_Scenario2_2$r1, r_SPDE_Barrier_Scenario2_3$r1 ), na.rm = TRUE ) 

#### Produce a version of the map with a legend and save it 
new.graph( 0.5 )
par( mar = c( 4.5, 4.5, 1, 1 ) )
plot( Mapping_parameters[["xlim_map"]], Mapping_parameters[["ylim_map"]], 
	xlim = Mapping_parameters[["xlim_map"]], ylim = Mapping_parameters[["ylim_map"]], 
	type = "n", ylab = "Latitude (°N)", xlab = "Longitude (°E)", cex.lab = 1.25, font.lab = 2 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, add = TRUE ) 
plot( map_sf, col = makeTransparent( "gray88", alpha = 170 ), border = NA, add = TRUE )
plot( contour_200m, col = "white", lwd = 3.5, lty = 2, add = TRUE )
plot( contour_500m, col = "gray48", lwd = 3.5, lty = 2, add = TRUE )
box()
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
	cex = 1., align = align, gradient = gradient, col = "black" )
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_2, "_", 
	scenario, "_with a legend" ), path = MapsDir, type = "png" )
pause( 2 ); graphics.off()

#### Produce a version of the map without a legend and save it 
new.graph( 0.5 )
par( mar = c( 4.5, 4.5, 1, 1 ) )
plot( Mapping_parameters[["xlim_map"]], Mapping_parameters[["ylim_map"]], 
	xlim = Mapping_parameters[["xlim_map"]], ylim = Mapping_parameters[["ylim_map"]], 
	type = "n", ylab = "Latitude (°N)", xlab = "Longitude (°E)", cex.lab = 1.25, font.lab = 2 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, add = TRUE ) 
plot( map_sf, col = makeTransparent( "gray88", alpha = 170 ), border = NA, add = TRUE )
plot( contour_200m, col = "white", lwd = 3.5, lty = 2, add = TRUE )
plot( contour_500m, col = "gray48", lwd = 3.5, lty = 2, add = TRUE )
box()
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_2, "_", 
	scenario, "_without a legend" ), path = MapsDir, type = "png" )
pause( 2 ); graphics.off()

######## Produce a map of the approximated spatial correlation for the tinyVAST model using the SPDE method 
######## and assuming isotropy for the third location of interest and save it
scenario <- "SPDE_Isotropy"
plot_sf <- r_SPDE_Isotropy_3
zlim <- range( c( r_SPDE_Isotropy_1$r1, r_SPDE_Isotropy_2$r1, r_SPDE_Isotropy_3$r1, 
	r_SPDE_Barrier_Scenario1_1$r1, r_SPDE_Barrier_Scenario1_2$r1, r_SPDE_Barrier_Scenario1_3$r1, 
	r_SPDE_Barrier_Scenario2_1$r1, r_SPDE_Barrier_Scenario2_2$r1, r_SPDE_Barrier_Scenario2_3$r1 ), na.rm = TRUE ) 

#### Produce a version of the map with a legend and save it 
new.graph( 0.5 )
par( mar = c( 4.5, 4.5, 1, 1 ) )
plot( Mapping_parameters[["xlim_map"]], Mapping_parameters[["ylim_map"]], 
	xlim = Mapping_parameters[["xlim_map"]], ylim = Mapping_parameters[["ylim_map"]], 
	type = "n", ylab = "Latitude (°N)", xlab = "Longitude (°E)", cex.lab = 1.25, font.lab = 2 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, add = TRUE ) 
plot( map_sf, col = makeTransparent( "gray88", alpha = 170 ), border = NA, add = TRUE )
plot( contour_200m, col = "white", lwd = 3.5, lty = 2, add = TRUE )
plot( contour_500m, col = "gray48", lwd = 3.5, lty = 2, add = TRUE )
box()
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
	cex = 1., align = align, gradient = gradient, col = "black" )
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_3, "_", 
	scenario, "_with a legend" ), path = MapsDir, type = "png" )
pause( 2 ); graphics.off()

#### Produce a version of the map without a legend and save it 
new.graph( 0.5 )
par( mar = c( 4.5, 4.5, 1, 1 ) )
plot( Mapping_parameters[["xlim_map"]], Mapping_parameters[["ylim_map"]], 
	xlim = Mapping_parameters[["xlim_map"]], ylim = Mapping_parameters[["ylim_map"]], 
	type = "n", ylab = "Latitude (°N)", xlab = "Longitude (°E)", cex.lab = 1.25, font.lab = 2 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, add = TRUE ) 
plot( map_sf, col = makeTransparent( "gray88", alpha = 170 ), border = NA, add = TRUE )
plot( contour_200m, col = "white", lwd = 3.5, lty = 2, add = TRUE )
plot( contour_500m, col = "gray48", lwd = 3.5, lty = 2, add = TRUE )
box()
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_3, "_", 
	scenario, "_without a legend" ), path = MapsDir, type = "png" )
pause( 2 ); graphics.off()

######## Produce a map of the approximated spatial correlation for the tinyVAST model using the SPDE-Barrier method 
######## under Scenario 1 for the first location of interest and save it
scenario <- "SPDE_Barrier_Scenario1"
plot_sf <- r_SPDE_Barrier_Scenario1_1
zlim <- range( c( r_SPDE_Isotropy_1$r1, r_SPDE_Isotropy_2$r1, r_SPDE_Isotropy_3$r1, 
	r_SPDE_Barrier_Scenario1_1$r1, r_SPDE_Barrier_Scenario1_2$r1, r_SPDE_Barrier_Scenario1_3$r1, 
	r_SPDE_Barrier_Scenario2_1$r1, r_SPDE_Barrier_Scenario2_2$r1, r_SPDE_Barrier_Scenario2_3$r1 ), na.rm = TRUE ) 

#### Produce a version of the map with a legend and save it 
new.graph( 0.5 )
par( mar = c( 4.5, 4.5, 1, 1 ) )
plot( Mapping_parameters[["xlim_map"]], Mapping_parameters[["ylim_map"]], 
	xlim = Mapping_parameters[["xlim_map"]], ylim = Mapping_parameters[["ylim_map"]], 
	type = "n", ylab = "Latitude (°N)", xlab = "Longitude (°E)", cex.lab = 1.25, font.lab = 2 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, add = TRUE ) 
plot( map_sf, col = makeTransparent( "gray88", alpha = 170 ), border = NA, add = TRUE )
plot( contour_200m, col = "white", lwd = 3.5, lty = 2, add = TRUE )
plot( contour_500m, col = "gray48", lwd = 3.5, lty = 2, add = TRUE )
box()
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
	cex = 1., align = align, gradient = gradient, col = "black" )
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_1, "_", 
	scenario, "_with a legend" ), path = MapsDir, type = "png" )
pause( 2 ); graphics.off()

#### Produce a version of the map without a legend and save it 
new.graph( 0.5 )
par( mar = c( 4.5, 4.5, 1, 1 ) )
plot( Mapping_parameters[["xlim_map"]], Mapping_parameters[["ylim_map"]], 
	xlim = Mapping_parameters[["xlim_map"]], ylim = Mapping_parameters[["ylim_map"]], 
	type = "n", ylab = "Latitude (°N)", xlab = "Longitude (°E)", cex.lab = 1.25, font.lab = 2 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, add = TRUE ) 
plot( map_sf, col = makeTransparent( "gray88", alpha = 170 ), border = NA, add = TRUE )
plot( contour_200m, col = "white", lwd = 3.5, lty = 2, add = TRUE )
plot( contour_500m, col = "gray48", lwd = 3.5, lty = 2, add = TRUE )
box()
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_1, "_", 
	scenario, "_without a legend" ), path = MapsDir, type = "png" )
pause( 2 ); graphics.off()

######## Produce a map of the approximated spatial correlation for the tinyVAST model using the SPDE-Barrier method 
######## under Scenario 1 for the second location of interest and save it
scenario <- "SPDE_Barrier_Scenario1"
plot_sf <- r_SPDE_Barrier_Scenario1_2
zlim <- range( c( r_SPDE_Isotropy_1$r1, r_SPDE_Isotropy_2$r1, r_SPDE_Isotropy_3$r1, 
	r_SPDE_Barrier_Scenario1_1$r1, r_SPDE_Barrier_Scenario1_2$r1, r_SPDE_Barrier_Scenario1_3$r1, 
	r_SPDE_Barrier_Scenario2_1$r1, r_SPDE_Barrier_Scenario2_2$r1, r_SPDE_Barrier_Scenario2_3$r1 ), na.rm = TRUE ) 

#### Produce a version of the map with a legend and save it 
new.graph( 0.5 )
par( mar = c( 4.5, 4.5, 1, 1 ) )
plot( Mapping_parameters[["xlim_map"]], Mapping_parameters[["ylim_map"]], 
	xlim = Mapping_parameters[["xlim_map"]], ylim = Mapping_parameters[["ylim_map"]], 
	type = "n", ylab = "Latitude (°N)", xlab = "Longitude (°E)", cex.lab = 1.25, font.lab = 2 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, add = TRUE ) 
plot( map_sf, col = makeTransparent( "gray88", alpha = 170 ), border = NA, add = TRUE )
plot( contour_200m, col = "white", lwd = 3.5, lty = 2, add = TRUE )
plot( contour_500m, col = "gray48", lwd = 3.5, lty = 2, add = TRUE )
box()
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
	cex = 1., align = align, gradient = gradient, col = "black" )
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_2, "_", 
	scenario, "_with a legend" ), path = MapsDir, type = "png" )
pause( 2 ); graphics.off()

#### Produce a version of the map without a legend and save it 
new.graph( 0.5 )
par( mar = c( 4.5, 4.5, 1, 1 ) )
plot( Mapping_parameters[["xlim_map"]], Mapping_parameters[["ylim_map"]], 
	xlim = Mapping_parameters[["xlim_map"]], ylim = Mapping_parameters[["ylim_map"]], 
	type = "n", ylab = "Latitude (°N)", xlab = "Longitude (°E)", cex.lab = 1.25, font.lab = 2 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, add = TRUE ) 
plot( map_sf, col = makeTransparent( "gray88", alpha = 170 ), border = NA, add = TRUE )
plot( contour_200m, col = "white", lwd = 3.5, lty = 2, add = TRUE )
plot( contour_500m, col = "gray48", lwd = 3.5, lty = 2, add = TRUE )
box()
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_2, "_", 
	scenario, "_without a legend" ), path = MapsDir, type = "png" )
pause( 2 ); graphics.off()

######## Produce a map of the approximated spatial correlation for the tinyVAST model using the SPDE-Barrier method 
######## under Scenario 1 for the third location of interest and save it
scenario <- "SPDE_Barrier_Scenario1"
plot_sf <- r_SPDE_Barrier_Scenario1_3
zlim <- range( c( r_SPDE_Isotropy_1$r1, r_SPDE_Isotropy_2$r1, r_SPDE_Isotropy_3$r1, 
	r_SPDE_Barrier_Scenario1_1$r1, r_SPDE_Barrier_Scenario1_2$r1, r_SPDE_Barrier_Scenario1_3$r1, 
	r_SPDE_Barrier_Scenario2_1$r1, r_SPDE_Barrier_Scenario2_2$r1, r_SPDE_Barrier_Scenario2_3$r1 ), na.rm = TRUE ) 

#### Produce a version of the map with a legend and save it 
new.graph( 0.5 )
par( mar = c( 4.5, 4.5, 1, 1 ) )
plot( Mapping_parameters[["xlim_map"]], Mapping_parameters[["ylim_map"]], 
	xlim = Mapping_parameters[["xlim_map"]], ylim = Mapping_parameters[["ylim_map"]], 
	type = "n", ylab = "Latitude (°N)", xlab = "Longitude (°E)", cex.lab = 1.25, font.lab = 2 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, add = TRUE ) 
plot( map_sf, col = makeTransparent( "gray88", alpha = 170 ), border = NA, add = TRUE )
plot( contour_200m, col = "white", lwd = 3.5, lty = 2, add = TRUE )
plot( contour_500m, col = "gray48", lwd = 3.5, lty = 2, add = TRUE )
box()
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
	cex = 1., align = align, gradient = gradient, col = "black" )
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_3, "_", 
	scenario, "_with a legend" ), path = MapsDir, type = "png" )
pause( 2 ); graphics.off()

#### Produce a version of the map without a legend and save it 
new.graph( 0.5 )
par( mar = c( 4.5, 4.5, 1, 1 ) )
plot( Mapping_parameters[["xlim_map"]], Mapping_parameters[["ylim_map"]], 
	xlim = Mapping_parameters[["xlim_map"]], ylim = Mapping_parameters[["ylim_map"]], 
	type = "n", ylab = "Latitude (°N)", xlab = "Longitude (°E)", cex.lab = 1.25, font.lab = 2 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, add = TRUE ) 
plot( map_sf, col = makeTransparent( "gray88", alpha = 170 ), border = NA, add = TRUE )
plot( contour_200m, col = "white", lwd = 3.5, lty = 2, add = TRUE )
plot( contour_500m, col = "gray48", lwd = 3.5, lty = 2, add = TRUE )
box()
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_3, "_", 
	scenario, "_without a legend" ), path = MapsDir, type = "png" )
pause( 2 ); graphics.off()

######## Produce a map of the approximated spatial correlation for the tinyVAST model using the SPDE-Barrier method 
######## under Scenario 2 for the first location of interest and save it
scenario <- "SPDE_Barrier_Scenario2"
plot_sf <- r_SPDE_Barrier_Scenario2_1
zlim <- range( c( r_SPDE_Isotropy_1$r1, r_SPDE_Isotropy_2$r1, r_SPDE_Isotropy_3$r1, 
	r_SPDE_Barrier_Scenario1_1$r1, r_SPDE_Barrier_Scenario1_2$r1, r_SPDE_Barrier_Scenario1_3$r1, 
	r_SPDE_Barrier_Scenario2_1$r1, r_SPDE_Barrier_Scenario2_2$r1, r_SPDE_Barrier_Scenario2_3$r1 ), na.rm = TRUE ) 

#### Produce a version of the map with a legend and save it 
new.graph( 0.5 )
par( mar = c( 4.5, 4.5, 1, 1 ) )
plot( Mapping_parameters[["xlim_map"]], Mapping_parameters[["ylim_map"]], 
	xlim = Mapping_parameters[["xlim_map"]], ylim = Mapping_parameters[["ylim_map"]], 
	type = "n", ylab = "Latitude (°N)", xlab = "Longitude (°E)", cex.lab = 1.25, font.lab = 2 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, add = TRUE ) 
plot( map_sf, col = makeTransparent( "gray88", alpha = 170 ), border = NA, add = TRUE )
plot( contour_200m, col = "white", lwd = 3.5, lty = 2, add = TRUE )
plot( contour_500m, col = "gray48", lwd = 3.5, lty = 2, add = TRUE )
box()
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
	cex = 1., align = align, gradient = gradient, col = "black" )
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_1, "_", 
	scenario, "_with a legend" ), path = MapsDir, type = "png" )
pause( 2 ); graphics.off()

#### Produce a version of the map without a legend and save it 
new.graph( 0.5 )
par( mar = c( 4.5, 4.5, 1, 1 ) )
plot( Mapping_parameters[["xlim_map"]], Mapping_parameters[["ylim_map"]], 
	xlim = Mapping_parameters[["xlim_map"]], ylim = Mapping_parameters[["ylim_map"]], 
	type = "n", ylab = "Latitude (°N)", xlab = "Longitude (°E)", cex.lab = 1.25, font.lab = 2 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, add = TRUE ) 
plot( map_sf, col = makeTransparent( "gray88", alpha = 170 ), border = NA, add = TRUE )
plot( contour_200m, col = "white", lwd = 3.5, lty = 2, add = TRUE )
plot( contour_500m, col = "gray48", lwd = 3.5, lty = 2, add = TRUE )
box()
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_1, "_", 
	scenario, "_without a legend" ), path = MapsDir, type = "png" )
pause( 2 ); graphics.off()

######## Produce a map of the approximated spatial correlation for the tinyVAST model using the SPDE-Barrier method 
######## under Scenario 2 for the second location of interest and save it
scenario <- "SPDE_Barrier_Scenario2"
plot_sf <- r_SPDE_Barrier_Scenario2_2
zlim <- range( c( r_SPDE_Isotropy_1$r1, r_SPDE_Isotropy_2$r1, r_SPDE_Isotropy_3$r1, 
	r_SPDE_Barrier_Scenario1_1$r1, r_SPDE_Barrier_Scenario1_2$r1, r_SPDE_Barrier_Scenario1_3$r1, 
	r_SPDE_Barrier_Scenario2_1$r1, r_SPDE_Barrier_Scenario2_2$r1, r_SPDE_Barrier_Scenario2_3$r1 ), na.rm = TRUE ) 

#### Produce a version of the map with a legend and save it 
new.graph( 0.5 )
par( mar = c( 4.5, 4.5, 1, 1 ) )
plot( Mapping_parameters[["xlim_map"]], Mapping_parameters[["ylim_map"]], 
	xlim = Mapping_parameters[["xlim_map"]], ylim = Mapping_parameters[["ylim_map"]], 
	type = "n", ylab = "Latitude (°N)", xlab = "Longitude (°E)", cex.lab = 1.25, font.lab = 2 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, add = TRUE ) 
plot( map_sf, col = makeTransparent( "gray88", alpha = 170 ), border = NA, add = TRUE )
plot( contour_200m, col = "white", lwd = 3.5, lty = 2, add = TRUE )
plot( contour_500m, col = "gray48", lwd = 3.5, lty = 2, add = TRUE )
box()
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
	cex = 1., align = align, gradient = gradient, col = "black" )
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_2, "_", 
	scenario, "_with a legend" ), path = MapsDir, type = "png" )
pause( 2 ); graphics.off()

#### Produce a version of the map without a legend and save it 
new.graph( 0.5 )
par( mar = c( 4.5, 4.5, 1, 1 ) )
plot( Mapping_parameters[["xlim_map"]], Mapping_parameters[["ylim_map"]], 
	xlim = Mapping_parameters[["xlim_map"]], ylim = Mapping_parameters[["ylim_map"]], 
	type = "n", ylab = "Latitude (°N)", xlab = "Longitude (°E)", cex.lab = 1.25, font.lab = 2 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, add = TRUE ) 
plot( map_sf, col = makeTransparent( "gray88", alpha = 170 ), border = NA, add = TRUE )
plot( contour_200m, col = "white", lwd = 3.5, lty = 2, add = TRUE )
plot( contour_500m, col = "gray48", lwd = 3.5, lty = 2, add = TRUE )
box()
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_2, "_", 
	scenario, "_without a legend" ), path = MapsDir, type = "png" )
pause( 2 ); graphics.off()

######## Produce a map of the approximated spatial correlation for the tinyVAST model using the SPDE-Barrier method 
######## under Scenario 2 for the third location of interest and save it
scenario <- "SPDE_Barrier_Scenario2"
plot_sf <- r_SPDE_Barrier_Scenario2_3
zlim <- range( c( r_SPDE_Isotropy_1$r1, r_SPDE_Isotropy_2$r1, r_SPDE_Isotropy_3$r1, 
	r_SPDE_Barrier_Scenario1_1$r1, r_SPDE_Barrier_Scenario1_2$r1, r_SPDE_Barrier_Scenario1_3$r1, 
	r_SPDE_Barrier_Scenario2_1$r1, r_SPDE_Barrier_Scenario2_2$r1, r_SPDE_Barrier_Scenario2_3$r1 ), na.rm = TRUE ) 

#### Produce a version of the map with a legend and save it 
new.graph( 0.5 )
par( mar = c( 4.5, 4.5, 1, 1 ) )
plot( Mapping_parameters[["xlim_map"]], Mapping_parameters[["ylim_map"]], 
	xlim = Mapping_parameters[["xlim_map"]], ylim = Mapping_parameters[["ylim_map"]], 
	type = "n", ylab = "Latitude (°N)", xlab = "Longitude (°E)", cex.lab = 1.25, font.lab = 2 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, add = TRUE ) 
plot( map_sf, col = makeTransparent( "gray88", alpha = 170 ), border = NA, add = TRUE )
plot( contour_200m, col = "white", lwd = 3.5, lty = 2, add = TRUE )
plot( contour_500m, col = "gray48", lwd = 3.5, lty = 2, add = TRUE )
box()
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
	cex = 1., align = align, gradient = gradient, col = "black" )
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_3, "_", 
	scenario, "_with a legend" ), path = MapsDir, type = "png" )
pause( 2 ); graphics.off()

#### Produce a version of the map without a legend and save it 
new.graph( 0.5 )
par( mar = c( 4.5, 4.5, 1, 1 ) )
plot( Mapping_parameters[["xlim_map"]], Mapping_parameters[["ylim_map"]], 
	xlim = Mapping_parameters[["xlim_map"]], ylim = Mapping_parameters[["ylim_map"]], 
	type = "n", ylab = "Latitude (°N)", xlab = "Longitude (°E)", cex.lab = 1.25, font.lab = 2 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, add = TRUE ) 
plot( map_sf, col = makeTransparent( "gray88", alpha = 170 ), border = NA, add = TRUE )
plot( contour_200m, col = "white", lwd = 3.5, lty = 2, add = TRUE )
plot( contour_500m, col = "gray48", lwd = 3.5, lty = 2, add = TRUE )
box()
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_3, "_", 
	scenario, "_without a legend" ), path = MapsDir, type = "png" )
pause( 2 ); graphics.off()


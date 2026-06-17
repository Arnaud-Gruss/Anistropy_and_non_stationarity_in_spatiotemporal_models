###############################################################################################################################
##
##  (01) This script explores the issues of barriers and spatially stratified anisotropy for the North Pacific case
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
library( gridExtra )
library( surveyjoin )

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
MapsDir <- paste0( DIR$Maps, "North Pacific case" )
if ( !dir.exists( MapsDir ) ) { dir.create( MapsDir ) } 

######## Define some parameters/settings
case <- "North Pacific" 
location_1 <- "Gulf of Alaska"
location_2 <- "Bering Sea"
map_sf <- st_union( ne_countries( country = c( "united states of america", "canada" ) , scale = 10 ) )
location_of_interest_1 <- c( -151.756347, 58.820505 )
location_of_interest_2 <- c( -170.960448, 60.550198 )
Mapping_parameters <- list()
Mapping_parameters[["xlim_map"]] <- c( -178, -128 )
Mapping_parameters[["ylim_map"]] <- c( 42, 68.37320 )
Mapping_parameters[["land_color"]] <- makeTransparent( "gray88", alpha = 170 )
Mapping_parameters[["legend_x"]] <- c( 0.86, 0.91 )
Mapping_parameters[["legend_y"]] <- c( 0.72, 0.97 )
cutoff <- 0.5
barrier_stiffness <- 1e-16
cellsize <- 0.5 * c( 1, 1 )

######## Load and process data
d <- get_data( "sablefish" )
Data <- as.data.frame( d )[,c( 'survey_name', 'lat_start', 'lon_start', 'effort', 'catch_weight' )]
Data <- na.omit( Data )
Data$Lon <- Data$lon_start
Data$Lat <- Data$lat_start
Survey_data <- Data
Data <- expand.grid( "Lon" = seq( Min( Data$Lon ), Max( Data$Lon ), length.out  = 200 ), 
	"Lat" = seq( Min( Data$Lat ), Max( Data$Lat ), length.out  = 200 ) )
Data$catch_weight <- 0
Data$catch_weight[seq( 2, length( Data$catch_weight ), by = 2 )] <- rlnorm( nrow( Data ) / 2, meanlog = 0, sdlog = 1 )
Data$var <- "density"

######## Generate model settings  
formula <- catch_weight ~ 1
space_term <- "
"
spacetime_term <- "
" 
time_term <- "
" 

######## Generate a mesh for the tinyVAST model using the SPDE method and assuming isotropy
mesh <- fm_mesh_2d( Data[,c( 'Lon', 'Lat' )], cutoff = cutoff )

plot( mesh )
plot( map_sf, add = TRUE )
pause( 3 ); graphics.off()

######## Prepare the tinyVAST model using the SPDE method and assuming isotropy 
fit_norun <- tinyVAST( 
	data = droplevels( Data ), 
	formula = formula, 
	space_term = space_term,
  	spacetime_term = spacetime_term,
  	time_term = time_term,
  	spatial_domain = mesh,
  	family = tweedie( "log" ),
	space_columns = c( "Lon", "Lat" ),
  	variable_column = "var",
	control = tinyVASTcontrol( 
		verbose = TRUE, 
		opt_loops = 0,
		newton_loops = 0, 
		trace = 1,
    		use_anisotropy = FALSE, 
		getsd = FALSE,
    		calculate = FALSE, 
	)
)
tmb_par <- fit_norun$tmb_inputs$tmb_par

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
  	variable_column = "var",
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

######## Approximate spatial correlation for the tinyVAST model using the SPDE method and assuming isotropy
######## for the first location of interest
fit <- SPDE_Isotropy
mesh_sf <- fm_as_sfc( mesh )
domain_sf <- st_union( mesh_sf )
grid_sf <- st_make_grid( domain_sf, cellsize = cellsize )
grid_sf <- st_intersection( grid_sf, domain_sf )
pred <- st_coordinates( st_centroid( grid_sf ) )
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
fit <- SPDE_Isotropy
mesh_sf <- fm_as_sfc( mesh )
domain_sf <- st_union( mesh_sf )
grid_sf <- st_make_grid( domain_sf, cellsize = cellsize )
grid_sf <- st_intersection( grid_sf, domain_sf )
pred <- st_coordinates( st_centroid( grid_sf ) )
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

######## Generate a mesh for a tinyVAST model using the SPDE method and covariate-based regional anisotropy 
mesh <- fm_mesh_2d( Data[,c( "Lon", "Lat" )], cutoff = cutoff )
nn <- RANN::nn2( data = Survey_data[,c( 'Lon', 'Lat' )], query = mesh$loc[,1:2], k = 1 )
survey_s <- Survey_data[nn$nn.idx, 'survey_name']
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
  	x_goa = ifelse( survey_s %in% c( "Gulf of Alaska", "Aleutian Islands" ), mesh$loc[,1], NA ),
  	y_goa = ifelse( survey_s %in% c( "Gulf of Alaska", "Aleutian Islands" ), mesh$loc[,2], NA ),
  	x_bs = ifelse( survey_s %in% c( "Bering Sea Slope", "eastern Bering Sea", "northern Bering Sea" ), 
		mesh$loc[,1], NA ),
  	y_bs = ifelse( survey_s %in% c( "Bering Sea Slope", "eastern Bering Sea", "northern Bering Sea"), 
		mesh$loc[,2], NA ),
  	x_wc = ifelse( survey_s %in% c( "Gulf of Alaska", "Aleutian Islands", "Bering Sea Slope", "eastern Bering Sea",
		"northern Bering Sea"), NA, mesh$loc[,1] ),
  	y_wc = ifelse( survey_s %in% c( "Gulf of Alaska", "Aleutian Islands", "Bering Sea Slope", "eastern Bering Sea",
		"northern Bering Sea" ), NA, mesh$loc[,2] )
)

plot( mesh_cov )
pause( 3 ); graphics.off()

######## Prepare the tinyVAST model using SPDE method and covariate-based regional anisotropy 
fit_norun <- tinyVAST( 
	data = droplevels( Data ), 
	formula = formula, 
	space_term = space_term,
  	spacetime_term = spacetime_term,
  	time_term = time_term,
  	spatial_domain = mesh_cov,
  	family = tweedie( "log" ),
	space_columns = c( "Lon", "Lat" ),
  	variable_column = "var",
  	development = list(
		vertex_formula = ~
			I(x_goa + y_goa) + y_goa +
      			I(x_bs + y_bs) + y_bs +
			0
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
tmb_par <- fit_norun$tmb_inputs$tmb_par
tmb_par$ln_H_input[3] <- 0
tmb_par$ln_H_input[4] <- 50
tmb_par$ln_H_input[5] <- 50
tmb_par$ln_H_input[6] <- 50

######## Fit the tinyVAST model using the SPDE method and covariate-based regional anisotropy 
Regional_Anisotropy <- tinyVAST( 
	data = droplevels( Data ), 
	formula = formula, 
	space_term = space_term,
  	spacetime_term = spacetime_term,
  	time_term = time_term,
  	spatial_domain = mesh_cov,
  	family = tweedie( "log" ),
	space_columns = c( "Lon", "Lat" ),
  	variable_column = "var",
  	development = list(
		vertex_formula = ~
			I(x_goa + y_goa) + y_goa +
      			I(x_bs + y_bs) + y_bs +
			0
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
######## for the first location of interest
fit <- Regional_Anisotropy
mesh_sf <- fm_as_sfc( mesh )
domain_sf <- st_union( mesh_sf )
grid_sf <- st_make_grid( domain_sf, cellsize = cellsize )
grid_sf <- st_intersection( grid_sf, domain_sf )
pred <- st_coordinates( st_centroid( grid_sf ) )
r1 <- spatial_cor(
	fit$rep$Q,
  	mesh = mesh,
  	coord = location_of_interest_1,
  	pred = pred
)
plot_sf <- st_sf( st_geometry( grid_sf ), r1 = r1 )
r_Regional_Anisotropy_1 <- plot_sf

plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
pause( 3 ); graphics.off()

######## Approximate spatial correlation for the tinyVAST model using the SPDE method and covariate-based regional anisotropy 
######## for the second location of interest
fit <- Regional_Anisotropy
mesh_sf <- fm_as_sfc( mesh )
domain_sf <- st_union( mesh_sf )
grid_sf <- st_make_grid( domain_sf, cellsize = cellsize )
grid_sf <- st_intersection( grid_sf, domain_sf )
pred <- st_coordinates( st_centroid( grid_sf ) )
r1 <- spatial_cor(
	fit$rep$Q,
  	mesh = mesh,
  	coord = location_of_interest_2,
  	pred = pred
)
plot_sf <- st_sf( st_geometry( grid_sf ), r1 = r1 )
r_Regional_Anisotropy_2 <- plot_sf

plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
pause( 3 ); graphics.off()

######## Generate a mesh for a tinyVAST model using the SPDE-Barrier method and covariate-based regional anisotropy 
mesh <- fm_mesh_2d( Data[,c( "Lon", "Lat" )], cutoff = cutoff )
nn <- RANN::nn2( data = Survey_data[,c( 'Lon', 'Lat' )], query = mesh$loc[,1:2], k = 1 )
survey_s <- Survey_data[nn$nn.idx, 'survey_name']
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
  	x_goa = ifelse( survey_s %in% c( "Gulf of Alaska", "Aleutian Islands" ), mesh$loc[,1], NA ),
  	y_goa = ifelse( survey_s %in% c( "Gulf of Alaska", "Aleutian Islands" ), mesh$loc[,2], NA ),
  	x_bs = ifelse( survey_s %in% c( "Bering Sea Slope", "eastern Bering Sea", "northern Bering Sea" ), 
		mesh$loc[,1], NA ),
  	y_bs = ifelse( survey_s %in% c( "Bering Sea Slope", "eastern Bering Sea", "northern Bering Sea"), 
		mesh$loc[,2], NA ),
  	x_wc = ifelse( survey_s %in% c( "Gulf of Alaska", "Aleutian Islands", "Bering Sea Slope", "eastern Bering Sea",
		"northern Bering Sea"), NA, mesh$loc[,1] ),
  	y_wc = ifelse( survey_s %in% c( "Gulf of Alaska", "Aleutian Islands", "Bering Sea Slope", "eastern Bering Sea",
		"northern Bering Sea" ), NA, mesh$loc[,2] )
)
intersected <- st_intersects( st_make_valid( st_set_crs( fm_as_sfc( mesh ), st_crs( map_sf ) ) ), st_make_valid( map_sf ) )
mesh_cov$triangle_covariates$barrier_polygon <- ifelse( is.na( as.numeric( intersected ) ), 0, 1 )

plot( mesh_cov )
pause( 3 ); graphics.off()

map_data_sp <- as( map_sf, 'Spatial' )
n_triangles <- length( mesh_cov$graph$tv[,1] )
posTri <- matrix( NA, nrow = n_triangles, ncol = 2 )
for ( tri_index in 1 : n_triangles ) {
	temp <- mesh_cov$loc[mesh_cov$graph$tv[tri_index,], ]
	posTri[tri_index,] <- colMeans( temp )[c( 1, 2 )]
}
posTri <- sp::SpatialPoints( posTri, proj4string = map_data_sp@proj4string )
mesh_triangles_over_land <- unlist( sp::over( map_data_sp, posTri, returnList = TRUE ) )
posTri_plot <- as.data.frame( posTri@coords[,1 : 2] )
plot( posTri_plot[,1], posTri_plot[,2], 
	type = "n", ylab = "Latitude (°N)", xlab = "Longitude (°E)", cex.lab = 1.25, font.lab = 2 )
sp::plot( map_data_sp, col = "grey88", add = T )
points( posTri_plot[mesh_cov$triangle_covariates$barrier_polygon == 1,1], 
	posTri_plot[mesh_cov$triangle_covariates$barrier_polygon == 1,2], 
	pch = 16, col = "red", cex = 1.2 )
points( posTri_plot[!mesh_cov$triangle_covariates$barrier_polygon == 1,1], 
	posTri_plot[!mesh_cov$triangle_covariates$barrier_polygon ==1,2], 
	pch = 16, col = "green", cex = 1.2 )
sp::plot( mesh_cov, add = T )
pause( 3 ); graphics.off()

######## Prepare the tinyVAST model using SPDE-Barrier method and covariate-based regional anisotropy 
fit_norun <- tinyVAST( 
	data = droplevels( Data ), 
	formula = formula, 
	space_term = space_term,
  	spacetime_term = spacetime_term,
  	time_term = time_term,
  	spatial_domain = mesh_cov,
  	family = tweedie( "log" ),
	space_columns = c( "Lon", "Lat" ),
  	variable_column = "var",
  	development = list(
		triangle_formula = ~ offset(barrier_polygon),
		vertex_formula = ~
			I(x_goa + y_goa) + y_goa +
      			I(x_bs + y_bs) + y_bs +
			0
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
tmb_par <- fit_norun$tmb_inputs$tmb_par
tmb_par$ln_H_input[3] <- 0
tmb_par$ln_H_input[4] <- 50
tmb_par$ln_H_input[5] <- 50
tmb_par$ln_H_input[6] <- 50

######## Fit the tinyVAST model using the SPDE-Barrier method and covariate-based regional anisotropy 
Regional_Anisotropy_Barrier <- tinyVAST( 
	data = droplevels( Data ), 
	formula = formula, 
	space_term = space_term,
  	spacetime_term = spacetime_term,
  	time_term = time_term,
  	spatial_domain = mesh_cov,
  	family = tweedie( "log" ),
	space_columns = c( "Lon", "Lat" ),
  	variable_column = "var",
  	development = list(
		triangle_formula = ~ offset(barrier_polygon),
		vertex_formula = ~
			I(x_goa + y_goa) + y_goa +
      			I(x_bs + y_bs) + y_bs +
			0
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

######## Approximate spatial correlation for the tinyVAST model using the SPDE-Barrier method and 
######## covariate-based regional anisotropy for the first location of interest
fit <- Regional_Anisotropy_Barrier
mesh_sf <- fm_as_sfc( mesh )
domain_sf <- st_union( mesh_sf )
grid_sf <- st_make_grid( domain_sf, cellsize = cellsize )
grid_sf <- st_intersection( grid_sf, domain_sf )
pred <- st_coordinates( st_centroid( grid_sf ) )
r1 <- spatial_cor(
	fit$rep$Q,
  	mesh = mesh,
  	coord = location_of_interest_1,
  	pred = pred
)
plot_sf <- st_sf( st_geometry( grid_sf ), r1 = r1 )
r_Regional_Anisotropy_Barrier_1 <- plot_sf

plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
pause( 3 ); graphics.off()

######## Approximate spatial correlation for the tinyVAST model using the SPDE-Barrier method and 
######## covariate-based regional anisotropy for the second location of interest
fit <- Regional_Anisotropy_Barrier
mesh_sf <- fm_as_sfc( mesh )
domain_sf <- st_union( mesh_sf )
grid_sf <- st_make_grid( domain_sf, cellsize = cellsize )
grid_sf <- st_intersection( grid_sf, domain_sf )
pred <- st_coordinates( st_centroid( grid_sf ) )
r1 <- spatial_cor(
	fit$rep$Q,
  	mesh = mesh,
  	coord = location_of_interest_2,
  	pred = pred
)
plot_sf <- st_sf( st_geometry( grid_sf ), r1 = r1 )
r_Regional_Anisotropy_Barrier_2 <- plot_sf

plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
pause( 3 ); graphics.off()

######## Produce a map showing the approximated spatial correlation for the first location of interest estimated by the
######## tinyVAST model using the SPDE method and assuming isotropy and save it 
scenario <- "SPDE_Isotropy"
plot_sf <- r_SPDE_Isotropy_1
zlim <- range( c( r_SPDE_Isotropy_1$r1, r_SPDE_Isotropy_2$r1, 
	r_Regional_Anisotropy_1$r1, r_Regional_Anisotropy_2$r1,
	r_Regional_Anisotropy_Barrier_1$r1, r_Regional_Anisotropy_Barrier_2$r1 ), na.rm = T )

#### Produce a version of the map with a legend and save it
new.graph( 0.5 )
par( mar = c( 4.5, 4.5, 1, 1.5 ), cex = 1.3 )
plot( Mapping_parameters[["xlim_map"]], Mapping_parameters[["ylim_map"]], 
	xlim = Mapping_parameters[["xlim_map"]], ylim = Mapping_parameters[["ylim_map"]], 
	type = "n", ylab = "Latitude (°N)", xlab = "Longitude (°E)", cex.lab = 1.1, font.lab = 2 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, add = TRUE ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
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
	cex = 1.2, align = align, gradient = gradient, col = "black" )
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_1, "_", 
	scenario, "_with a legend" ), path = MapsDir, type = "png" )
pause( 3 ); graphics.off()

#### Produce a version of the map without a legend and save it
new.graph( 0.5 )
par( mar = c( 4.5, 4.5, 1, 1.5 ), cex = 1.3 )
plot( Mapping_parameters[["xlim_map"]], Mapping_parameters[["ylim_map"]], 
	xlim = Mapping_parameters[["xlim_map"]], ylim = Mapping_parameters[["ylim_map"]], 
	type = "n", ylab = "Latitude (°N)", xlab = "Longitude (°E)", cex.lab = 1.1, font.lab = 2 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, add = TRUE ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
box()
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_1, "_", 
	scenario, "_without a legend" ), path = MapsDir, type = "png" )
pause( 3 ); graphics.off()

######## Produce a map showing the approximated spatial correlation for the second location of interest estimated by the
######## tinyVAST model using the SPDE method and assuming isotropy and save it 
scenario <- "SPDE_Isotropy"
plot_sf <- r_SPDE_Isotropy_2
zlim <- range( c( r_SPDE_Isotropy_1$r1, r_SPDE_Isotropy_2$r1, 
	r_Regional_Anisotropy_1$r1, r_Regional_Anisotropy_2$r1,
	r_Regional_Anisotropy_Barrier_1$r1, r_Regional_Anisotropy_Barrier_2$r1 ), na.rm = T )

#### Produce a version of the map with a legend and save it
new.graph( 0.5 )
par( mar = c( 4.5, 4.5, 1, 1.5 ), cex = 1.3 )
plot( Mapping_parameters[["xlim_map"]], Mapping_parameters[["ylim_map"]], 
	xlim = Mapping_parameters[["xlim_map"]], ylim = Mapping_parameters[["ylim_map"]], 
	type = "n", ylab = "Latitude (°N)", xlab = "Longitude (°E)", cex.lab = 1.1, font.lab = 2 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, add = TRUE ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
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
	cex = 1.2, align = align, gradient = gradient, col = "black" )
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_2, "_", 
	scenario, "_with a legend" ), path = MapsDir, type = "png" )
pause( 3 ); graphics.off()

#### Produce a version of the map without a legend and save it
new.graph( 0.5 )
par( mar = c( 4.5, 4.5, 1, 1.5 ), cex = 1.3 )
plot( Mapping_parameters[["xlim_map"]], Mapping_parameters[["ylim_map"]], 
	xlim = Mapping_parameters[["xlim_map"]], ylim = Mapping_parameters[["ylim_map"]], 
	type = "n", ylab = "Latitude (°N)", xlab = "Longitude (°E)", cex.lab = 1.1, font.lab = 2 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, add = TRUE ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
box()
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_2, "_", 
	scenario, "_without a legend" ), path = MapsDir, type = "png" )
pause( 3 ); graphics.off()

######## Produce a map showing the approximated spatial correlation for the first location of interest estimated by the
######## tinyVAST model using the SPDE method and covariate-based regional anisotropy and save it 
scenario <- "Regional_Anisotropy"
plot_sf <- r_Regional_Anisotropy_1
zlim <- range( c( r_SPDE_Isotropy_1$r1, r_SPDE_Isotropy_2$r1, 
	r_Regional_Anisotropy_1$r1, r_Regional_Anisotropy_2$r1,
	r_Regional_Anisotropy_Barrier_1$r1, r_Regional_Anisotropy_Barrier_2$r1 ), na.rm = T )

#### Produce a version of the map with a legend and save it
new.graph( 0.5 )
par( mar = c( 4.5, 4.5, 1, 1.5 ), cex = 1.3 )
plot( Mapping_parameters[["xlim_map"]], Mapping_parameters[["ylim_map"]], 
	xlim = Mapping_parameters[["xlim_map"]], ylim = Mapping_parameters[["ylim_map"]], 
	type = "n", ylab = "Latitude (°N)", xlab = "Longitude (°E)", cex.lab = 1.1, font.lab = 2 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, add = TRUE ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
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
	cex = 1.2, align = align, gradient = gradient, col = "black" )
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_1, "_", 
	scenario, "_with a legend" ), path = MapsDir, type = "png" )
pause( 3 ); graphics.off()

#### Produce a version of the map without a legend and save it
new.graph( 0.5 )
par( mar = c( 4.5, 4.5, 1, 1.5 ), cex = 1.3 )
plot( Mapping_parameters[["xlim_map"]], Mapping_parameters[["ylim_map"]], 
	xlim = Mapping_parameters[["xlim_map"]], ylim = Mapping_parameters[["ylim_map"]], 
	type = "n", ylab = "Latitude (°N)", xlab = "Longitude (°E)", cex.lab = 1.1, font.lab = 2 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, add = TRUE ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
box()
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_1, "_", 
	scenario, "_without a legend" ), path = MapsDir, type = "png" )
pause( 3 ); graphics.off()

######## Produce a map showing the approximated spatial correlation for the second location of interest estimated by the
######## tinyVAST model using the SPDE method and assuming isotropy and save it 
scenario <- "Regional_Anisotropy"
plot_sf <- r_Regional_Anisotropy_2
zlim <- range( c( r_SPDE_Isotropy_1$r1, r_SPDE_Isotropy_2$r1, 
	r_Regional_Anisotropy_1$r1, r_Regional_Anisotropy_2$r1,
	r_Regional_Anisotropy_Barrier_1$r1, r_Regional_Anisotropy_Barrier_2$r1 ), na.rm = T )

#### Produce a version of the map with a legend and save it
new.graph( 0.5 )
par( mar = c( 4.5, 4.5, 1, 1.5 ), cex = 1.3 )
plot( Mapping_parameters[["xlim_map"]], Mapping_parameters[["ylim_map"]], 
	xlim = Mapping_parameters[["xlim_map"]], ylim = Mapping_parameters[["ylim_map"]], 
	type = "n", ylab = "Latitude (°N)", xlab = "Longitude (°E)", cex.lab = 1.1, font.lab = 2 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, add = TRUE ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
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
	cex = 1.2, align = align, gradient = gradient, col = "black" )
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_2, "_", 
	scenario, "_with a legend" ), path = MapsDir, type = "png" )
pause( 3 ); graphics.off()

#### Produce a version of the map without a legend and save it
new.graph( 0.5 )
par( mar = c( 4.5, 4.5, 1, 1.5 ), cex = 1.3 )
plot( Mapping_parameters[["xlim_map"]], Mapping_parameters[["ylim_map"]], 
	xlim = Mapping_parameters[["xlim_map"]], ylim = Mapping_parameters[["ylim_map"]], 
	type = "n", ylab = "Latitude (°N)", xlab = "Longitude (°E)", cex.lab = 1.1, font.lab = 2 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, add = TRUE ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
box()
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_2, "_", 
	scenario, "_without a legend" ), path = MapsDir, type = "png" )
pause( 3 ); graphics.off()

######## Produce a map showing the approximated spatial correlation for the first location of interest estimated by the
######## tinyVAST model using the SPDE-Barrier method and covariate-based regional anisotropy and save it 
scenario <- "Regional_Anisotropy_Barrier"
plot_sf <- r_Regional_Anisotropy_Barrier_1
zlim <- range( c( r_SPDE_Isotropy_1$r1, r_SPDE_Isotropy_2$r1, 
	r_Regional_Anisotropy_1$r1, r_Regional_Anisotropy_2$r1,
	r_Regional_Anisotropy_Barrier_1$r1, r_Regional_Anisotropy_Barrier_2$r1 ), na.rm = T )

#### Produce a version of the map with a legend and save it
new.graph( 0.5 )
par( mar = c( 4.5, 4.5, 1, 1.5 ), cex = 1.3 )
plot( Mapping_parameters[["xlim_map"]], Mapping_parameters[["ylim_map"]], 
	xlim = Mapping_parameters[["xlim_map"]], ylim = Mapping_parameters[["ylim_map"]], 
	type = "n", ylab = "Latitude (°N)", xlab = "Longitude (°E)", cex.lab = 1.1, font.lab = 2 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, add = TRUE ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
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
	cex = 1.2, align = align, gradient = gradient, col = "black" )
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_1, "_", 
	scenario, "_with a legend" ), path = MapsDir, type = "png" )
pause( 3 ); graphics.off()

#### Produce a version of the map without a legend and save it
new.graph( 0.5 )
par( mar = c( 4.5, 4.5, 1, 1.5 ), cex = 1.3 )
plot( Mapping_parameters[["xlim_map"]], Mapping_parameters[["ylim_map"]], 
	xlim = Mapping_parameters[["xlim_map"]], ylim = Mapping_parameters[["ylim_map"]], 
	type = "n", ylab = "Latitude (°N)", xlab = "Longitude (°E)", cex.lab = 1.1, font.lab = 2 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, add = TRUE ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
box()
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_1, "_", 
	scenario, "_without a legend" ), path = MapsDir, type = "png" )
pause( 3 ); graphics.off()

######## Produce a map showing the approximated spatial correlation for the second location of interest estimated by the
######## tinyVAST model using the SPDE-Barrier method and assuming isotropy and save it 
scenario <- "Regional_Anisotropy_Barrier"
plot_sf <- r_Regional_Anisotropy_Barrier_2
zlim <- range( c( r_SPDE_Isotropy_1$r1, r_SPDE_Isotropy_2$r1, 
	r_Regional_Anisotropy_1$r1, r_Regional_Anisotropy_2$r1,
	r_Regional_Anisotropy_Barrier_1$r1, r_Regional_Anisotropy_Barrier_2$r1 ), na.rm = T )

#### Produce a version of the map with a legend and save it
new.graph( 0.5 )
par( mar = c( 4.5, 4.5, 1, 1.5 ), cex = 1.3 )
plot( Mapping_parameters[["xlim_map"]], Mapping_parameters[["ylim_map"]], 
	xlim = Mapping_parameters[["xlim_map"]], ylim = Mapping_parameters[["ylim_map"]], 
	type = "n", ylab = "Latitude (°N)", xlab = "Longitude (°E)", cex.lab = 1.1, font.lab = 2 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, add = TRUE ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
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
	cex = 1.2, align = align, gradient = gradient, col = "black" )
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_2, "_", 
	scenario, "_with a legend" ), path = MapsDir, type = "png" )
pause( 3 ); graphics.off()

#### Produce a version of the map without a legend and save it
new.graph( 0.5 )
par( mar = c( 4.5, 4.5, 1, 1.5 ), cex = 1.3 )
plot( Mapping_parameters[["xlim_map"]], Mapping_parameters[["ylim_map"]], 
	xlim = Mapping_parameters[["xlim_map"]], ylim = Mapping_parameters[["ylim_map"]], 
	type = "n", ylab = "Latitude (°N)", xlab = "Longitude (°E)", cex.lab = 1.1, font.lab = 2 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, add = TRUE ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
box()
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_2, "_", 
	scenario, "_without a legend" ), path = MapsDir, type = "png" )
pause( 3 ); graphics.off()


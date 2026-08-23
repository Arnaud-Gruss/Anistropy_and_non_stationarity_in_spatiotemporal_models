###############################################################################################################################
##
##  (04) Demonstration 4
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
library( Kendall )

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
MapsDir <- paste0( DIR$Maps, "West Coast South Island case" )
if ( !dir.exists( MapsDir ) ) { dir.create( MapsDir ) } 

######## Load mapping parameters 
load( file = make.filename( "Mapping_parameters_SNA7.RData", DIR$Input ) )

######## Define some parameters/settings and modify others
case <- "West Coast South Island" 
location_1 <- "Tasman Bay and Golden Bay"
location_2 <- "West Coast"
map_sf <- rnaturalearth::ne_countries( country = c( "new zealand" ), scale = 10 )
location_of_interest_1 <- c( 172.97, -40.67 )
location_of_interest_2 <- c( 171, -42.3 )
cutoff <- 0.2
barrier_stiffness <- 0.001
cellsize <- 0.025 * c( 1, 1 )
Mapping_parameters[["land_color"]] <- makeTransparent( "gray88", alpha = 170 )
Mapping_parameters[["xaxis"]] <- seq( from = 168, to = 175, by = 1 )
Mapping_parameters[["xaxis_label"]] <- as.character( Mapping_parameters[["xaxis"]] ) 
Mapping_parameters[["yaxis"]] <- seq( from = -44, to = -40, by = 1 )
Mapping_parameters[["yaxis_label"]] <- as.character( Mapping_parameters[["yaxis"]] ) 
Mapping_parameters[["legend_x"]] <- c( 0.06, 0.11 )
Mapping_parameters[["legend_y"]] <- c( 0.62, 0.87 )

######## Load and process data 
load( file = make.filename( "SNA_groomed_research_survey_data.RData", DIR$Input ) )
Data <- TRAWL_data
Data$Year <- Data$fyear 
Data <- Data[,c( 'Year', 'Lon', 'Lat', 'Catch_rate' )]
Data <- subset( Data, Year >= 2007 ) 
Data <- as.data.frame( Data )
Data$var <- "density"
Mesh_data <- expand.grid( "Lon" = seq( Min( Data$Lon ), Max( Data$Lon ), length.out  = 25 ), 
	"Lat" = seq( Min( Data$Lat ), Max( Data$Lat ), length.out  = 25 ) )
Mesh_data <- Data 

######## Define settings for the tinyVAST model
formula <- Catch_rate ~ 1
space_term <- "
	density <-> density, sd_space
"
spacetime_term <- "
	density -> density, 1, ar_st
  	density <-> density, 0, sd_st
" 
time_term <- "
	density -> density, 1, ar_t
  	density <-> density, 0, sd_t
" 

######## Generate a mesh for the tinyVAST model using the SPDE method and assuming isotropy
mesh <- fm_mesh_2d( Mesh_data[,c( 'Lon', 'Lat' )], cutoff = cutoff )

plot( mesh )
plot( map_sf, add = TRUE )
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
		newton_loops = 2, 
		trace = 1,
    		use_anisotropy = FALSE
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
mesh <- fm_mesh_2d( Mesh_data[,c( "Lon", "Lat" )], cutoff = cutoff )
nn <- RANN::nn2( data = Mesh_data[,c( 'Lon', 'Lat' )], query = mesh$loc[,1 : 2], k = 1 )
Mesh_data$Region <- ifelse( Mesh_data$Lon > 172.6, "TBGB", "WEST_COAST" )
region_s <- Mesh_data[nn$nn.idx, 'Region']
mesh_cov <- add_mesh_covariates(
  	mesh = mesh,
  	data = Mesh_data,
  	covariates = c(),
  	coords = c( "Lon", "Lat" )
)
mesh_cov$vertex_covariates <- data.frame(
  	x = ( mesh$loc[,1] ),
  	y = ( mesh$loc[,2] )
)
mesh_cov$vertex_covariates <- cbind(
  	mesh_cov$vertex_covariates,
  	x_tbgb = ifelse( region_s %in% c( "TBGB" ), mesh$loc[,1], NA ),
  	y_tbgb = ifelse( region_s %in% c( "TBGB" ), mesh$loc[,2], NA ),
  	x_wc = ifelse( region_s %in% c( "WEST_COAST" ), mesh$loc[,1], NA ),
  	y_wc = ifelse( region_s %in% c( "WEST_COAST" ), mesh$loc[,2], NA )
)

plot( mesh_cov )
pause( 3 ); graphics.off()

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
  	time_column = "Year",
  	times = min( Data$Year ) : max( Data$Year ),
  	variable_column = "var",
	development = list(
		vertex_formula = ~
			I(x_tbgb + y_tbgb) + y_tbgb +
				I(x_wc + y_wc) + y_wc +
				0 
  	),
	control = tinyVASTcontrol( 
		verbose = TRUE, 
		newton_loops = 2, 
		trace = 1,
    		use_anisotropy = FALSE
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
mesh <- fm_mesh_2d( Mesh_data[,c( "Lon", "Lat" )], cutoff = cutoff )
nn <- RANN::nn2( data = Mesh_data[,c( 'Lon', 'Lat' )], query = mesh$loc[,1 : 2], k = 1 )
Mesh_data$Region <- ifelse( Mesh_data$Lon > 172.6, "TBGB", "WEST_COAST" )
region_s <- Mesh_data[nn$nn.idx, 'Region']
mesh_cov <- add_mesh_covariates(
  	mesh = mesh,
  	data = Mesh_data,
  	covariates = c(),
  	coords = c( "Lon", "Lat" )
)
mesh_cov$vertex_covariates <- data.frame(
  	x = ( mesh$loc[,1] ),
  	y = ( mesh$loc[,2] )
)
mesh_cov$vertex_covariates <- cbind(
  	mesh_cov$vertex_covariates,
  	x_tbgb = ifelse( region_s %in% c( "TBGB" ), mesh$loc[,1], NA ),
  	y_tbgb = ifelse( region_s %in% c( "TBGB" ), mesh$loc[,2], NA ),
  	x_wc = ifelse( region_s %in% c( "WEST_COAST" ), mesh$loc[,1], NA ),
  	y_wc = ifelse( region_s %in% c( "WEST_COAST" ), mesh$loc[,2], NA )
)
map.sf <- ne_countries( country = c( "new zealand", scale = 10 ) )
intersected <- st_intersects( st_make_valid( st_set_crs( fm_as_sfc( mesh ), 
	st_crs( map.sf ) ) ), st_make_valid( map.sf ) )
mesh_cov$triangle_covariates$barrier_polygon <- ifelse( rowSums( as.matrix.sgbp( intersected ) ) == 0, 0, 1 )

plot( mesh_cov )
pause( 3 ); graphics.off()

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
  	time_column = "Year",
  	times = min( Data$Year ) : max( Data$Year ),
  	variable_column = "var",
	development = list(
		triangle_formula = ~ offset(barrier_polygon),
		vertex_formula = ~
			I(x_tbgb + y_tbgb) + y_tbgb +
				I(x_wc + y_wc) + y_wc +
				0 
  	),
	control = tinyVASTcontrol( 
		verbose = TRUE, 
		newton_loops = 2, 
		trace = 1,
    		use_anisotropy = FALSE, 
		barrier_stiffness = barrier_stiffness
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
par( oma = c( 3, 3, 0, 0 ), cex = 1.1 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, main = "", axes = FALSE, 
	xlim = st_bbox(plot_sf)[c( 1, 3 )], ylim = st_bbox(plot_sf)[c( 2, 4 )] ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
box()
axis( side = 1, at = Mapping_parameters[["xaxis"]], labels = Mapping_parameters[["xaxis_label"]] )
axis( side = 2, at = Mapping_parameters[["yaxis"]], labels = Mapping_parameters[["yaxis_label"]] )
minor.tick( nx = 1, ny = 1 )
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
pause( 3 ); graphics.off()

#### Produce a version of the map with a legend and save it
new.graph( 0.5 )
par( oma = c( 3, 3, 0, 0 ), cex = 1.1 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, main = "", axes = FALSE, 
	xlim = st_bbox(plot_sf)[c( 1, 3 )], ylim = st_bbox(plot_sf)[c( 2, 4 )] ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
box()
axis( side = 1, at = Mapping_parameters[["xaxis"]], labels = Mapping_parameters[["xaxis_label"]] )
axis( side = 2, at = Mapping_parameters[["yaxis"]], labels = Mapping_parameters[["yaxis_label"]] )
minor.tick( nx = 1, ny = 1 )
mtext( side = 1, outer = TRUE, "Longitude (°E)", cex = 1.25, line = par()$oma[1]/2, padj = 0 )
mtext( side = 2, outer = TRUE, "Latitude (°N)", cex = 1.25, line = par()$oma[2]/2, padj = 0, las = 0 )
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
par( oma = c( 3, 3, 0, 0 ), cex = 1.1 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, main = "", axes = FALSE, 
	xlim = st_bbox(plot_sf)[c( 1, 3 )], ylim = st_bbox(plot_sf)[c( 2, 4 )] ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
box()
axis( side = 1, at = Mapping_parameters[["xaxis"]], labels = Mapping_parameters[["xaxis_label"]] )
axis( side = 2, at = Mapping_parameters[["yaxis"]], labels = Mapping_parameters[["yaxis_label"]] )
minor.tick( nx = 1, ny = 1 )
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
pause( 3 ); graphics.off()

#### Produce a version of the map with a legend and save it
new.graph( 0.5 )
par( oma = c( 3, 3, 0, 0 ), cex = 1.1 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, main = "", axes = FALSE, 
	xlim = st_bbox(plot_sf)[c( 1, 3 )], ylim = st_bbox(plot_sf)[c( 2, 4 )] ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
box()
axis( side = 1, at = Mapping_parameters[["xaxis"]], labels = Mapping_parameters[["xaxis_label"]] )
axis( side = 2, at = Mapping_parameters[["yaxis"]], labels = Mapping_parameters[["yaxis_label"]] )
minor.tick( nx = 1, ny = 1 )
mtext( side = 1, outer = TRUE, "Longitude (°E)", cex = 1.25, line = par()$oma[1]/2, padj = 0 )
mtext( side = 2, outer = TRUE, "Latitude (°N)", cex = 1.25, line = par()$oma[2]/2, padj = 0, las = 0 )
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
par( oma = c( 3, 3, 0, 0 ), cex = 1.1 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, main = "", axes = FALSE, 
	xlim = st_bbox(plot_sf)[c( 1, 3 )], ylim = st_bbox(plot_sf)[c( 2, 4 )] ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
box()
axis( side = 1, at = Mapping_parameters[["xaxis"]], labels = Mapping_parameters[["xaxis_label"]] )
axis( side = 2, at = Mapping_parameters[["yaxis"]], labels = Mapping_parameters[["yaxis_label"]] )
minor.tick( nx = 1, ny = 1 )
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
pause( 3 ); graphics.off()

#### Produce a version of the map with a legend and save it
new.graph( 0.5 )
par( oma = c( 3, 3, 0, 0 ), cex = 1.1 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, main = "", axes = FALSE, 
	xlim = st_bbox(plot_sf)[c( 1, 3 )], ylim = st_bbox(plot_sf)[c( 2, 4 )] ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
box()
axis( side = 1, at = Mapping_parameters[["xaxis"]], labels = Mapping_parameters[["xaxis_label"]] )
axis( side = 2, at = Mapping_parameters[["yaxis"]], labels = Mapping_parameters[["yaxis_label"]] )
minor.tick( nx = 1, ny = 1 )
mtext( side = 1, outer = TRUE, "Longitude (°E)", cex = 1.25, line = par()$oma[1]/2, padj = 0 )
mtext( side = 2, outer = TRUE, "Latitude (°N)", cex = 1.25, line = par()$oma[2]/2, padj = 0, las = 0 )
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
par( oma = c( 3, 3, 0, 0 ), cex = 1.1 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, main = "", axes = FALSE, 
	xlim = st_bbox(plot_sf)[c( 1, 3 )], ylim = st_bbox(plot_sf)[c( 2, 4 )] ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
box()
axis( side = 1, at = Mapping_parameters[["xaxis"]], labels = Mapping_parameters[["xaxis_label"]] )
axis( side = 2, at = Mapping_parameters[["yaxis"]], labels = Mapping_parameters[["yaxis_label"]] )
minor.tick( nx = 1, ny = 1 )
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
pause( 3 ); graphics.off()

#### Produce a version of the map with a legend and save it
new.graph( 0.5 )
par( oma = c( 3, 3, 0, 0 ), cex = 1.1 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, main = "", axes = FALSE, 
	xlim = st_bbox(plot_sf)[c( 1, 3 )], ylim = st_bbox(plot_sf)[c( 2, 4 )] ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
box()
axis( side = 1, at = Mapping_parameters[["xaxis"]], labels = Mapping_parameters[["xaxis_label"]] )
axis( side = 2, at = Mapping_parameters[["yaxis"]], labels = Mapping_parameters[["yaxis_label"]] )
minor.tick( nx = 1, ny = 1 )
mtext( side = 1, outer = TRUE, "Longitude (°E)", cex = 1.25, line = par()$oma[1]/2, padj = 0 )
mtext( side = 2, outer = TRUE, "Latitude (°N)", cex = 1.25, line = par()$oma[2]/2, padj = 0, las = 0 )
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
par( oma = c( 3, 3, 0, 0 ), cex = 1.1 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, main = "", axes = FALSE, 
	xlim = st_bbox(plot_sf)[c( 1, 3 )], ylim = st_bbox(plot_sf)[c( 2, 4 )] ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
box()
axis( side = 1, at = Mapping_parameters[["xaxis"]], labels = Mapping_parameters[["xaxis_label"]] )
axis( side = 2, at = Mapping_parameters[["yaxis"]], labels = Mapping_parameters[["yaxis_label"]] )
minor.tick( nx = 1, ny = 1 )
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
pause( 3 ); graphics.off()

#### Produce a version of the map with a legend and save it
new.graph( 0.5 )
par( oma = c( 3, 3, 0, 0 ), cex = 1.1 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, main = "", axes = FALSE, 
	xlim = st_bbox(plot_sf)[c( 1, 3 )], ylim = st_bbox(plot_sf)[c( 2, 4 )] ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
box()
axis( side = 1, at = Mapping_parameters[["xaxis"]], labels = Mapping_parameters[["xaxis_label"]] )
axis( side = 2, at = Mapping_parameters[["yaxis"]], labels = Mapping_parameters[["yaxis_label"]] )
minor.tick( nx = 1, ny = 1 )
mtext( side = 1, outer = TRUE, "Longitude (°E)", cex = 1.25, line = par()$oma[1]/2, padj = 0 )
mtext( side = 2, outer = TRUE, "Latitude (°N)", cex = 1.25, line = par()$oma[2]/2, padj = 0, las = 0 )
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_1, "_", 
	scenario, "_without a legend" ), path = MapsDir, type = "png" )
pause( 3 ); graphics.off()

######## Produce a map showing the approximated spatial correlation for the second location of interest estimated by the
######## tinyVAST model using the SPDE-Barrier method and covariate-based regional anisotropy and save it 
scenario <- "Regional_Anisotropy_Barrier"
plot_sf <- r_Regional_Anisotropy_Barrier_2
zlim <- range( c( r_SPDE_Isotropy_1$r1, r_SPDE_Isotropy_2$r1, 
	r_Regional_Anisotropy_1$r1, r_Regional_Anisotropy_2$r1,
	r_Regional_Anisotropy_Barrier_1$r1, r_Regional_Anisotropy_Barrier_2$r1 ), na.rm = T )

#### Produce a version of the map with a legend and save it
new.graph( 0.5 )
par( oma = c( 3, 3, 0, 0 ), cex = 1.1 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, main = "", axes = FALSE, 
	xlim = st_bbox(plot_sf)[c( 1, 3 )], ylim = st_bbox(plot_sf)[c( 2, 4 )] ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
box()
axis( side = 1, at = Mapping_parameters[["xaxis"]], labels = Mapping_parameters[["xaxis_label"]] )
axis( side = 2, at = Mapping_parameters[["yaxis"]], labels = Mapping_parameters[["yaxis_label"]] )
minor.tick( nx = 1, ny = 1 )
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
pause( 3 ); graphics.off()

#### Produce a version of the map with a legend and save it
new.graph( 0.5 )
par( oma = c( 3, 3, 0, 0 ), cex = 1.1 )
plot( plot_sf, border = NA, reset = FALSE, key.pos = NULL, main = "", axes = FALSE, 
	xlim = st_bbox(plot_sf)[c( 1, 3 )], ylim = st_bbox(plot_sf)[c( 2, 4 )] ) 
plot( map_sf, col = Mapping_parameters[["land_color"]], border = NA, add = TRUE )
box()
axis( side = 1, at = Mapping_parameters[["xaxis"]], labels = Mapping_parameters[["xaxis_label"]] )
axis( side = 2, at = Mapping_parameters[["yaxis"]], labels = Mapping_parameters[["yaxis_label"]] )
minor.tick( nx = 1, ny = 1 )
mtext( side = 1, outer = TRUE, "Longitude (°E)", cex = 1.25, line = par()$oma[1]/2, padj = 0 )
mtext( side = 2, outer = TRUE, "Latitude (°N)", cex = 1.25, line = par()$oma[2]/2, padj = 0, las = 0 )
SavePlot( filename = paste0( "Approximated spatial correlation for the ", location_2, "_", 
	scenario, "_without a legend" ), path = MapsDir, type = "png" )
pause( 3 ); graphics.off()

######## Save the estimated tinyVAST models in a .RData file 
save( SPDE_Isotropy, Regional_Anisotropy, Regional_Anisotropy_Barrier, 
	file = make.filename( "Estimated_tinyVAST_models_West_Coast_South_Island.RData", DIR$Output ) )

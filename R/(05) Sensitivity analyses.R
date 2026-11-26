###############################################################################################################################
##
##  (05) Sensitivity analyses
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
MapsDir <- paste0( DIR$Maps, "Sensitivity analyses" )
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
cutoff_1 <- 0.2
cutoff_2 <- 0.4
cutoff_3 <- 0.1
barrier_stiffness_1 <- 0.001
barrier_stiffness_2 <- 1e-16
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
	density -> density, 1, NA, 1
  	density <-> density, 0, sd_st
" 
time_term <- "
	density -> density, 1, NA, 1
  	density <-> density, 0, sd_t
" 

######## Generate a mesh for a tinyVAST model using the SPDE method and covariate-based regional anisotropy 
######## under Scenario 1 (cutoff_1 + barrier_stiffness_1) 
cutoff <- cutoff_1
barrier_stiffness <- barrier_stiffness_1
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
######## under Scenario 1 (cutoff_1 + barrier_stiffness_1) 
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

######## Examine the fitted tinyVAST model using the SPDE method and covariate-based regional anisotropy 
######## under Scenario 1 (cutoff_1 + barrier_stiffness_1) 
cAIC( Regional_Anisotropy_Scenario1 ) #### 3265.985

######## Generate a mesh for a tinyVAST model using the SPDE method and covariate-based regional anisotropy 
######## under Scenario 2 (cutoff_2 + barrier_stiffness_1) 
cutoff <- cutoff_2
barrier_stiffness <- barrier_stiffness_1
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
######## under Scenario 2 (cutoff_2 + barrier_stiffness_1) 
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

######## Examine the fitted tinyVAST model using the SPDE method and covariate-based regional anisotropy 
######## under Scenario 2 (cutoff_2 + barrier_stiffness_1) 
cAIC( Regional_Anisotropy_Scenario2 ) #### 3320.586

######## Generate a mesh for a tinyVAST model using the SPDE method and covariate-based regional anisotropy 
######## under Scenario 3 (cutoff_3 + barrier_stiffness_1) 
cutoff <- cutoff_3
barrier_stiffness <- barrier_stiffness_1
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
######## under Scenario 3 (cutoff_3 + barrier_stiffness_1) 
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

######## Examine the fitted tinyVAST model using the SPDE method and covariate-based regional anisotropy 
######## under Scenario 3 (cutoff_3 + barrier_stiffness_1) 
cAIC( Regional_Anisotropy_Scenario3 ) #### 3256.208

######## Generate a mesh for a tinyVAST model using the SPDE method and covariate-based regional anisotropy 
######## under Scenario 4 (cutoff_1 + barrier_stiffness_2) 
cutoff <- cutoff_1
barrier_stiffness <- barrier_stiffness_2
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
######## under Scenario 4 (cutoff_1 + barrier_stiffness_2) 
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

######## Examine the fitted tinyVAST model using the SPDE method and covariate-based regional anisotropy 
######## under Scenario 4 (cutoff_1 + barrier_stiffness_2) 
cAIC( Regional_Anisotropy_Scenario4 ) #### 3265.985

######## Generate a mesh for a tinyVAST model using the SPDE method and covariate-based regional anisotropy 
######## under Scenario 5 (cutoff_2 + barrier_stiffness_2) 
cutoff <- cutoff_2
barrier_stiffness <- barrier_stiffness_2
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
######## under Scenario 5 (cutoff_2 + barrier_stiffness_2) 
Regional_Anisotropy_Scenario5 <- tinyVAST( 
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

######## Examine the fitted tinyVAST model using the SPDE method and covariate-based regional anisotropy 
######## under Scenario 5 (cutoff_2 + barrier_stiffness_2) 
cAIC( Regional_Anisotropy_Scenario5 ) #### 3320.586

######## Generate a mesh for a tinyVAST model using the SPDE method and covariate-based regional anisotropy 
######## under Scenario 6 (cutoff_3 + barrier_stiffness_2) 
cutoff <- cutoff_3
barrier_stiffness <- barrier_stiffness_2
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
######## under Scenario 6 (cutoff_3 + barrier_stiffness_2) 
Regional_Anisotropy_Scenario6 <- tinyVAST( 
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

######## Examine the fitted tinyVAST model using the SPDE method and covariate-based regional anisotropy 
######## under Scenario 6 (cutoff_3 + barrier_stiffness_2) 
cAIC( Regional_Anisotropy_Scenario6 ) #### 3256.208

######## Save the fits of the tinyVAST models using the SPDE method and covariate-based regional anisotropy 
######## under Scenarios 1-6 
save( Regional_Anisotropy_Scenario1, Regional_Anisotropy_Scenario2, Regional_Anisotropy_Scenario3,
	Regional_Anisotropy_Scenario4, Regional_Anisotropy_Scenario5, Regional_Anisotropy_Scenario6, 
	file = make.filename( "Regional_Anisotropy_models_West_Coast_South_Island.RData", DIR$Output ) )

######## Generate predictions for each fishing year of the time period with the estimated tinyVAST models 
######## using the SPDE method and covariate-based regional anisotropy under Scenarios 1-6, 
######## and save these predictions in a .RData file
#### Load the estimated tinyVAST models and define prediction grids
load( file = make.filename( "Regional_Anisotropy_models_West_Coast_South_Island.RData", DIR$Output ) )
load( file = make.filename( "Prediction_grid_SNA7_5kmx5km.RData", DIR$Input ) )

#### Generate predictions for each fishing year of the time period with the estimated models
Years <- as.numeric( min( Data$Year ) : max( Data$Year ) )
areas_km2 <- 25
index_Regional_Anisotropy_Scenario1 <- index_Regional_Anisotropy_Scenario2 <- index_Regional_Anisotropy_Scenario3 <-
	index_Regional_Anisotropy_Scenario4 <- index_Regional_Anisotropy_Scenario5 <- 
	index_Regional_Anisotropy_Scenario6 <- NULL

pb <- txtProgressBar( min = 0, max = length( Years ), style = 3 )
yr <- 0

for( Year in Years ) {

	yr <- yr + 1
	setTxtProgressBar( pb, yr )

	## Produce predictive dataframes
	newdata <- data.frame( "Year" = Years[yr], 
		"Lon" = Prediction_grid[,'Lon'], 
		"Lat" = Prediction_grid[,'Lat'],
		"var" = "density" )

  	## Generate predictions 
	total_Regional_Anisotropy_Scenario1 <- integrate_output( Regional_Anisotropy_Scenario1, 
		newdata = newdata, area = areas_km2 )
	index_Regional_Anisotropy_Scenario1 <- cbind( index_Regional_Anisotropy_Scenario1, 
		total_Regional_Anisotropy_Scenario1[c( 'Est. (bias.correct)', 'Std. Error' )] )
	total_Regional_Anisotropy_Scenario2 <- integrate_output( Regional_Anisotropy_Scenario2, 
		newdata = newdata, area = areas_km2 )
	index_Regional_Anisotropy_Scenario2 <- cbind( index_Regional_Anisotropy_Scenario2, 
		total_Regional_Anisotropy_Scenario2[c( 'Est. (bias.correct)', 'Std. Error' )] )
	total_Regional_Anisotropy_Scenario3 <- integrate_output( Regional_Anisotropy_Scenario3, 
		newdata = newdata, area = areas_km2 )
	index_Regional_Anisotropy_Scenario3 <- cbind( index_Regional_Anisotropy_Scenario3, 
		total_Regional_Anisotropy_Scenario3[c( 'Est. (bias.correct)', 'Std. Error' )] )
	total_Regional_Anisotropy_Scenario4 <- integrate_output( Regional_Anisotropy_Scenario4, 
		newdata = newdata, area = areas_km2 )
	index_Regional_Anisotropy_Scenario4 <- cbind( index_Regional_Anisotropy_Scenario4, 
		total_Regional_Anisotropy_Scenario4[c( 'Est. (bias.correct)', 'Std. Error' )] )
	total_Regional_Anisotropy_Scenario5 <- integrate_output( Regional_Anisotropy_Scenario5, 
		newdata = newdata, area = areas_km2 )
	index_Regional_Anisotropy_Scenario5 <- cbind( index_Regional_Anisotropy_Scenario5, 
		total_Regional_Anisotropy_Scenario5[c( 'Est. (bias.correct)', 'Std. Error' )] )
	total_Regional_Anisotropy_Scenario6 <- integrate_output( Regional_Anisotropy_Scenario6, 
		newdata = newdata, area = areas_km2 )
	index_Regional_Anisotropy_Scenario6 <- cbind( index_Regional_Anisotropy_Scenario6, 
		total_Regional_Anisotropy_Scenario6[c( 'Est. (bias.correct)', 'Std. Error' )] )
	
}

close( pb )	

#### Save predictions in a .RData file 
save( Years, index_Regional_Anisotropy_Scenario1, index_Regional_Anisotropy_Scenario2, index_Regional_Anisotropy_Scenario3, 
	index_Regional_Anisotropy_Scenario4, index_Regional_Anisotropy_Scenario5, index_Regional_Anisotropy_Scenario6, 
	file = make.filename( "Predictions_Regional_Anisotropy_models_West_Coast_South_Island.RData", DIR$Output ) )

######## Produce figures showing the indices of relative biomass predicted by tinyVAST models 
######## using the SPDE method and covariate-based regional anisotropy under Scenarios 1-6, and save this figure
#### Load estimated models and their predictions
load( file = make.filename( "Regional_Anisotropy_models_West_Coast_South_Island.RData", DIR$Output ) )
load( file = make.filename( "Predictions_Regional_Anisotropy_models_West_Coast_South_Island.RData", DIR$Output ) )

#### Produce the figures and save them 
Biomasses_Regional_Anisotropy_Scenario1 <- data.frame( "year" = Years, Biomass = index_Regional_Anisotropy_Scenario1[1,], 
	SE.Biomass = index_Regional_Anisotropy_Scenario1[2,] )
Biomasses_Regional_Anisotropy_Scenario2 <- data.frame( "year" = Years, Biomass = index_Regional_Anisotropy_Scenario2[1,], 
	SE.Biomass = index_Regional_Anisotropy_Scenario2[2,] )
Biomasses_Regional_Anisotropy_Scenario3 <- data.frame( "year" = Years, Biomass = index_Regional_Anisotropy_Scenario3[1,], 
	SE.Biomass = index_Regional_Anisotropy_Scenario3[2,] )
Biomasses_Regional_Anisotropy_Scenario4 <- data.frame( "year" = Years, Biomass = index_Regional_Anisotropy_Scenario4[1,], 
	SE.Biomass = index_Regional_Anisotropy_Scenario4[2,] )
Biomasses_Regional_Anisotropy_Scenario5 <- data.frame( "year" = Years, Biomass = index_Regional_Anisotropy_Scenario5[1,], 
	SE.Biomass = index_Regional_Anisotropy_Scenario5[2,] )
Biomasses_Regional_Anisotropy_Scenario6 <- data.frame( "year" = Years, Biomass = index_Regional_Anisotropy_Scenario6[1,], 
	SE.Biomass = index_Regional_Anisotropy_Scenario6[2,] )
Mean_Biomass <- Mean( c( Biomasses_Regional_Anisotropy_Scenario1$Biomass, Biomasses_Regional_Anisotropy_Scenario2$Biomass, 
	Biomasses_Regional_Anisotropy_Scenario3$Biomass, Biomasses_Regional_Anisotropy_Scenario4$Biomass, 
	Biomasses_Regional_Anisotropy_Scenario5$Biomass, Biomasses_Regional_Anisotropy_Scenario6$Biomass ) ) 
Biomass_Regional_Anisotropy_Scenario1 <- Biomasses_Regional_Anisotropy_Scenario1$Biomass / Mean_Biomass 
Biomass_Regional_Anisotropy_Scenario2 <- Biomasses_Regional_Anisotropy_Scenario2$Biomass / Mean_Biomass 
Biomass_Regional_Anisotropy_Scenario3 <- Biomasses_Regional_Anisotropy_Scenario3$Biomass / Mean_Biomass 
Biomass_Regional_Anisotropy_Scenario4 <- Biomasses_Regional_Anisotropy_Scenario4$Biomass / Mean_Biomass 
Biomass_Regional_Anisotropy_Scenario5 <- Biomasses_Regional_Anisotropy_Scenario5$Biomass / Mean_Biomass 
Biomass_Regional_Anisotropy_Scenario6 <- Biomasses_Regional_Anisotropy_Scenario6$Biomass / Mean_Biomass 
SE.Biomass_Regional_Anisotropy_Scenario1 <- ( Biomasses_Regional_Anisotropy_Scenario1$Biomass %o% c( 1, 1 ) + 
	Biomasses_Regional_Anisotropy_Scenario1$SE.Biomass %o% c( -1, 1 ) ) / Mean_Biomass 
SE.Biomass_Regional_Anisotropy_Scenario2 <- ( Biomasses_Regional_Anisotropy_Scenario2$Biomass %o% c( 1, 1 ) + 
	Biomasses_Regional_Anisotropy_Scenario2$SE.Biomass %o% c( -1, 1 ) ) / Mean_Biomass 
SE.Biomass_Regional_Anisotropy_Scenario3 <- ( Biomasses_Regional_Anisotropy_Scenario3$Biomass %o% c( 1, 1 ) + 
	Biomasses_Regional_Anisotropy_Scenario3$SE.Biomass %o% c( -1, 1 ) ) / Mean_Biomass 
SE.Biomass_Regional_Anisotropy_Scenario4 <- ( Biomasses_Regional_Anisotropy_Scenario4$Biomass %o% c( 1, 1 ) + 
	Biomasses_Regional_Anisotropy_Scenario4$SE.Biomass %o% c( -1, 1 ) ) / Mean_Biomass 
SE.Biomass_Regional_Anisotropy_Scenario5 <- ( Biomasses_Regional_Anisotropy_Scenario5$Biomass %o% c( 1, 1 ) + 
	Biomasses_Regional_Anisotropy_Scenario5$SE.Biomass %o% c( -1, 1 ) ) / Mean_Biomass 
SE.Biomass_Regional_Anisotropy_Scenario6 <- ( Biomasses_Regional_Anisotropy_Scenario6$Biomass %o% c( 1, 1 ) + 
	Biomasses_Regional_Anisotropy_Scenario6$SE.Biomass %o% c( -1, 1 ) ) / Mean_Biomass 
ylim_max <- max( c( SE.Biomass_Regional_Anisotropy_Scenario1[,2], SE.Biomass_Regional_Anisotropy_Scenario2[,2], 
	SE.Biomass_Regional_Anisotropy_Scenario3[,2], SE.Biomass_Regional_Anisotropy_Scenario4[,2], 
	SE.Biomass_Regional_Anisotropy_Scenario5[,2], SE.Biomass_Regional_Anisotropy_Scenario6[,2] ) ) * 1.01
Colour_palette <- viridisLite::viridis( 7 )

new.graph( 0.5 )
par( cex = 1.5 )
par( mar = c( 4.5, 4.5, 1, 1 ) )
plot( Years, Biomass_Regional_Anisotropy_Scenario1, type = "n", xlab = "Year", ylab = "Relative biomass", cex.axis = 1, 
	xlim = c( min( Years ) - 1, max( Years ) + 1 ), 
	ylim = c( 0, ylim_max ) )
polygon( x = c( Years, rev( Years ) ), 
	y = c( SE.Biomass_Regional_Anisotropy_Scenario2[,1], rev( SE.Biomass_Regional_Anisotropy_Scenario2[,2] ) ), 
	col = makeTransparent( Colour_palette[1], alpha = 30 ), border = NA, lty = "solid" )
lines( Years, Biomass_Regional_Anisotropy_Scenario2, col = Colour_palette[1], type = "l", lwd = 4.5, 
	ylab = "", xlab = "", cex.axis = 0.7 )
polygon( x = c( Years, rev( Years ) ), 
	y = c( SE.Biomass_Regional_Anisotropy_Scenario1[,1], rev( SE.Biomass_Regional_Anisotropy_Scenario1[,2] ) ), 
	col = makeTransparent( Colour_palette[4], alpha = 30 ), border = NA, lty = "solid" )
lines( Years, Biomass_Regional_Anisotropy_Scenario1, col = Colour_palette[4], type = "l", lwd = 4.5, 
	ylab = "", xlab = "", cex.axis = 0.7 )
polygon( x = c( Years, rev( Years ) ), 
	y = c( SE.Biomass_Regional_Anisotropy_Scenario3[,1], rev( SE.Biomass_Regional_Anisotropy_Scenario3[,2] ) ), 
	col = makeTransparent( Colour_palette[7], alpha = 30 ), border = NA, lty = "solid" )
lines( Years, Biomass_Regional_Anisotropy_Scenario3, col = Colour_palette[7], type = "l", lwd = 4.5, 
	ylab = "", xlab = "", cex.axis = 0.7 )
abline( h = 1, col = "black", lwd = 1.3, lty = 2 )
minor.tick( nx = 5, ny = 2 )
legend( "topleft", c( "cutoff = 0.4", 
	"cutoff = 0.2", "cutoff = 0.1" ), 
	text.col = c( Colour_palette[1], Colour_palette[4], Colour_palette[7] ), 
	bty = "n", cex = 1, horiz = FALSE )
SavePlot( filename = "Indices_Regional_Anisotropy_models_Sensitivity_analysis_1", path = DIR$Figures, type = "png" )
pause( 3 ); graphics.off()

Mean( index_Regional_Anisotropy_Scenario2[2,] / index_Regional_Anisotropy_Scenario2[1,] ) #### 0.2981686
Mean( index_Regional_Anisotropy_Scenario1[2,] / index_Regional_Anisotropy_Scenario1[1,] ) #### 0.3744417
Mean( index_Regional_Anisotropy_Scenario3[2,] / index_Regional_Anisotropy_Scenario3[1,] ) #### 0.2713811

sqrt( sum( ( Regional_Anisotropy_Scenario2$rep$mu_i - Data$Catch_rate ) ^ 2 ) / length( Data$Catch_rate ) ) #### 244.4778
sqrt( sum( ( Regional_Anisotropy_Scenario1$rep$mu_i - Data$Catch_rate ) ^ 2 ) / length( Data$Catch_rate ) ) #### 217.1276
sqrt( sum( ( Regional_Anisotropy_Scenario3$rep$mu_i - Data$Catch_rate ) ^ 2 ) / length( Data$Catch_rate ) ) #### 209.0368

new.graph( 0.5 )
par( cex = 1.5 )
par( mar = c( 4.5, 4.5, 1, 1 ) )
plot( Years, Biomass_Regional_Anisotropy_Scenario1, type = "n", xlab = "Year", ylab = "Relative biomass", cex.axis = 1, 
	xlim = c( min( Years ) - 1, max( Years ) + 1 ), 
	ylim = c( 0, ylim_max ) )
polygon( x = c( Years, rev( Years ) ), 
	y = c( SE.Biomass_Regional_Anisotropy_Scenario5[,1], rev( SE.Biomass_Regional_Anisotropy_Scenario5[,2] ) ), 
	col = makeTransparent( Colour_palette[1], alpha = 30 ), border = NA, lty = "solid" )
lines( Years, Biomass_Regional_Anisotropy_Scenario5, col = Colour_palette[1], type = "l", lwd = 4.5, 
	ylab = "", xlab = "", cex.axis = 0.7 )
polygon( x = c( Years, rev( Years ) ), 
	y = c( SE.Biomass_Regional_Anisotropy_Scenario4[,1], rev( SE.Biomass_Regional_Anisotropy_Scenario4[,2] ) ), 
	col = makeTransparent( Colour_palette[4], alpha = 30 ), border = NA, lty = "solid" )
lines( Years, Biomass_Regional_Anisotropy_Scenario4, col = Colour_palette[4], type = "l", lwd = 4.5, 
	ylab = "", xlab = "", cex.axis = 0.7 )
polygon( x = c( Years, rev( Years ) ), 
	y = c( SE.Biomass_Regional_Anisotropy_Scenario6[,1], rev( SE.Biomass_Regional_Anisotropy_Scenario6[,2] ) ), 
	col = makeTransparent( Colour_palette[7], alpha = 30 ), border = NA, lty = "solid" )
lines( Years, Biomass_Regional_Anisotropy_Scenario6, col = Colour_palette[7], type = "l", lwd = 4.5, 
	ylab = "", xlab = "", cex.axis = 0.7 )
abline( h = 1, col = "black", lwd = 1.3, lty = 2 )
minor.tick( nx = 5, ny = 2 )
legend( "topleft", c( "cutoff = 0.4", 
	"cutoff = 0.2", "cutoff = 0.1" ), 
	text.col = c( Colour_palette[1], Colour_palette[4], Colour_palette[7] ), 
	bty = "n", cex = 1, horiz = FALSE )
SavePlot( filename = "Indices_Regional_Anisotropy_models_Sensitivity_analysis_2", path = DIR$Figures, type = "png" )
pause( 3 ); graphics.off()
 
Mean( index_Regional_Anisotropy_Scenario5[2,] / index_Regional_Anisotropy_Scenario5[1,] ) #### 0.2981686
Mean( index_Regional_Anisotropy_Scenario4[2,] / index_Regional_Anisotropy_Scenario4[1,] ) #### 0.3744417
Mean( index_Regional_Anisotropy_Scenario6[2,] / index_Regional_Anisotropy_Scenario6[1,] ) #### 0.2713811

sqrt( sum( ( Regional_Anisotropy_Scenario5$rep$mu_i - Data$Catch_rate ) ^ 2 ) / length( Data$Catch_rate ) ) #### 244.4778
sqrt( sum( ( Regional_Anisotropy_Scenario4$rep$mu_i - Data$Catch_rate ) ^ 2 ) / length( Data$Catch_rate ) ) #### 217.1276
sqrt( sum( ( Regional_Anisotropy_Scenario6$rep$mu_i - Data$Catch_rate ) ^ 2 ) / length( Data$Catch_rate ) ) #### 209.0368

######## Generate a mesh for a tinyVAST model using the SPDE-Barrier method and covariate-based regional anisotropy 
######## under Scenario 1 (cutoff_1 + barrier_stiffness_1) 
cutoff <- cutoff_1
barrier_stiffness <- barrier_stiffness_1
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
######## under Scenario 1 (cutoff_1 + barrier_stiffness_1) 
Regional_Anisotropy_Barrier_Scenario1 <- tinyVAST( 
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

######## Examine the fitted tinyVAST model using the SPDE-Barrier method and covariate-based regional anisotropy 
######## under Scenario 1 (cutoff_1 + barrier_stiffness_1) 
cAIC( Regional_Anisotropy_Barrier_Scenario1 ) #### 3267.455

######## Generate a mesh for a tinyVAST model using the SPDE-Barrier method and covariate-based regional anisotropy 
######## under Scenario 2 (cutoff_2 + barrier_stiffness_1) 
cutoff <- cutoff_2
barrier_stiffness <- barrier_stiffness_1
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
######## under Scenario 2 (cutoff_2 + barrier_stiffness_1) 
Regional_Anisotropy_Barrier_Scenario2 <- tinyVAST( 
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

######## Examine the fitted tinyVAST model using the SPDE-Barrier method and covariate-based regional anisotropy 
######## under Scenario 2 (cutoff_2 + barrier_stiffness_1) 
cAIC( Regional_Anisotropy_Barrier_Scenario2 ) #### 3309.503

######## Generate a mesh for a tinyVAST model using the SPDE-Barrier method and covariate-based regional anisotropy 
######## under Scenario 3 (cutoff_3 + barrier_stiffness_1) 
cutoff <- cutoff_3
barrier_stiffness <- barrier_stiffness_1
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
######## under Scenario 3 (cutoff_3 + barrier_stiffness_1) 
Regional_Anisotropy_Barrier_Scenario3 <- tinyVAST( 
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

######## Examine the fitted tinyVAST model using the SPDE-Barrier method and covariate-based regional anisotropy 
######## under Scenario 3 (cutoff_3 + barrier_stiffness_1) 
cAIC( Regional_Anisotropy_Barrier_Scenario3 ) #### 3281.003

######## Generate a mesh for a tinyVAST model using the SPDE-Barrier method and covariate-based regional anisotropy 
######## under Scenario 4 (cutoff_1 + barrier_stiffness_2) 
cutoff <- cutoff_1
barrier_stiffness <- barrier_stiffness_2
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
######## under Scenario 4 (cutoff_1 + barrier_stiffness_2) 
Regional_Anisotropy_Barrier_Scenario4 <- tinyVAST( 
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

######## Examine the fitted tinyVAST model using the SPDE-Barrier method and covariate-based regional anisotropy 
######## under Scenario 4 (cutoff_1 + barrier_stiffness_2) 
cAIC( Regional_Anisotropy_Barrier_Scenario4 ) #### 3268.307

######## Generate a mesh for a tinyVAST model using the SPDE-Barrier method and covariate-based regional anisotropy 
######## under Scenario 5 (cutoff_2 + barrier_stiffness_2) 
cutoff <- cutoff_2
barrier_stiffness <- barrier_stiffness_2
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
######## under Scenario 5 (cutoff_2 + barrier_stiffness_2) 
Regional_Anisotropy_Barrier_Scenario5 <- tinyVAST( 
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

######## Examine the fitted tinyVAST model using the SPDE-Barrier method and covariate-based regional anisotropy 
######## under Scenario 5 (cutoff_2 + barrier_stiffness_2) 
cAIC( Regional_Anisotropy_Barrier_Scenario5 ) #### 3309.411

######## Generate a mesh for a tinyVAST model using the SPDE-Barrier method and covariate-based regional anisotropy 
######## under Scenario 6 (cutoff_3 + barrier_stiffness_2) 
cutoff <- cutoff_3
barrier_stiffness <- barrier_stiffness_2
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
######## under Scenario 6 (cutoff_3 + barrier_stiffness_2) 
Regional_Anisotropy_Barrier_Scenario6 <- tinyVAST( 
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

######## Examine the fitted tinyVAST model using the SPDE-Barrier method and covariate-based regional anisotropy 
######## under Scenario 6 (cutoff_3 + barrier_stiffness_2) 
cAIC( Regional_Anisotropy_Barrier_Scenario6 ) #### 3281.949

######## Save the fits of the tinyVAST models using the SPDE-Barrier method and covariate-based regional anisotropy 
######## under Scenarios 1-6 
save( Regional_Anisotropy_Barrier_Scenario1, Regional_Anisotropy_Barrier_Scenario2, Regional_Anisotropy_Barrier_Scenario3,
	Regional_Anisotropy_Barrier_Scenario4, Regional_Anisotropy_Barrier_Scenario5, Regional_Anisotropy_Barrier_Scenario6, 
	file = make.filename( "Regional_Anisotropy_Barrier_models_West_Coast_South_Island.RData", DIR$Output ) )

######## Generate predictions for each fishing year of the time period with the estimated tinyVAST models 
######## using the SPDE-Barrier method and covariate-based regional anisotropy under Scenarios 1-6, 
######## and save these predictions in a .RData file
#### Load the estimated tinyVAST models and define prediction grids
load( file = make.filename( "Regional_Anisotropy_Barrier_models_West_Coast_South_Island.RData", DIR$Output ) )
load( file = make.filename( "Prediction_grid_SNA7_5kmx5km.RData", DIR$Input ) )

#### Generate predictions for each fishing year of the time period with the estimated models
Years <- as.numeric( min( Data$Year ) : max( Data$Year ) )
areas_km2 <- 25
index_Regional_Anisotropy_Barrier_Scenario1 <- index_Regional_Anisotropy_Barrier_Scenario2 <- 
	index_Regional_Anisotropy_Barrier_Scenario3 <- index_Regional_Anisotropy_Barrier_Scenario4 <- 
	index_Regional_Anisotropy_Barrier_Scenario5 <- index_Regional_Anisotropy_Barrier_Scenario6 <- NULL

pb <- txtProgressBar( min = 0, max = length( Years ), style = 3 )
yr <- 0

for( Year in Years ) {

	yr <- yr + 1
	setTxtProgressBar( pb, yr )

	## Produce predictive dataframes
	newdata <- data.frame( "Year" = Years[yr], 
		"Lon" = Prediction_grid[,'Lon'], 
		"Lat" = Prediction_grid[,'Lat'],
		"var" = "density" )

  	## Generate predictions 
	total_Regional_Anisotropy_Barrier_Scenario1 <- integrate_output( Regional_Anisotropy_Barrier_Scenario1, 
		newdata = newdata, area = areas_km2 )
	index_Regional_Anisotropy_Barrier_Scenario1 <- cbind( index_Regional_Anisotropy_Barrier_Scenario1, 
		total_Regional_Anisotropy_Barrier_Scenario1[c( 'Est. (bias.correct)', 'Std. Error' )] )
	total_Regional_Anisotropy_Barrier_Scenario2 <- integrate_output( Regional_Anisotropy_Barrier_Scenario2, 
		newdata = newdata, area = areas_km2 )
	index_Regional_Anisotropy_Barrier_Scenario2 <- cbind( index_Regional_Anisotropy_Barrier_Scenario2, 
		total_Regional_Anisotropy_Barrier_Scenario2[c( 'Est. (bias.correct)', 'Std. Error' )] )
	total_Regional_Anisotropy_Barrier_Scenario3 <- integrate_output( Regional_Anisotropy_Barrier_Scenario3, 
		newdata = newdata, area = areas_km2 )
	index_Regional_Anisotropy_Barrier_Scenario3 <- cbind( index_Regional_Anisotropy_Barrier_Scenario3, 
		total_Regional_Anisotropy_Barrier_Scenario3[c( 'Est. (bias.correct)', 'Std. Error' )] )
	total_Regional_Anisotropy_Barrier_Scenario4 <- integrate_output( Regional_Anisotropy_Barrier_Scenario4, 
		newdata = newdata, area = areas_km2 )
	index_Regional_Anisotropy_Barrier_Scenario4 <- cbind( index_Regional_Anisotropy_Barrier_Scenario4, 
		total_Regional_Anisotropy_Barrier_Scenario4[c( 'Est. (bias.correct)', 'Std. Error' )] )
	total_Regional_Anisotropy_Barrier_Scenario5 <- integrate_output( Regional_Anisotropy_Barrier_Scenario5, 
		newdata = newdata, area = areas_km2 )
	index_Regional_Anisotropy_Barrier_Scenario5 <- cbind( index_Regional_Anisotropy_Barrier_Scenario5, 
		total_Regional_Anisotropy_Barrier_Scenario5[c( 'Est. (bias.correct)', 'Std. Error' )] )
	total_Regional_Anisotropy_Barrier_Scenario6 <- integrate_output( Regional_Anisotropy_Barrier_Scenario6, 
		newdata = newdata, area = areas_km2 )
	index_Regional_Anisotropy_Barrier_Scenario6 <- cbind( index_Regional_Anisotropy_Barrier_Scenario6, 
		total_Regional_Anisotropy_Barrier_Scenario6[c( 'Est. (bias.correct)', 'Std. Error' )] )
	
}

close( pb )	

#### Save predictions in a .RData file 
save( Years, index_Regional_Anisotropy_Barrier_Scenario1, index_Regional_Anisotropy_Barrier_Scenario2, 
	index_Regional_Anisotropy_Barrier_Scenario3, index_Regional_Anisotropy_Barrier_Scenario4, 
	index_Regional_Anisotropy_Barrier_Scenario5, index_Regional_Anisotropy_Barrier_Scenario6, 
	file = make.filename( "Predictions_Regional_Anisotropy_Barrier_models_West_Coast_South_Island.RData", DIR$Output ) )

######## Produce figures showing the indices of relative biomass predicted by tinyVAST models 
######## using the SPDE-Barrier method and covariate-based regional anisotropy under Scenarios 1-6, and save this figure
#### Load estimated models and their predictions
load( file = make.filename( "Regional_Anisotropy_Barrier_models_West_Coast_South_Island.RData", DIR$Output ) )
load( file = make.filename( "Predictions_Regional_Anisotropy_Barrier_models_West_Coast_South_Island.RData", DIR$Output ) )

#### Produce the figures and save them 
Biomasses_Regional_Anisotropy_Barrier_Scenario1 <- data.frame( "year" = Years, Biomass = index_Regional_Anisotropy_Barrier_Scenario1[1,], 
	SE.Biomass = index_Regional_Anisotropy_Barrier_Scenario1[2,] )
Biomasses_Regional_Anisotropy_Barrier_Scenario2 <- data.frame( "year" = Years, Biomass = index_Regional_Anisotropy_Barrier_Scenario2[1,], 
	SE.Biomass = index_Regional_Anisotropy_Barrier_Scenario2[2,] )
Biomasses_Regional_Anisotropy_Barrier_Scenario3 <- data.frame( "year" = Years, Biomass = index_Regional_Anisotropy_Barrier_Scenario3[1,], 
	SE.Biomass = index_Regional_Anisotropy_Barrier_Scenario3[2,] )
Biomasses_Regional_Anisotropy_Barrier_Scenario4 <- data.frame( "year" = Years, Biomass = index_Regional_Anisotropy_Barrier_Scenario4[1,], 
	SE.Biomass = index_Regional_Anisotropy_Barrier_Scenario4[2,] )
Biomasses_Regional_Anisotropy_Barrier_Scenario5 <- data.frame( "year" = Years, Biomass = index_Regional_Anisotropy_Barrier_Scenario5[1,], 
	SE.Biomass = index_Regional_Anisotropy_Barrier_Scenario5[2,] )
Biomasses_Regional_Anisotropy_Barrier_Scenario6 <- data.frame( "year" = Years, Biomass = index_Regional_Anisotropy_Barrier_Scenario6[1,], 
	SE.Biomass = index_Regional_Anisotropy_Barrier_Scenario6[2,] )
Mean_Biomass <- Mean( c( Biomasses_Regional_Anisotropy_Barrier_Scenario1$Biomass, Biomasses_Regional_Anisotropy_Barrier_Scenario2$Biomass, 
	Biomasses_Regional_Anisotropy_Barrier_Scenario3$Biomass, Biomasses_Regional_Anisotropy_Barrier_Scenario4$Biomass, 
	Biomasses_Regional_Anisotropy_Barrier_Scenario5$Biomass, Biomasses_Regional_Anisotropy_Barrier_Scenario6$Biomass ) ) 
Biomass_Regional_Anisotropy_Barrier_Scenario1 <- Biomasses_Regional_Anisotropy_Barrier_Scenario1$Biomass / Mean_Biomass 
Biomass_Regional_Anisotropy_Barrier_Scenario2 <- Biomasses_Regional_Anisotropy_Barrier_Scenario2$Biomass / Mean_Biomass 
Biomass_Regional_Anisotropy_Barrier_Scenario3 <- Biomasses_Regional_Anisotropy_Barrier_Scenario3$Biomass / Mean_Biomass 
Biomass_Regional_Anisotropy_Barrier_Scenario4 <- Biomasses_Regional_Anisotropy_Barrier_Scenario4$Biomass / Mean_Biomass 
Biomass_Regional_Anisotropy_Barrier_Scenario5 <- Biomasses_Regional_Anisotropy_Barrier_Scenario5$Biomass / Mean_Biomass 
Biomass_Regional_Anisotropy_Barrier_Scenario6 <- Biomasses_Regional_Anisotropy_Barrier_Scenario6$Biomass / Mean_Biomass 
SE.Biomass_Regional_Anisotropy_Barrier_Scenario1 <- ( Biomasses_Regional_Anisotropy_Barrier_Scenario1$Biomass %o% c( 1, 1 ) + 
	Biomasses_Regional_Anisotropy_Barrier_Scenario1$SE.Biomass %o% c( -1, 1 ) ) / Mean_Biomass 
SE.Biomass_Regional_Anisotropy_Barrier_Scenario2 <- ( Biomasses_Regional_Anisotropy_Barrier_Scenario2$Biomass %o% c( 1, 1 ) + 
	Biomasses_Regional_Anisotropy_Barrier_Scenario2$SE.Biomass %o% c( -1, 1 ) ) / Mean_Biomass 
SE.Biomass_Regional_Anisotropy_Barrier_Scenario3 <- ( Biomasses_Regional_Anisotropy_Barrier_Scenario3$Biomass %o% c( 1, 1 ) + 
	Biomasses_Regional_Anisotropy_Barrier_Scenario3$SE.Biomass %o% c( -1, 1 ) ) / Mean_Biomass 
SE.Biomass_Regional_Anisotropy_Barrier_Scenario4 <- ( Biomasses_Regional_Anisotropy_Barrier_Scenario4$Biomass %o% c( 1, 1 ) + 
	Biomasses_Regional_Anisotropy_Barrier_Scenario4$SE.Biomass %o% c( -1, 1 ) ) / Mean_Biomass 
SE.Biomass_Regional_Anisotropy_Barrier_Scenario5 <- ( Biomasses_Regional_Anisotropy_Barrier_Scenario5$Biomass %o% c( 1, 1 ) + 
	Biomasses_Regional_Anisotropy_Barrier_Scenario5$SE.Biomass %o% c( -1, 1 ) ) / Mean_Biomass 
SE.Biomass_Regional_Anisotropy_Barrier_Scenario6 <- ( Biomasses_Regional_Anisotropy_Barrier_Scenario6$Biomass %o% c( 1, 1 ) + 
	Biomasses_Regional_Anisotropy_Barrier_Scenario6$SE.Biomass %o% c( -1, 1 ) ) / Mean_Biomass 
ylim_max <- max( c( SE.Biomass_Regional_Anisotropy_Barrier_Scenario1[,2], SE.Biomass_Regional_Anisotropy_Barrier_Scenario2[,2], 
	SE.Biomass_Regional_Anisotropy_Barrier_Scenario3[,2], SE.Biomass_Regional_Anisotropy_Barrier_Scenario4[,2], 
	SE.Biomass_Regional_Anisotropy_Barrier_Scenario5[,2], SE.Biomass_Regional_Anisotropy_Barrier_Scenario6[,2] ) ) * 1.01
Colour_palette <- viridisLite::viridis( 7 )

new.graph( 0.5 )
par( cex = 1.5 )
par( mar = c( 4.5, 4.5, 1, 1 ) )
plot( Years, Biomass_Regional_Anisotropy_Barrier_Scenario1, type = "n", xlab = "Year", ylab = "Relative biomass", cex.axis = 1, 
	xlim = c( min( Years ) - 1, max( Years ) + 1 ), 
	ylim = c( 0, ylim_max ) )
polygon( x = c( Years, rev( Years ) ), 
	y = c( SE.Biomass_Regional_Anisotropy_Barrier_Scenario2[,1], rev( SE.Biomass_Regional_Anisotropy_Barrier_Scenario2[,2] ) ), 
	col = makeTransparent( Colour_palette[1], alpha = 30 ), border = NA, lty = "solid" )
lines( Years, Biomass_Regional_Anisotropy_Barrier_Scenario2, col = Colour_palette[1], type = "l", lwd = 4.5, 
	ylab = "", xlab = "", cex.axis = 0.7 )
polygon( x = c( Years, rev( Years ) ), 
	y = c( SE.Biomass_Regional_Anisotropy_Barrier_Scenario1[,1], rev( SE.Biomass_Regional_Anisotropy_Barrier_Scenario1[,2] ) ), 
	col = makeTransparent( Colour_palette[4], alpha = 30 ), border = NA, lty = "solid" )
lines( Years, Biomass_Regional_Anisotropy_Barrier_Scenario1, col = Colour_palette[4], type = "l", lwd = 4.5, 
	ylab = "", xlab = "", cex.axis = 0.7 )
polygon( x = c( Years, rev( Years ) ), 
	y = c( SE.Biomass_Regional_Anisotropy_Barrier_Scenario3[,1], rev( SE.Biomass_Regional_Anisotropy_Barrier_Scenario3[,2] ) ), 
	col = makeTransparent( Colour_palette[7], alpha = 30 ), border = NA, lty = "solid" )
lines( Years, Biomass_Regional_Anisotropy_Barrier_Scenario3, col = Colour_palette[7], type = "l", lwd = 4.5, 
	ylab = "", xlab = "", cex.axis = 0.7 )
abline( h = 1, col = "black", lwd = 1.3, lty = 2 )
minor.tick( nx = 5, ny = 2 )
legend( "topleft", c( "cutoff = 0.4", 
	"cutoff = 0.2", "cutoff = 0.1" ), 
	text.col = c( Colour_palette[1], Colour_palette[4], Colour_palette[7] ), 
	bty = "n", cex = 1, horiz = FALSE )
SavePlot( filename = "Indices_Regional_Anisotropy_Barrier_models_Sensitivity_analysis_1", path = DIR$Figures, type = "png" )
pause( 3 ); graphics.off()

Mean( index_Regional_Anisotropy_Barrier_Scenario2[2,] / index_Regional_Anisotropy_Barrier_Scenario2[1,] ) #### 0.3038319
Mean( index_Regional_Anisotropy_Barrier_Scenario1[2,] / index_Regional_Anisotropy_Barrier_Scenario1[1,] ) #### 0.2915658
Mean( index_Regional_Anisotropy_Barrier_Scenario3[2,] / index_Regional_Anisotropy_Barrier_Scenario3[1,] ) #### 0.2694122

sqrt( sum( ( Regional_Anisotropy_Barrier_Scenario2$rep$mu_i - Data$Catch_rate ) ^ 2 ) / length( Data$Catch_rate ) ) #### 217.2727
sqrt( sum( ( Regional_Anisotropy_Barrier_Scenario1$rep$mu_i - Data$Catch_rate ) ^ 2 ) / length( Data$Catch_rate ) ) #### 202.673
sqrt( sum( ( Regional_Anisotropy_Barrier_Scenario3$rep$mu_i - Data$Catch_rate ) ^ 2 ) / length( Data$Catch_rate ) ) #### 187.6917

new.graph( 0.5 )
par( cex = 1.5 )
par( mar = c( 4.5, 4.5, 1, 1 ) )
plot( Years, Biomass_Regional_Anisotropy_Barrier_Scenario1, type = "n", xlab = "Year", ylab = "Relative biomass", cex.axis = 1, 
	xlim = c( min( Years ) - 1, max( Years ) + 1 ), 
	ylim = c( 0, ylim_max ) )
polygon( x = c( Years, rev( Years ) ), 
	y = c( SE.Biomass_Regional_Anisotropy_Barrier_Scenario5[,1], rev( SE.Biomass_Regional_Anisotropy_Barrier_Scenario5[,2] ) ), 
	col = makeTransparent( Colour_palette[1], alpha = 30 ), border = NA, lty = "solid" )
lines( Years, Biomass_Regional_Anisotropy_Barrier_Scenario5, col = Colour_palette[1], type = "l", lwd = 4.5, 
	ylab = "", xlab = "", cex.axis = 0.7 )
polygon( x = c( Years, rev( Years ) ), 
	y = c( SE.Biomass_Regional_Anisotropy_Barrier_Scenario4[,1], rev( SE.Biomass_Regional_Anisotropy_Barrier_Scenario4[,2] ) ), 
	col = makeTransparent( Colour_palette[4], alpha = 30 ), border = NA, lty = "solid" )
lines( Years, Biomass_Regional_Anisotropy_Barrier_Scenario4, col = Colour_palette[4], type = "l", lwd = 4.5, 
	ylab = "", xlab = "", cex.axis = 0.7 )
polygon( x = c( Years, rev( Years ) ), 
	y = c( SE.Biomass_Regional_Anisotropy_Barrier_Scenario6[,1], rev( SE.Biomass_Regional_Anisotropy_Barrier_Scenario6[,2] ) ), 
	col = makeTransparent( Colour_palette[7], alpha = 30 ), border = NA, lty = "solid" )
lines( Years, Biomass_Regional_Anisotropy_Barrier_Scenario6, col = Colour_palette[7], type = "l", lwd = 4.5, 
	ylab = "", xlab = "", cex.axis = 0.7 )
abline( h = 1, col = "black", lwd = 1.3, lty = 2 )
minor.tick( nx = 5, ny = 2 )
legend( "topleft", c( "cutoff = 0.4", 
	"cutoff = 0.2", "cutoff = 0.1" ), 
	text.col = c( Colour_palette[1], Colour_palette[4], Colour_palette[7] ), 
	bty = "n", cex = 1, horiz = FALSE )
SavePlot( filename = "Indices_Regional_Anisotropy_Barrier_models_Sensitivity_analysis_2", path = DIR$Figures, type = "png" )
pause( 3 ); graphics.off()
 
Mean( index_Regional_Anisotropy_Barrier_Scenario5[2,] / index_Regional_Anisotropy_Barrier_Scenario5[1,] ) #### 0.3028409
Mean( index_Regional_Anisotropy_Barrier_Scenario4[2,] / index_Regional_Anisotropy_Barrier_Scenario4[1,] ) #### 0.2912545
Mean( index_Regional_Anisotropy_Barrier_Scenario6[2,] / index_Regional_Anisotropy_Barrier_Scenario6[1,] ) #### 0.2687169

sqrt( sum( ( Regional_Anisotropy_Barrier_Scenario5$rep$mu_i - Data$Catch_rate ) ^ 2 ) / length( Data$Catch_rate ) ) #### 217.197
sqrt( sum( ( Regional_Anisotropy_Barrier_Scenario4$rep$mu_i - Data$Catch_rate ) ^ 2 ) / length( Data$Catch_rate ) ) #### 202.4292
sqrt( sum( ( Regional_Anisotropy_Barrier_Scenario6$rep$mu_i - Data$Catch_rate ) ^ 2 ) / length( Data$Catch_rate ) ) #### 187.7346


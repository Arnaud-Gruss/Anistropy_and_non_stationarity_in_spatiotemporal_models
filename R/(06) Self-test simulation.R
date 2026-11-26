###############################################################################################################################
##
##  (06) Self-test simulation 
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
library( tictoc )
library( tidyr )

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

######## Define settings for the tinyVAST models
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

######## Generate a mesh for the tinyVAST operating model (OM) using the SPDE method and assuming isotropy
mesh <- fm_mesh_2d( Mesh_data[,c( 'Lon', 'Lat' )], cutoff = cutoff )

plot( mesh )
plot( map_sf, add = TRUE )
pause( 3 ); graphics.off()

######## Fit the tinyVAST OM using the SPDE method and assuming isotropy and save the model fit in a .RData file
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
save( SPDE_Isotropy, file = make.filename( "SPDE_Isotropy.RData", DIR$Output ) )

######## Simulate new data with the fitted OM using the SPDE method and assuming isotropy and save the
######## simulated data in a .RData file
Sim_Data <- Data
fit <- SPDE_Isotropy
Sim_SPDE_Isotropy <- simulate( fit, nsim = 100, type = "mle-mvn" )
save( Sim_Data, mesh, Sim_SPDE_Isotropy, file = make.filename( "Sim_SPDE_Isotropy.RData", DIR$Output ) )

######## Fit tinyVAST estimation models (EMs) using the SPDE method and assuming isotropy and save the model fits 
######## in .RData files
load( file = make.filename( "Sim_SPDE_Isotropy.RData", DIR$Input ) )

tic( "Estimation of 100 models" )

pgBar <- txtProgressBar( min = 1, max = 100, style = 3 )

for ( i in 1 : 100 ) {

	setTxtProgressBar( pgBar, i )

	Sim_Data$Catch_rate <- Sim_SPDE_Isotropy[,i]

	refit <- tinyVAST( 
		data = droplevels( Sim_Data ), 
		formula = formula, 
		space_term = space_term,
  		spacetime_term = spacetime_term,
  		time_term = time_term,
  		spatial_domain = mesh,
  		family = tweedie( "log" ),
		space_columns = c( "Lon", "Lat" ),
  		time_column = "Year",
  		times = min( Sim_Data$Year ) : max( Sim_Data$Year ),
  		variable_column = "var",
		control = tinyVASTcontrol( 
			verbose = FALSE,  
			newton_loops = 2, 
    			use_anisotropy = FALSE
		)
	)

	save( refit, file = make.filename( paste0( "Refit_SPDE_Isotropy_Sim", i, ".RData" ), DIR$Output ) )

}

close( pgBar )
Sys.sleep( 1.5 )
toc()

######## Generate a mesh for the tinyVAST OM using the SPDE method and covariate-based regional anisotropy 
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

######## Fit the tinyVAST OM using the SPDE method and covariate-based regional anisotropy 
######## and save the model fit in a .RData file
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
save( Regional_Anisotropy, file = make.filename( "Regional_Anisotropy.RData", DIR$Output ) )

######## Simulate new data with the fitted OM using the SPDE method and covariate-based regional anisotropy and save the
######## simulated data in a .RData file
Sim_Data <- Data
fit <- Regional_Anisotropy
Sim_Regional_Anisotropy <- simulate( fit, nsim = 100, type = "mle-mvn" )
save( Sim_Data, mesh_cov, Sim_Regional_Anisotropy, file = make.filename( "Sim_Regional_Anisotropy.RData", DIR$Output ) )

######## Fit tinyVAST EMs using the SPDE method and covariate-based regional anisotropy and save the model fits 
######## in .RData files
load( file = make.filename( "Sim_Regional_Anisotropy.RData", DIR$Input ) )

tic( "Estimation of 100 models" )

pgBar <- txtProgressBar( min = 1, max = 100, style = 3 )

for ( i in 1 : 100 ) {

	setTxtProgressBar( pgBar, i )

	Sim_Data$Catch_rate <- Sim_Regional_Anisotropy[,i]

	refit <- tinyVAST( 
		data = droplevels( Sim_Data ), 
		formula = formula, 
		space_term = space_term,
  		spacetime_term = spacetime_term,
  		time_term = time_term,
  		spatial_domain = mesh_cov,
  		family = tweedie( "log" ),
		space_columns = c( "Lon", "Lat" ),
  		time_column = "Year",
  		times = min( Sim_Data$Year ) : max( Sim_Data$Year ),
  		variable_column = "var",
		development = list(
			vertex_formula = ~
				I(x_tbgb + y_tbgb) + y_tbgb +
					I(x_wc + y_wc) + y_wc +
					0 
  		),
		control = tinyVASTcontrol( 
			verbose = FALSE, 
			newton_loops = 2, 
    			use_anisotropy = FALSE
		)
	)

	save( refit, file = make.filename( paste0( "Refit_Regional_Anisotropy_Sim", i, ".RData" ), DIR$Output ) )

}

close( pgBar )
Sys.sleep( 1.5 )
toc()

######## Generate a mesh for the tinyVAST OM using the SPDE-Barrier method and covariate-based regional anisotropy 
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

######## Fit the tinyVAST OM using the SPDE-Barrier method and covariate-based regional anisotropy 
######## and save the model fit in a .RData file
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
save( Regional_Anisotropy_Barrier, file = make.filename( "Regional_Anisotropy_Barrier.RData", DIR$Output ) )

######## Simulate new data with the fitted OM using the SPDE-Barrier method and covariate-based regional anisotropy 
######## and save the simulated data in a .RData file
Sim_Data <- Data
fit <- Regional_Anisotropy_Barrier
Sim_Regional_Anisotropy_Barrier <- simulate( fit, nsim = 100, type = "mle-mvn" )
save( Sim_Data, mesh_cov, Sim_Regional_Anisotropy_Barrier, 
	file = make.filename( "Sim_Regional_Anisotropy_Barrier.RData", DIR$Output ) )

######## Fit tinyVAST EMs using the SPDE-Barrier method and covariate-based regional anisotropy and save the model fits 
######## in .RData files
load( file = make.filename( "Sim_Regional_Anisotropy_Barrier.RData", DIR$Input ) )

tic( "Estimation of 100 models" )

pgBar <- txtProgressBar( min = 1, max = 100, style = 3 )

for ( i in 1 : 100 ) {

	setTxtProgressBar( pgBar, i )

	Sim_Data$Catch_rate <- Sim_Regional_Anisotropy_Barrier[,i]

	refit <- tinyVAST( 
		data = droplevels( Sim_Data ), 
		formula = formula, 
		space_term = space_term,
  		spacetime_term = spacetime_term,
  		time_term = time_term,
  		spatial_domain = mesh_cov,
  		family = tweedie( "log" ),
		space_columns = c( "Lon", "Lat" ),
  		time_column = "Year",
  		times = min( Sim_Data$Year ) : max( Sim_Data$Year ),
  		variable_column = "var",
		development = list(
			triangle_formula = ~ offset(barrier_polygon),
			vertex_formula = ~
				I(x_tbgb + y_tbgb) + y_tbgb +
					I(x_wc + y_wc) + y_wc +
					0 
  		),
		control = tinyVASTcontrol( 
			verbose = FALSE, 
			newton_loops = 2, 
    			use_anisotropy = FALSE, 
			barrier_stiffness = barrier_stiffness
		)
	)

	save( refit, file = make.filename( paste0( "Refit_Regional_Anisotropy_Barrier_Sim", i, ".RData" ), DIR$Output ) )

}

close( pgBar )
Sys.sleep( 1.5 )
toc()

######## Process the outputs of the self-test experimental for the tinyVAST model using the SPDE method and 
######## assuming isotropy, to produce relative error estimates for this model and save these 
######## relative error estimates in a .RData file
load( file = make.filename( "SPDE_Isotropy.RData", DIR$Input ) )
fit <- SPDE_Isotropy
Parameters_fit <- as.data.frame( t( fit$sdrep$par.fixed[1 : 7] ) )
names( Parameters_fit )[5 : 6] <- c( "log_sigma_1", "log_sigma_2" )
Parameters_fit <- dplyr::bind_rows( replicate( 100, Parameters_fit, simplify = FALSE ) )

Parameters_Refit <- data.frame()
pgBar <- txtProgressBar( min = 1, max = 100, style = 3 )
for ( i in 1 : 100 ) {
	setTxtProgressBar( pgBar, i )
	load( file = make.filename( paste0( "Refit_SPDE_Isotropy_Sim", i, ".RData" ), DIR$Output ) )
	if ( sanity( refit, silent = TRUE )$all_ok == TRUE ) {
		Parameters_Refit <- rbind( Parameters_Refit, t( refit$sdrep$par.fixed[1 : 7] ) )
	}
}
close( pgBar ) 
nrow( Parameters_Refit ) #### 100% of the EMs passed the sanity checks 
names( Parameters_Refit )[5 : 6] <- c( "log_sigma_1", "log_sigma_2" )

Par_fit <- as.matrix( Parameters_fit )
Par_refit <- as.matrix( Parameters_Refit )
Par_RE <- as.data.frame( ( Par_refit - Par_fit ) / Par_fit )
names( Par_RE ) <- names( Parameters_fit ) 

save( Par_RE, file = make.filename( "Par_RE_SPDE_Isotropy.RData", DIR$Output ) )

######## Process the outputs of the self-test experimental for the tinyVAST model using the SPDE method and 
######## covariate-based regional anisotropy, to produce relative error estimates for this model and save these 
######## relative error estimates in a .RData file
load( file = make.filename( "Regional_Anisotropy.RData", DIR$Input ) )
fit <- Regional_Anisotropy
Parameters_fit <- as.data.frame( t( fit$sdrep$par.fixed[1 : 7] ) )
names( Parameters_fit )[5 : 6] <- c( "log_sigma_1", "log_sigma_2" )
Parameters_fit <- dplyr::bind_rows( replicate( 100, Parameters_fit, simplify = FALSE ) )

Parameters_Refit <- data.frame()
Replicates_dropped <- c()
pgBar <- txtProgressBar( min = 1, max = 100, style = 3 )
for ( i in 1 : 100 ) {
	setTxtProgressBar( pgBar, i )
	load( file = make.filename( paste0( "Refit_Regional_Anisotropy_Sim", i, ".RData" ), DIR$Output ) )
	if ( sanity( refit, silent = TRUE )$all_ok == TRUE ) {
		Parameters_Refit <- rbind( Parameters_Refit, t( refit$sdrep$par.fixed[1 : 7] ) )
	} else {
		Replicates_dropped <- c( Replicates_dropped, i )
	}
}
close( pgBar ) 
nrow( Parameters_Refit ) #### 86% of the EMs passed the sanity checks 
names( Parameters_Refit )[5 : 6] <- c( "log_sigma_1", "log_sigma_2" )

Par_fit <- as.matrix( Parameters_fit[-Replicates_dropped,] )
Par_refit <- as.matrix( Parameters_Refit )
Par_RE <- as.data.frame( ( Par_refit - Par_fit ) / Par_fit )
names( Par_RE ) <- names( Parameters_fit )

save( Par_RE, file = make.filename( "Par_RE_Regional_Anisotropy.RData", DIR$Output ) )

######## Process the outputs of the self-test experimental for the tinyVAST model using the SPDE-Barrier method and 
######## covariate-based regional anisotropy, to produce relative error estimates for this model and save these 
######## relative error estimates in a .RData file
load( file = make.filename( "Regional_Anisotropy_Barrier.RData", DIR$Input ) )
fit <- Regional_Anisotropy_Barrier
Parameters_fit <- as.data.frame( t( fit$sdrep$par.fixed[1 : 7] ) )
names( Parameters_fit )[5 : 6] <- c( "log_sigma_1", "log_sigma_2" )
Parameters_fit <- dplyr::bind_rows( replicate( 100, Parameters_fit, simplify = FALSE ) )

Parameters_Refit <- data.frame()
Replicates_dropped <- c()
pgBar <- txtProgressBar( min = 1, max = 100, style = 3 )
for ( i in 1 : 100 ) {
	setTxtProgressBar( pgBar, i )
	load( file = make.filename( paste0( "Refit_Regional_Anisotropy_Barrier_Sim", i, ".RData" ), DIR$Output ) )
	if ( sanity( refit, silent = TRUE )$all_ok == TRUE ) {
		Parameters_Refit <- rbind( Parameters_Refit, t( refit$sdrep$par.fixed[1 : 7] ) )
	} else {
		Replicates_dropped <- c( Replicates_dropped, i )
	}
}
close( pgBar ) 
nrow( Parameters_Refit ) #### 81% of the EMs passed the sanity checks 
names( Parameters_Refit )[5 : 6] <- c( "log_sigma_1", "log_sigma_2" )

Par_fit <- as.matrix( Parameters_fit[-Replicates_dropped,] )
Par_refit <- as.matrix( Parameters_Refit )
Par_RE <- as.data.frame( ( Par_refit - Par_fit ) / Par_fit )
names( Par_RE ) <- names( Parameters_fit )

save( Par_RE, file = make.filename( "Par_RE_Regional_Anisotropy_Barrier.RData", DIR$Output ) )

######## Produce a figure summarising the results of the self-test simulation experiment in terms of relative error
######## and save this figure 
load( file = make.filename( "Par_RE_SPDE_Isotropy.RData", DIR$Output ) )
long_df <- Par_RE %>%
	pivot_longer( cols = everything(), names_to = "Parameter", values_to = "Value" ) %>%
	mutate( "Model" = "Isotropy" )
Parameters_RE <- long_df

load( file = make.filename( "Par_RE_Regional_Anisotropy.RData", DIR$Output ) )
long_df <- Par_RE %>% 
	pivot_longer( cols = everything(), names_to = "Parameter", values_to = "Value" ) %>%
	mutate( "Model" = "Anisotropy" )
Parameters_RE <- rbind( Parameters_RE, long_df )

load( file = make.filename( "Par_RE_Regional_Anisotropy_Barrier.RData", DIR$Output ) )
long_df <- Par_RE %>%
	pivot_longer( cols = everything(), names_to = "Parameter", values_to = "Value" ) %>%
	mutate( "Model" = "Anisotropy + barriers" )
Parameters_RE <- rbind( Parameters_RE, long_df )

Parameters_RE$Model <- factor( Parameters_RE$Model, 
	levels = c( "Isotropy", "Anisotropy", "Anisotropy + barriers" ) )

Parameters_Relative_Error_Figure <- ggplot( Parameters_RE, aes( x = Parameter, y = Value, fill = Model ) ) +
	geom_boxplot(  ) +
  	labs(
		y = "Relative error",
    		x = "Parameter\n",
    		fill = "Model" ) +
	geom_hline( yintercept = 0 ) +
  	theme_bw()+
  	coord_flip() +
	theme( legend.position = "bottom", 
		axis.title = element_text( size = rel( 1.5 ), face = "bold" ),
		axis.text.x = element_text( size = 12 ),
		axis.text.y = element_text( size = 12 ), 
		axis.ticks.length = unit( 7.5, "pt" ),
    		axis.minor.ticks.length = rel( 0.75 ), 
		legend.title = element_text( size = 15, face = "bold" ), 
		legend.text = element_text(size = 13 ) )

new.graph( 0.5 )
print( Parameters_Relative_Error_Figure ) 
SavePlot( filename = "Results of the self-test simulation experiment", path = DIR$Figures, type = "png" )
pause( 2 ); graphics.off()


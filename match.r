library("lpSolve") 

match.cc <- function(cases,controls,features,id.column=NA,controls_per_case=1,distance.measure="mahalanobis",match.method="gsa",caliper=NA) {
	# Get columns of a dataframe by column names
	subset_named_columns <- function(df,columns) {
		column_index_vector <- as.vector(sapply(features,function(f) { which(colnames(df) == f) }))
		df2 <- as.data.frame(df[,column_index_vector])
		names(df2) <- columns
		df2
	}


	#####################################################
	# Distance calculation using eucledian distances
	#####################################################
	
	std <- function(x){
		if(length(which(is.na(x)))==0)
	    	(x-mean(x))/sd(x)
		else
			(x-mean(x,na.rm=T))/sd(x,na.rm=T)
	}

	normalize.df <- function(df,colnames) {
 		for(colname in colnames) {
			c <- which(colnames(df)==colname)
			df[,c] <- std(df[,c])
	  	}
	  	df
	}

	distance_matrix_euclidean <- function() {
		print("Computing euclidean distance matrix (this can take a while)")
		unit_scale.both <- normalize.df(features.both,features)
		unit_scale.cases <- as.data.frame(unit_scale.both[1:nrow(features.cases),])
		unit_scale.controls <- as.data.frame(unit_scale.both[(nrow(features.cases)+1):nrow(features.both),])


		sapply(1:nrow(cases),function(x) { 
  			sapply(1:nrow(controls), function(y) {
				v1 <- unit_scale.cases[x,]
				v2 <- unit_scale.controls[y,]
    			dist(t(matrix(c(v1,v2),length(v1),2)))
  			})
		})
	}

	#####################################################
	# Distance calculation using mahalanobis distances
	#####################################################

	mahalanobis_dist <- function(x,y,inv_cov) {
		d <- as.matrix(x-y)
		(d %*% inv_cov) %*% t(d)
	}

	distance_matrix_mahalanobis <- function() {
		print("Computing mahalanobis distance matrix (this can take a while)")
		inv_cov <- solve(cov(features.both))

		sapply(1:nrow(cases),function(x) { 
  			sapply(1:nrow(controls), function(y) {
					v1 <- features.cases[x,]
					v2 <- features.controls[y,]
					mahalanobis_dist(v1,v2,inv_cov)
  			})
		})
	}

	############################################################
	# Nearest neighbor matching (random order) 
	############################################################
	match.nn <- function(distance_matrix) {
		taken <- c()
		assignment <- list()

		for(case in sample(1:ncol(distance_matrix))) {
			min <- 1000000
			select <- NA
			for(control in 1:nrow(distance_matrix)) {
				if (control %in% taken)
					next
				if (distance_matrix[control,case] < min) {
					min <- distance_matrix[control,case]
					assignment[[case]] <- control
				}
			}
			taken <- c(taken,assignment[[case]])
		}
		assignment
	}


	############################################################
	# Optimal matching (using lpSolve package) 
	############################################################

	match.opt <- function(distance_matrix) {
		print("matching with lpsolve")
		# FIXME: Make distance matrix quadratic with large cost padding..

		padding <- data.frame(matrix(rep(1000000,(nrow(distance_matrix) - ncol(distance_matrix))*nrow(distance_matrix)), nrow(distance_matrix), nrow(distance_matrix) - ncol(distance_matrix)))

		cost.mat <- as.matrix(cbind(distance_matrix,padding))
		
		solution <- lp.assign(cost.mat)$solution 

		assignment <- list()
		for(case in 1:ncol(distance_matrix)) {
			sol_idx <- as.vector(as.matrix(solution[,case]))
			assignment[[case]] <- which(sol_idx > 0)
		}
		print("done")

		assignment
	}

	############################################################
	# Stable-marriage matching (modified Gale-Shapley algorithm)
	############################################################
	match.gsa <- function(distance_matrix) {
		print("Gale-Shapley matching...")

		# Create preference lists
		case_preferences <- list()
		for(case in 1:ncol(distance_matrix)) {
  			preftab <- as.data.frame(distance_matrix[,case])
  			colnames(preftab) <- c("dist")
  			preftab$control_id <- 1:nrow(preftab)
  			preflist <- (preftab[with(preftab,order(dist)),])$control_id
  			case_preferences[[length(case_preferences)+1]] <- preflist
		}

		control_preferences <- list()
		for(control in 1:nrow(distance_matrix)) {
  			preftab <- as.data.frame(distance_matrix[control,])
  			colnames(preftab) <- c("dist")
  			preftab$case_id <- 1:nrow(preftab)
  			preflist <- (preftab[with(preftab,order(dist)),])$case_id
  			control_preferences[[length(control_preferences)+1]] <- preflist
		}

		unmatched_cases <- 1:(nrow(cases)*controls_per_case)
		unmatched_controls <- 1:nrow(controls)

		engagement_case <- list()
		for(i in unmatched_cases)
  	  	  engagement_case[[i]] <- 0

		engagement_control <- list()
		for(i in unmatched_controls)
  	  	  engagement_control[[i]] <- 0

		engaged <- function(x,y) { print(paste(x,paste("engaged to",y))) }

		# Adjusted indexing
		# This allows for multiple number of controls (fixed number per case)
		index <- function(idx) {
  	  	  if (idx %% nrow(cases) == 0)
    		nrow(cases)
  	  	  else
    		idx %% nrow(cases)
		}

		repeat {
  	  	  if(length(unmatched_cases) == 0)
    		break;
  	  	  match_case <- unmatched_cases[1]
  	  	  unmatched_cases <- unmatched_cases[-1]

  	  	  case_prefs <- case_preferences[[index(match_case)]]

  	  	  repeat { 
    		favourite <- case_prefs[1]
    		case_prefs <- case_prefs[-1]
    		
    		if (favourite %in% unmatched_controls) {
      	  	  engagement_case[[match_case]] <- favourite
      	  	  engagement_control[[favourite]] <- match_case
      	  	  unmatched_controls <- unmatched_controls[unmatched_controls != favourite]
      	  	  engaged(match_case,favourite)
      	  	  break
    		} else {
      	  	  # Does the control prefer this case to the one it is now engaged to?
      	  	  current_engagement <- match(index(engagement_control[[favourite]]),
                                        		control_preferences[[favourite]])

      	  	  tentive_engagement <- match(index(match_case),
                                        		control_preferences[[favourite]])

      	  	  if (tentive_engagement < current_engagement) {
        		engagement_case[[match_case]] <- favourite
        		now_single <- engagement_control[[favourite]]
        		engagement_control[[favourite]] <- match_case
        		unmatched_cases <- c(now_single,unmatched_cases)
        		engaged(match_case,favourite)        
        		break
      	  	  }
    		}
  	  	  }
		}
		print("finished matching!")
		engagement_case
	}

	############################################################
	# Caliper matching 
	############################################################
	match.caliper <- function(distance_matrix, assignment, caliper=10000) {
		for(case in sample(1:ncol(distance_matrix))) {
			for(control in 1:nrow(distance_matrix)) {
				if (distance_matrix[case,assignment[[case]]] > caliper) {
					assignment[[case]] <- NA
				}
			}
		}
		assignment
	}

	#########################################################
	# Main 
	#########################################################

	# Extract relevant features as separate dataframes
	features.cases <- subset_named_columns(cases,features)
	features.controls <- subset_named_columns(controls,features)
	features.both <- rbind(features.cases, features.controls) 

	# Compute distance matrix
	if (distance.measure=="mahalanobis") {
		if (length(features) >= 2) {
			distance_matrix <- distance_matrix_mahalanobis()
		} else {
			distance_matrix <- distance_matrix_euclidean()
		}
	} else if (distance.measure=="euclidean") {
		distance_matrix <- distance_matrix_euclidean()
	} else {
		stop(paste("Unknown distance measure:",distance.measure))
	}

	algorithm.time <- system.time({
		if (match.method=="gsa") 
			matchings <- match.gsa(distance_matrix)
		else if (match.method=="opt")
			matchings <- match.opt(distance_matrix)
		else if (match.method=="nn")
			matchings <- match.nn(distance_matrix)
		else
			stop(paste("Unknown match.method:",match.method))
	})

	# Subsequently do caliper matching
	#if(!is.na(caliper)) {
	#	matchings <- matchings
		#matchings <- match.caliper(distance_matrix,matchings,caliper)
	#}

	case_control_matches <- data.frame(matrix(as.vector(matchings),nrow(cases),controls_per_case))
	case_control_matches <- cbind(rownames(case_control_matches),case_control_matches) 
	controls_names <- sapply(1:controls_per_case, function(x) { paste0("control",x) } )
	names(case_control_matches) <- c("case",controls_names)

	# Make a data-frame for reporting differences
	differences <- data.frame(matrix(rep(0,length(features)*controls_per_case*nrow(cases)),nrow(cases),length(features)*controls_per_case))
	names(differences) <- sapply(controls_names, function(n) { sapply(features, function(f) { paste0(n,".",f) }) }) 

	print(case_control_matches)

	idx <- 1
	for(ctrl in 1:controls_per_case) {
		ctrl_ids <- unlist(case_control_matches[,1+ctrl])
		print(ctrl_ids)
		ctrl.feat <- as.data.frame(features.controls[ctrl_ids,])
		names(ctrl.feat) <- names(features.controls)
		print(str(ctrl.feat))
		print(str(features.cases))
		for(feature in features) {
			values.cases <- features.cases[,which(colnames(features.cases) == feature)]
			print(values.cases)
			values.controls <- ctrl.feat[,which(colnames(ctrl.feat) == feature)]
			print(values.controls)
			differences[,idx] <- values.cases - values.controls
			idx <- idx + 1
		} 
	}


	# Run statistical tests for differences
	test.metrics <- function() {
		test <- data.frame(matrix(rep(0,length(features)*controls_per_case),length(features),controls_per_case))
		names(test) <- controls_names
		rownames(test) <- features
		test
	}

# Calculate test statistics for result
table.ttest <- test.metrics()
	table.wilcox <- test.metrics()

	for(ctrl in 1:controls_per_case) {
		ctrl_ids <- unlist(case_control_matches[,1+ctrl])
		ctrl.feat <- as.data.frame(features.controls[ctrl_ids,])
		names(ctrl.feat) <- names(features.controls)

		for(feature in features) {
			values.cases <- features.cases[,which(colnames(features.cases) == feature)]
			values.controls <- ctrl.feat[,which(colnames(ctrl.feat) == feature)]
			wilcox <- wilcox.test(values.cases,values.controls,paired=T)
			ttest <- t.test(values.cases,values.controls,paired=T)
			table.ttest[feature,ctrl] <- ttest$p.value
			table.wilcox[feature,ctrl] <- wilcox$p.value
		} 
	}

	# Make a distances data frame 
	distances <- matrix(rep(0,nrow(cases)*controls_per_case), nrow(cases), controls_per_case)

	for(i in 1:nrow(cases))
		for(j in 1:controls_per_case) {
			control <- unlist(case_control_matches[i,j+1])
			distances[i,j] <- distance_matrix[control,i]
		}

	max.distance = max(distances,na.rm=T) 
	mean.distance = mean(distances,na.rm=T)
	median.distance = median(distances,na.rm=T)
	sd.distance = sd(distances,na.rm=T)
	unmatched = sum(is.na(distances))

	distances <- data.frame(distances)
	
	names(distances) <- controls_names  

	# Use named case and controls if specified 
	if (!is.na(id.column)) {
		rownames(differences) <- cases[,id.column]
		rownames(distances) <- cases[,id.column]
		case_control_matches$case <- cases[,id.column]
		for(ctrl in 1:controls_per_case) 
			case_control_matches[,ctrl+1] <- sapply(as.vector(case_control_matches[,ctrl+1]), function(x) { controls[x,id.column] })
	}

	# Return value is the list of matches associated tables and metrics
	list(matching=case_control_matches,
		 distances=distances,
		 differences=differences,
		 test.metrics=list(
				t.test=table.ttest,
				wilcox=table.wilcox,
				longest.distance = max.distance,
				mean.distance = mean.distance,
				median.distance = median.distance,
				sd.distance = sd.distance, 
				elapsed.time = algorithm.time[3],
				unmatched = unmatched
		)
	)
}

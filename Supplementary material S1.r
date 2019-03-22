# Supplementary material S1. R script of the model
# A framework for modeling the effects of global warming on the annual distribution of Lyme borreliosis incidences. International Journal of Global Warming


# example codes:
# temperature_reference <- c(0.525584416, 0.968051948, 1.212337662, -0.253506494, 1.268311688, 2.850909091, 2.368441558, 3.6, 4.442467532, 5.247012987, 6.846883117, 6.840909091, 9.565194805, 10.84272727, 11.22155844, 13.02675325, 15.01051948, 16.47974026, 17.08727273, 17.12324675, 18.25545455, 19.46337662, 20.19649351, 21.11896104, 21.6525974, 21.94922078, 21.93792208, 22.11051948, 23.11506494, 23.34194805, 23.14233766, 22.03246753, 22.98857143, 22.15090909, 19.44142857, 18.54480519, 17.44441558, 16.10922078, 15.52792208, 15.36818182, 13.18675325, 10.67922078, 11.31467532, 9.993376623, 7.455844156, 6.307272727, 4.848831169, 3.842467532, 3.107662338, 0.614545455, -0.035974026, -1.518181818)
# temperature_future <- c(4.143121517, 2.533840788, 5.02514879, 4.98937175, 6.341839651, 7.347490081, 7.563740551, 8.864998364, 9.798232517, 10.13040328, 12.34854132, 13.072441, 14.47406425, 15.56969005, 18.32036731, 18.63026555, 19.02474706, 19.64323935, 21.88566403, 21.60486091, 22.31980831, 23.46785052, 24.08197117, 25.32103779, 25.43385909, 25.06923442, 26.53963844, 26.87934416, 25.76484221, 25.25323961, 26.11458727, 24.13661623, 22.27884857, 20.6605821, 19.72977626, 18.71996009, 18.72795275, 17.85945584, 15.39848206, 13.68586849, 14.67704456, 12.11088471, 10.31582156, 8.988419002, 8.340396404, 7.032454971, 5.100380457, 4.525348073, 1.585883952, 2.852230669, 6.94218, 7.049697143)
# plot(Lyme_model(T = temperature_reference)$LB, xlab = "week number", ylab = "relative Lyme borreliosis incidence")
# prediction_reference <- data.frame(Lyme_model(T = temperature_reference))
# prediction_future <- data.frame(Lyme_model(T = temperature_future))
# plot(prediction_reference$BA, ylim = c(0, 8), xlab = "week number", ylab = "relative biting activity and Lyme borreliosis incidence")
# points(prediction_future$BA, pch = 4)
# points(prediction_reference$LB, col = "red")
# points(prediction_future$LB, col = "red", pch = 4)
# lines(prediction_reference$BA)
# lines(prediction_future$BA)
# lines(prediction_reference$LB, col = "red")
# lines(prediction_future$LB, col = "red")

# Inputs (arguments):
# 	T: numeric vector with 52 or n * 52 elements, weekly mean temperature values (°C), starting from the 1st week of January
# 	HM: numeric vector with 52 elements or NULL/NA, holiday multiplier of the weeks. If NULL/NA, default is used.
# 	delta: numeric vector with one element (>=0) or NULL/NA, weight parameter of Model I., season 1. If NULL/NA, default is used.
# 	weight_factors: numeric vector with 52 elements or NULL/NA, weight factors of Model II, starting from the -52th week, ending with the -1st week. If NULL/NA, default (model C) is used.
# 	alpha1, alpha2: numeric vectors with one element or NULL/NA, alpha (axis) parameter of Model I., season 1 or 2. If NULL/NA, default is used.
# 	mu1, mu2: numeric vectors with one element or NULL/NA, mu (mean) parameter of Model I., season 1 or 2. If NULL/NA, default is used.
# 	sigma1, sigma2: numeric vectors with one element (>=0) or NULL/NA, sigma (standard deviation) parameter of Model I., season 1 or 2. If NULL/NA, default is used.
# 	c1, c2: numeric vectors with one element or NULL/NA, c (multiplier) parameter of Model I., season 1 or 2. If NULL/NA, default is used.
# 	d1, d2: numeric vectors with one element or NULL/NA, d (difference) parameter of Model I., season 1 or 2. If NULL/NA, default is used.
# Output (returned value): named list of two vectors
# 	BA: numeric vector of the same length as that of T, weekly relative tick biting activity, % of the annual activity, starting from the 1st week of January
# 	LB: numeric vector of the same length as that of T, weekly relative Lyme borreliosis incidences, % of the annual incidence, starting from the 1st week of January
Lyme_model <- function(T, HM, delta, weight_factors, alpha1, alpha2, mu1, mu2, sigma1, sigma2, c1, c2, d1, d2) {
	number_of_weeks <- 52 # number of weeks per one year
	
	# Checking argument T
	if (missing(T)) { # if T is not set
		stop("Argument T should be set.") # throw error and stop the function
	} # if
	if (identical(NA, T)) { # if T is NA
		stop("Argument T must not be NA.") # throw error and stop the function
	} # if
	if (!(length(T) %% number_of_weeks == 0)) { # if length of T is not the multiple of 52
		stop("Length of argument T should be 52 or n*52.") # throw error and stop the function
	} # if
	number_of_years <- length(T) / number_of_weeks # number of years calculated from the length of T
	if (!(is.numeric(T))) { # if T is not numeric vector
		stop("Argument T should contain numeric values.") # throw error and stop the function
	} # if
	if (any(is.na(T))) { # if T contains NA value(s)
		stop("Argument T must not contain NA values.") # throw error and stop the function
	} # if
	
	# Checking argument HM
	if (missing(HM)) { # if HM is not set
		use_default_HM <- TRUE # logical flag set to TRUE
	} else if (identical(NA, HM)) { # if HM is NA
		use_default_HM <- TRUE # logical flag set to TRUE
	} else { # if HM is interpretable
		use_default_HM <- FALSE # logical flag set to FALSE
	} # else
	if (use_default_HM) { # if the default holiday multiplier should be used
		HM <- c(rep(x = 1, times = 24), 1.29, 1.59, 1.91, 2.24, 2.36, 2.49, 2.74, 2.99, 2.61, 2.22, 1.65, 1.08, rep(x = 1, times = number_of_weeks - 36)) # use default holiday multiplier values
	} else { # if argument HM should be interpreted
		if (length(HM) != number_of_weeks) { # if length of HM is not 52
			stop("Length of argument HM should be 52.") # throw error and stop the function
		} # if
		if (!(is.numeric(HM))) { # if HM is not numeric vector
			stop("Argument HM should contain numeric values.") # throw error and stop the function
		} # if
		if (any(is.na(HM))) { # if HM contains NA value(s)
			HM[is.na(HM)] <- 1 # impute 1 where HM is NA
			warning("Argument HM contained NA values. Ones were imputed.") # throw warning
		} # if
	} # else
	
	# Checking argument delta
	if (missing(delta)) { # if delta is not set
		use_default_delta <- TRUE # logical flag set to TRUE
	} else if (identical(NA, delta)) { # if delta is NA
		use_default_delta <- TRUE # logical flag set to TRUE
	} else { # if delta is interpretable
		use_default_delta <- FALSE # logical flag set to FALSE
	} # else
	if (use_default_delta) { # if the default weight parameter should be used
		delta <- 0.0078 # set the weight parameter to default
	} else { # if argument delta should be interpreted
		if (!(is.numeric(delta))) { # if delta is not numeric vector
			stop("Argument delta should contain a numeric value.") # throw error and stop the function
		} # if
		if (length(delta) != 1) { # if argument delta contains 0 or >1 values
			stop("Argument delta should contain exactly one numeric value.") # throw error and stop the function
		} # if
		if (delta < 0) { # if argument delta is lower than zero
			stop("Argument delta must not be lower than zero.") # throw error and stop the function
		} # if
	} # else
	
	# Checking argument weight_factors
	if (missing(weight_factors)) { # if weight_factors is not set
		use_default_weight_factors <- TRUE # logical flag set to TRUE
	} else if (identical(NA, weight_factors)) { # if weight_factors is NA
		use_default_weight_factors <- TRUE # logical flag set to TRUE
	} else { # if weight_factors is interpretable
		use_default_weight_factors <- FALSE # logical flag set to FALSE
	} # else
	if (use_default_weight_factors) { # if the default weight factors should be used
		weight_factors <- c(rep(x = 0, times = 38), 0.0025, 0.0029, 0, 0, 0.0991, 0, 0.0001, 0.1328, 0.0403, 0.1130, 0.0514, 0.1476, 0.1813, 0.2213) # use default weight factors as defined by Model C
	} else { # if argument weight_factors should be interpreted
		if (length(weight_factors) != number_of_weeks) { # if length of weight_factors is not 52
			stop("Length of argument weight_factors should be 52.") # throw error and stop the function
		} # if
		if (!(is.numeric(weight_factors))) { # if weight_factors is not numeric vector
			stop("Argument weight_factors should contain numeric values.") # throw error and stop the function
		} # if
		if (any(is.na(weight_factors))) { # if HM contains NA value(s)
			weight_factors[is.na(weight_factors)] <- 0 # impute 0 where HM is NA
			warning("Argument weight_factors contained NA values. Zeros were imputed.") # throw warning
		} # if
	} # else
	
	# Checking argument alpha1
	if (missing(alpha1)) { # if alpha1 is not set
		use_default_alpha1 <- TRUE # logical flag set to TRUE
	} else if (identical(NA, alpha1)) { # if alpha1 is NA
		use_default_alpha1 <- TRUE # logical flag set to TRUE
	} else { # if alpha1 is interpretable
		use_default_alpha1 <- FALSE # logical flag set to FALSE
	} # else
	if (use_default_alpha1) { # if the alpha1 parameter should be used
		alpha1 <- 26.3302 # set the alpha1 parameter to default
	} else { # if argument alpha1 should be interpreted
		if (!(is.numeric(alpha1))) { # if alpha1 is not numeric vector
			stop("Argument alpha1 should contain a numeric value.") # throw error and stop the function
		} # if
		if (length(alpha1) != 1) { # if argument alpha1 contains 0 or >1 values
			stop("Argument alpha1 should contain exactly one numeric value.") # throw error and stop the function
		} # if
	} # else
	
	# Checking argument alpha2
	if (missing(alpha2)) { # if alpha2 is not set
		use_default_alpha2 <- TRUE # logical flag set to TRUE
	} else if (identical(NA, alpha2)) { # if alpha2 is NA
		use_default_alpha2 <- TRUE # logical flag set to TRUE
	} else { # if alpha2 is interpretable
		use_default_alpha2 <- FALSE # logical flag set to FALSE
	} # else
	if (use_default_alpha2) { # if the alpha2 parameter should be used
		alpha2 <- 123.8382 # set the alpha2 parameter to default
	} else { # if argument alpha2 should be interpreted
		if (!(is.numeric(alpha2))) { # if alpha2 is not numeric vector
			stop("Argument alpha2 should contain a numeric value.") # throw error and stop the function
		} # if
		if (length(alpha2) != 1) { # if argument alpha2 contains 0 or >1 values
			stop("Argument alpha2 should contain exactly one numeric value.") # throw error and stop the function
		} # if
	} # else
	
	# Checking argument mu1
	if (missing(mu1)) { # if mu1 is not set
		use_default_mu1 <- TRUE # logical flag set to TRUE
	} else if (identical(NA, mu1)) { # if mu1 is NA
		use_default_mu1 <- TRUE # logical flag set to TRUE
	} else { # if mu1 is interpretable
		use_default_mu1 <- FALSE # logical flag set to FALSE
	} # else
	if (use_default_mu1) { # if the mu1 parameter should be used
		mu1 <- 2.0337 # set the mu1 parameter to default
	} else { # if argument mu1 should be interpreted
		if (!(is.numeric(mu1))) { # if mu1 is not numeric vector
			stop("Argument mu1 should contain a numeric value.") # throw error and stop the function
		} # if
		if (length(mu1) != 1) { # if argument mu1 contains 0 or >1 values
			stop("Argument mu1 should contain exactly one numeric value.") # throw error and stop the function
		} # if
	} # else
	
	# Checking argument mu2
	if (missing(mu2)) { # if mu2 is not set
		use_default_mu2 <- TRUE # logical flag set to TRUE
	} else if (identical(NA, mu2)) { # if mu2 is NA
		use_default_mu2 <- TRUE # logical flag set to TRUE
	} else { # if mu2 is interpretable
		use_default_mu2 <- FALSE # logical flag set to FALSE
	} # else
	if (use_default_mu2) { # if the mu2 parameter should be used
		mu2 <- 4.7124 # set the mu2 parameter to default
	} else { # if argument mu2 should be interpreted
		if (!(is.numeric(mu2))) { # if mu2 is not numeric vector
			stop("Argument mu2 should contain a numeric value.") # throw error and stop the function
		} # if
		if (length(mu2) != 1) { # if argument mu2 contains 0 or >1 values
			stop("Argument mu2 should contain exactly one numeric value.") # throw error and stop the function
		} # if
	} # else
	
	# Checking argument sigma1
	if (missing(sigma1)) { # if sigma1 is not set
		use_default_sigma1 <- TRUE # logical flag set to TRUE
	} else if (identical(NA, sigma1)) { # if sigma1 is NA
		use_default_sigma1 <- TRUE # logical flag set to TRUE
	} else { # if sigma1 is interpretable
		use_default_sigma1 <- FALSE # logical flag set to FALSE
	} # else
	if (use_default_sigma1) { # if the sigma1 parameter should be used
		sigma1 <- 0.4804 # set the sigma1 parameter to default
	} else { # if argument sigma1 should be interpreted
		if (!(is.numeric(sigma1))) { # if sigma1 is not numeric vector
			stop("Argument sigma1 should contain a numeric value.") # throw error and stop the function
		} # if
		if (length(sigma1) != 1) { # if argument sigma1 contains 0 or >1 values
			stop("Argument sigma1 should contain exactly one numeric value.") # throw error and stop the function
		} # if
		if (sigma1 < 0) { # if argument sigma1 is lower than zero
			stop("Argument sigma1 must not be lower than zero.") # throw error and stop the function
		} # if
	} # else
	
	# Checking argument sigma2
	if (missing(sigma2)) { # if sigma2 is not set
		use_default_sigma2 <- TRUE # logical flag set to TRUE
	} else if (identical(NA, sigma2)) { # if sigma2 is NA
		use_default_sigma2 <- TRUE # logical flag set to TRUE
	} else { # if sigma2 is interpretable
		use_default_sigma2 <- FALSE # logical flag set to FALSE
	} # else
	if (use_default_sigma2) { # if the sigma2 parameter should be used
		sigma2 <- 0.0145 # set the sigma2 parameter to default
	} else { # if argument sigma2 should be interpreted
		if (!(is.numeric(sigma2))) { # if sigma2 is not numeric vector
			stop("Argument sigma2 should contain a numeric value.") # throw error and stop the function
		} # if
		if (length(sigma2) != 1) { # if argument sigma2 contains 0 or >1 values
			stop("Argument sigma2 should contain exactly one numeric value.") # throw error and stop the function
		} # if
		if (sigma2 < 0) { # if argument sigma2 is lower than zero
			stop("Argument sigma2 must not be lower than zero.") # throw error and stop the function
		} # if
	} # else
	
	# Checking argument c1
	if (missing(c1)) { # if c1 is not set
		use_default_c1 <- TRUE # logical flag set to TRUE
	} else if (identical(NA, c1)) { # if c1 is NA
		use_default_c1 <- TRUE # logical flag set to TRUE
	} else { # if c1 is interpretable
		use_default_c1 <- FALSE # logical flag set to FALSE
	} # else
	if (use_default_c1) { # if the c1 parameter should be used
		c1 <- 82.7165 # set the c1 parameter to default
	} else { # if argument c1 should be interpreted
		if (!(is.numeric(c1))) { # if c1 is not numeric vector
			stop("Argument c1 should contain a numeric value.") # throw error and stop the function
		} # if
		if (length(c1) != 1) { # if argument c1 contains 0 or >1 values
			stop("Argument c1 should contain exactly one numeric value.") # throw error and stop the function
		} # if
	} # else
	
	# Checking argument c2
	if (missing(c2)) { # if c2 is not set
		use_default_c2 <- TRUE # logical flag set to TRUE
	} else if (identical(NA, c2)) { # if c2 is NA
		use_default_c2 <- TRUE # logical flag set to TRUE
	} else { # if c2 is interpretable
		use_default_c2 <- FALSE # logical flag set to FALSE
	} # else
	if (use_default_c2) { # if the c2 parameter should be used
		c2 <- 6.8409 # set the c2 parameter to default
	} else { # if argument c2 should be interpreted
		if (!(is.numeric(c2))) { # if c2 is not numeric vector
			stop("Argument c2 should contain a numeric value.") # throw error and stop the function
		} # if
		if (length(c2) != 1) { # if argument c2 contains 0 or >1 values
			stop("Argument c2 should contain exactly one numeric value.") # throw error and stop the function
		} # if
	} # else
	
	# Checking argument d1
	if (missing(d1)) { # if d1 is not set
		use_default_d1 <- TRUE # logical flag set to TRUE
	} else if (identical(NA, d1)) { # if d1 is NA
		use_default_d1 <- TRUE # logical flag set to TRUE
	} else { # if d1 is interpretable
		use_default_d1 <- FALSE # logical flag set to FALSE
	} # else
	if (use_default_d1) { # if the d1 parameter should be used
		d1 <- 0L # set the d1 parameter to default
	} else { # if argument d1 should be interpreted
		if (!(is.numeric(d1))) { # if d1 is not numeric vector
			stop("Argument d1 should contain a numeric value.") # throw error and stop the function
		} # if
		if (length(d1) != 1) { # if argument d1 contains 0 or >1 values
			stop("Argument d1 should contain exactly one numeric value.") # throw error and stop the function
		} # if
	} # else
	
	# Checking argument d2
	if (missing(d2)) { # if d2 is not set
		use_default_d2 <- TRUE # logical flag set to TRUE
	} else if (identical(NA, d2)) { # if d2 is NA
		use_default_d2 <- TRUE # logical flag set to TRUE
	} else { # if d2 is interpretable
		use_default_d2 <- FALSE # logical flag set to FALSE
	} # else
	if (use_default_d2) { # if the d2 parameter should be used
		d2 <- 0.4931 # set the d2 parameter to default
	} else { # if argument d2 should be interpreted
		if (!(is.numeric(d2))) { # if d2 is not numeric vector
			stop("Argument d2 should contain a numeric value.") # throw error and stop the function
		} # if
		if (length(d2) != 1) { # if argument d2 contains 0 or >1 values
			stop("Argument d2 should contain exactly one numeric value.") # throw error and stop the function
		} # if
	} # else
	
	# Calculation
	S <- P <- TDA <- A1 <- A2 <- A <- BA <- LB <- rep(x = NA, times = number_of_years * number_of_weeks) # create initial (empty) vectors
	for (year in 1:number_of_years) { # iterate through the years
		max_T_of_year <- max(T[((year - 1) * number_of_weeks + 1):(year * number_of_weeks)]) # maximum temperature of the studied year
		for (week_of_year in 1:number_of_weeks) { # iterate through the weeks of the year
			week_of_studied_period <- ((year - 1) * number_of_weeks) + week_of_year # week number of the studied period (index of vectors)
			
			# Equation 6:
			if (T[week_of_studied_period] >= alpha1 | T[week_of_studied_period] <= 5 | ifelse(week_of_studied_period == 1, FALSE, T[week_of_studied_period - 1] <= 5)) { # in case of too low or too high temperature
				TDA[week_of_studied_period] <- 0 # temperature-dependent activity is set to zero
			} else { # in case of optimal temperature
				TDA[week_of_studied_period] <- c1 * 1 / ((alpha1 - T[week_of_studied_period]) * sqrt(2 * pi) * sigma1) * exp(1) ^ ((-1) * ((log(alpha1 - T[week_of_studied_period]) - mu1) ^ 2) / (2 * sigma1 ^ 2)) + d1 # calculate temperature-dependent activity
			} # else
			
			# Equation 5:
			S[week_of_studied_period] <- TDA[week_of_studied_period] * delta # calculate the subtrahend as the product of temperature dependent activity and the weight parameter
			
			# Equation 3:
			if (week_of_year == 1) { # first week of the year
				P[week_of_studied_period] <- 1 # active population is set to 1
			
			# Equation 4:
			} else if (P[week_of_studied_period - 1] - S[week_of_studied_period - 1] < 0) { # other weeks, if there is no more active population
				P[week_of_studied_period] <- 0 # active population is set to 0
			} else { # other weeks, if there is still active population
				P[week_of_studied_period] <- P[week_of_studied_period - 1] - S[week_of_studied_period - 1] # calculate active population
			} # else
			
			# Equation 7:
			A1[week_of_studied_period] <- P[week_of_studied_period] * TDA[week_of_studied_period] # calculate tick activity of season 1
			
			# Equation 8:
			if (T[week_of_studied_period] >= alpha2 | T[week_of_studied_period] <= 5 | !(max_T_of_year %in% T[((year - 1) * number_of_weeks + 1):week_of_studied_period]) | week_of_year <= 28 || all(T[((year - 1) * number_of_weeks + 29):week_of_studied_period] >= 20)) { # in case of season2
				A2[week_of_studied_period] <- 0 # tick activity is set to zero
			} else { # in case of out of season2
				A2[week_of_studied_period] <- c2 * 1 / ((alpha2 - T[week_of_studied_period]) * sqrt(2 * pi) * sigma2) * exp(1) ^ ((-1) * ((log(alpha2 - T[week_of_studied_period]) - mu2) ^ 2) / (2 * sigma2 ^ 2)) + d2 # calculate tick activity of season 2
			} # else
			
			# Equation 1:
			A[week_of_studied_period] <- A1[week_of_studied_period] + A2[week_of_studied_period] # calculate tick activity as the sum of tick activity of season 1 and season 2
			
			# Equation 2:
			BA[week_of_studied_period] <- A[week_of_studied_period] * HM[week_of_year] # calculate biting activity as the product of tick activity and holiday multiplier
		} # for week_of_year
	} # for year
	BA_reordered <- c(BA[2:52], BA, BA[1]) # create a copy of the first year (to formally be the previous year), and shift with one week
	for (year in 1:number_of_years) { # iterate through the years
		for (week_of_year in 1:number_of_weeks) { # iterate through the weeks of the year
			
			# Equation 9:
			LB[((year - 1) * number_of_weeks) + week_of_year] <- sum(BA_reordered[(((year - 1) * number_of_weeks) + week_of_year):((year * number_of_weeks) + week_of_year - 1)] * weight_factors) # calculate relative Lyme borreliosis incidences
		} # for week_of_year
	} # for year
	return(list(BA = BA, LB = LB)) # return the result
} # Lyme_model()


qtrapezoid <- function(p, min = 0, mode1 = 1/3, mode2 = 2/3, max = 1, n1 = 2, 
		n3 = 2, alpha = 1, lower.tail = TRUE, log.p = FALSE)
{
	# Test parameters are well-formed
	max.length <- max(length(p), length(min), length(mode1), length(mode2), 
			length(max), length(n1), length(n3), length(alpha))
	
	if(log.p) p <- exp(p)
	if(!lower.tail) p <- 1-p
	
	parameters <- matrix(nrow = max.length, ncol = 8)
	tryCatch(
			{
				parameters[,1] <- p
				parameters[,2] <- min
				parameters[,3] <- mode1
				parameters[,4] <- mode2
				parameters[,5] <- max
				parameters[,6] <- n1
				parameters[,7] <- n3
				parameters[,8] <- alpha
			}, error = function(X) {
				stop(paste(" -- Aruments min, mode1, mode2, max, n1, n3, and alpha
										must be of equal length", sep=""))
			})
	
	# Function to calculate quantiles
	qFunction <- function(p, a, b, c, d, m, n, alpha)
	{
		# Normalizing constant
		normalizing.constant <- 
				(2 * m * n) /
				((2 * alpha * (b - a) * n) + 
					((alpha + 1) * (c - b) * m * n) +
					(2 * (d - c) * m))
		
		# Calculate pi1, pi2, pi3
		pi1 <- (2 * alpha * (b - a) * n) /
				((2 * alpha * (b - a) * n) + ((alpha + 1) * (c - b) * m * n) +
					(2 * (d - c) * m))
		pi2 <- ((alpha + 1) * (c - b) * m * n) /
				((2 * alpha * (b - a) * n) +
					((alpha + 1) * (c - b) * m * n) +
					(2 * (d - c) * m))
		pi3 <- (2 * (d - c) * m) / 
				((2 * alpha * (b - a) * n) +
					((alpha + 1) * (c - b) * m * n) +
					(2 * (d - c) * m))
		
		# Calculate and return quantile conditional on p
		if (0 <= p & p <= pi1)
		{
			out <- ((p*(2*alpha*(b-a)*n+(alpha+1)*(c-b)*m*n+2*(d-c)*m))/(2*alpha*(b-a)*n))^(1/m)*(b-a)+a
		} else if (pi1 < p & p <= 1 - pi3 & alpha != 1)
		{
			#out <- 2*c-b-(((p*(2*alpha*(b-a)*n+(alpha+1)*(c-b)*m*n+2*(d-c)*m)-2*alpha*(b-a)*n)/(2*(q-b)*m*n)-1)*2*(c-b))/(alpha-1)
			#out <- b+((p*(2*alpha*(b-a)*n+(alpha+1)*(c-b)*m*n+2*(d-c)*m)-2*alpha*(b-a)*n)*2*(c-b))/(2*(2*c-b-((-2)*(c-b))/(alpha-1))*(alpha-1)*n*m)
			#out <- 2*c-b-(((p*(2*alpha*(b-a)*n+(alpha+1)*(c-b)*m*n+2*(d-c)*m)-2*alpha*(b-a)*n)/(2*(q-b)*m*n)-1)*2*(c-b))/(alpha-1)
			#out <- ((Sqrt((((-2)*b*m*n*(1-alpha))/(2*(c-b))+2*m*n*(((2*c-b)*(alpha-1))/(2*(c-b))+1))^2/(2*alpha*(b-a)*n+(alpha+1)*(c-b)*m*n+2*(d-c)*m)^2-(((2*alpha*(b-a)*n+(-2)*b*m*n*(((2*c-b)*(alpha-1))/(2*(c-b))+1))/(2*alpha*(b-a)*n+(alpha+1)*(c-b)*m*n+2*(d-c)*m)-p)*4*2*m*n*(1-alpha))/(2*(c-b)*(2*alpha*(b-a)*n+(alpha+1)*(c-b)*m*n+2*(d-c)*m)))-(((-2)*b*m*n*(1-alpha))/(2*(c-b))+2*m*n*(((2*c-b)*(alpha-1))/(2*(c-b))+1))/(2*alpha*(b-a)*n+(alpha+1)*(c-b)*m*n+2*(d-c)*m))*2*(c-b)*(2*alpha*(b-a)*n+(alpha+1)*(c-b)*m*n+2*(d-c)*m))/(2*2*m*n*(1-alpha)),q==(-((((-2)*b*m*n*(1-alpha))/(2*(c-b))+2*m*n*(((2*c-b)*(alpha-1))/(2*(c-b))+1))/(2*alpha*(b-a)*n+(alpha+1)*(c-b)*m*n+2*(d-c)*m)+Sqrt((((-2)*b*m*n*(1-alpha))/(2*(c-b))+2*m*n*(((2*c-b)*(alpha-1))/(2*(c-b))+1))^2/(2*alpha*(b-a)*n+(alpha+1)*(c-b)*m*n+2*(d-c)*m)^2-(((2*alpha*(b-a)*n+(-2)*b*m*n*(((2*c-b)*(alpha-1))/(2*(c-b))+1))/(2*alpha*(b-a)*n+(alpha+1)*(c-b)*m*n+2*(d-c)*m)-p)*4*2*m*n*(1-alpha))/(2*(c-b)*(2*alpha*(b-a)*n+(alpha+1)*(c-b)*m*n+2*(d-c)*m))))*2*(c-b)*(2*alpha*(b-a)*n+(alpha+1)*(c-b)*m*n+2*(d-c)*m))/(2*2*m*n*(1-alpha))
			out <- ((sqrt((((-2)*b*m*n*(1-alpha))/(2*(c-b))+2*m*n*(((2*c-b)*(alpha-1))/(2*(c-b))+1))^2/(2*alpha*(b-a)*n+(alpha+1)*(c-b)*m*n+2*(d-c)*m)^2-(((2*alpha*(b-a)*n+(-2)*b*m*n*(((2*c-b)*(alpha-1))/(2*(c-b))+1))/(2*alpha*(b-a)*n+(alpha+1)*(c-b)*m*n+2*(d-c)*m)-p)*4*2*m*n*(1-alpha))/(2*(c-b)*(2*alpha*(b-a)*n+(alpha+1)*(c-b)*m*n+2*(d-c)*m)))-(((-2)*b*m*n*(1-alpha))/(2*(c-b))+2*m*n*(((2*c-b)*(alpha-1))/(2*(c-b))+1))/(2*alpha*(b-a)*n+(alpha+1)*(c-b)*m*n+2*(d-c)*m))*2*(c-b)*(2*alpha*(b-a)*n+(alpha+1)*(c-b)*m*n+2*(d-c)*m))/(2*2*m*n*(1-alpha))
		} else if (pi1 < p & p <= 1 - pi3 & alpha == 1)
		{
			out <- b + (((p - pi1) / (pi2)) * (c - b))
		} else if (1 - pi3 < p & p <= 1)
		{
			out <- d-(((1-p)*(2*alpha*(b-a)*n+(alpha+1)*(c-b)*m*n+2*(d-c)*m))/(2*(d-c)*m))^(1/n)*(d-c)
		}	else
		{
			out <- NaN
		}
		
		if(any(is.nan(out))) warning("NaN in qtrapezoid")
		if(any(is.na(out))) warning("NA in qtrapezoid")
		
		return(out)
	}
	
	# Calculate quantiles
	out <- apply(parameters, 1, function(x)
			{
				qFunction(p = x[1], a = x[2], b = x[3], c = x[4], d = x[5], 
						m = x[6], n = x[7], alpha = x[8])
			})
	
	return(out)
}

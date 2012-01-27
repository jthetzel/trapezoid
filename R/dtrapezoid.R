dtrapezoid <- function(x, min = 0, mode1 = 1/3, mode2 = 2/3, max = 1, n1 = 2, 
		n3 = 2, alpha = 1, log = FALSE)
{
	# Test parameters are well-formed
	max.length <- max(length(x), length(min), length(mode1), length(mode2), 
			length(max), length(n1), length(n3), length(alpha))
	
	parameters <- matrix(nrow = max.length, ncol = 8)
	tryCatch(
			{
				parameters[,1] <- x
				parameters[,2] <- min
				parameters[,3] <- mode1
				parameters[,4] <- mode2
				parameters[,5] <- max
				parameters[,6] <- n1
				parameters[,7] <- n3
				parameters[,8] <- alpha
			}, error = function(X) {
				stop(paste(" -- Arguments min, mode1, mode2, max, n1, n3, and alpha
										must be of equal length", sep=""))
			})
	
	# Function to calculate density
	dFunction <- function(x, a, b, c, d, m, n, alpha)
	{
		normalizing.constant <- 
				(2 * m * n) /
				((2 * alpha * (b - a) * n) + 
					((alpha + 1) * (c - b) * m * n) +
					(2 * (d - c) * m))
		
		if (a <= x & x < b)
		{
			normalizing.constant * 
					alpha * ((x - a) / (b - a))^(m - 1)
		} else if (b <= x & x < c)
		{
			normalizing.constant *
					(((1 - alpha) * ((x - b) / (c - b))) + alpha)
		} else if (c <= x & x <= d)
		{
			normalizing.constant *
					((d - x) / (d - c))^(n - 1)
		} else
		{
			0
		}
	}
	
	# Calculate density
	out <- apply(parameters, 1, function(x)
			{
				dFunction(x = x[1], a = x[2], b = x[3], c = x[4], d = x[5], 
						m = x[6], n = x[7], alpha = x[8])
			})
	
	if(log) out <- log(out)
	if(any(is.nan(out))) warning("NaN in dtrapezoid")
	if(any(is.na(out))) warning("NA in dtrapezoid")
	
	return(out)
}

ptrapezoid <- function(q, min = 0, mode1 = 1/3, mode2 = 2/3, max = 1, n1 = 2, 
		n3 = 2, alpha = 1, lower.tail = TRUE, log.p = FALSE)
{
	# Test parameters are well-formed
	max.length <- max(length(q), length(min), length(mode1), length(mode2), 
			length(max), length(n1), length(n3), length(alpha))
	
	parameters <- matrix(nrow = max.length, ncol = 8)
	tryCatch(
			{
				parameters[,1] <- q
				parameters[,2] <- min
				parameters[,3] <- mode1
				parameters[,4] <- mode2
				parameters[,5] <- max
				parameters[,6] <- n1
				parameters[,7] <- n3
				parameters[,8] <- alpha
			}, error = function(X) 
			{
				stop(paste(" -- min, mode1, mode2, max, n1, n3, and alpha
										must be of equal length", sep=""))
			})
	
	# Function to calculate cumulative probability
	pFunction <- function(q, a, b, c, d, m, n, alpha)
	{
		if (q < a)
		{
			0
		} else if (a <= q & q < b)
		{
			(2 * alpha * (b - a) * n) /
					((2 * alpha * (b - a) * n) + ((alpha + 1) * (c - b) * m * n) +
						(2 * (d - c) * m)) *
					((q - a) / (b - a))^(m)
		} else if (b <= q & q < c)
		{
			((2 * alpha * (b - a) * n) + (2 * (q - b) * m * n *
							(1 + ((alpha - 1) / 2) * ((2*c - b - q) / (c - b))))) / 
					((2 * alpha * (b - a) * n) + ((alpha + 1) * (c - b) * m * 
							n) + (2 * (d - c) * m))
		} else if (c <= q & q < d)
		{
			1 - ((2 * (d - c) * m) / 
						((2 * alpha * (b - a) * n) + ((alpha + 1) * (c - b) * m * 
								n) + (2 * (d - c) * m))) * ((d - q) / (d - c))^(n)
		} else if (q >= d)
		{
			1
		}
	}
	
	# Calculate cumulative probability
	out <- apply(parameters, 1, function(x)
			{
				pFunction(q = x[1], a = x[2], b = x[3], c = x[4], d = x[5], 
						m = x[6], n = x[7], alpha = x[8])
			})
	
	out[q < min] <- 0
	out[q > max] <- 1
	if(!lower.tail) out <- 1 - out
	if(log.p) out <- log(out)
	if(any(is.nan(out))) warning("NaN in ptrapezoid")
	if(any(is.na(out))) warning("NA in ptrapezoid")
	
	return(out)
	
}

% Copyright 2011 Jeremy Thoms Hetzel
	% Distributed under GPL 3 or later\name{trapezoid}
	
	\name{trapezoid}
\alias{Trapezoid}
\alias{trapezoid}
\alias{dtrapezoid}
\alias{ptrapezoid}
\alias{qtrapezoid}
\alias{rtrapezoid}
\title{The Trapezoidal Distribution}
\description{
	Density function, distribution function, quantile function, and random
	generation for the trapezoidal distribution with minimum equal to
	\samp{min}, lower mode equal to \samp{mode1}, upper mode equal  to \samp{mode2}, and maximum 
	equal to \samp{max}. For the generalized trapezoidal distribution, \samp{n1},
	\samp{n3}, and \samp{alpha} may optionally be specified.
}
\usage{
	dtrapezoid(x, min = 0, mode1 = 1/3, mode2 = 2/3, max = 1, n1 = 2, n3 = 2, 
		alpha = 1, log = FALSE)
	ptrapezoid(q, min = 0, mode1 = 1/3, mode2 = 2/3, max = 1, n1 = 2, n3 = 2, 
		alpha = 1, lower.tail = TRUE, log.p = FALSE)
	qtrapezoid(p, min = 0, mode1 = 1/3, mode2 = 2/3, max = 1, n1 = 2, n3 = 2, 
		alpha = 1, lower.tail = TRUE, log.p = FALSE)
	rtrapezoid(n, min = 0, mode1 = 1/3, mode2 = 2/3, max = 1, n1 = 2, n3 = 2, 
		alpha = 1)
}
\arguments{
	\item{x,q}{vector of quantiles.}
	\item{p}{vector of probabilities.}
	\item{n}{number of observations. If length(n) \eqn{>} 1, the length
		is taken to be the number required.}
	\item{min}{vector of minima.}
	\item{mode1}{vector of lower modes.}
	\item{mode2}{vector of upper modes.}
	\item{max}{vector of maxima.}
	\item{n1}{vector of growth parameters.}
	\item{n3}{vector of decay parameters.}
	\item{alpha}{vector of boundary ratio parameters.}
	\item{log, log.p}{logical; if \samp{TRUE}, probabilities \samp{p} are
		given as \samp{log(p)}.}
	\item{lower.tail}{logical; if \samp{TRUE} (default), probabilities
		are \samp{P[X <= x]}, otherwise, \samp{P[X > x]}.}
}
\value{
	\samp{dtrapezoid} gives the density function, \samp{ptrapezoid} gives the
	distribution function, \samp{qtrapezoid} gives the quantile function,
	and \samp{rtrapezoid} generates random deviates.
}
\details{
	The generalized trapezoidal distribution is described by van Dorp and Kotz (2003)
	and van Dorp and colleagues (2007). 
	With \samp{n1}, \samp{n3}, and \samp{alpha} equal to the default values of 
	2, 2, and 1, respectively, the distribution shape is of a quadrilateral trapezoid.
	Altering \samp{n1}, \samp{n3}, or \samp{alpha} changes the growth rate, decay rate, and
	boundary ratio parameters, respectively, as demonstrated in the examples below.
}
\seealso{
	\link{Distributions} for standard distributions.
}
\references{
	van Dorp, J. R. and Kotz, S. (2003)
	Generalized trapezoidal distributions.
	\emph{Metrika}. \bold{58}(1):85--97.
	Preprint available: \url{http://www.seas.gwu.edu/~dorpjr/Publications/JournalPapers/Metrika2003VanDorp.pdf}
	
	van Dorp, J. R., Rambaud, S.C., Perez, J. G., and Pleguezuelo, R. H. (2007)
	An elicitation procedure for the generalized trapezoidal distribution with a uniform central stage.
	\emph{Decision Analysis Journal}. \bold{4}:156--166.
	Preprint available: \url{http://www.seas.gwu.edu/~dorpjr/Publications/JournalPapers/DA2007.pdf}
}
\author{Jeremy Thoms Hetzel \email{jthetzel@gmail.com}}
\examples{
	## Plot default trapezoid distribution
	curve(dtrapezoid(x, min=0, mode1=1/3, mode2=2/3, max=1, n1=2, n3=2, alpha=1),
		from=0, to=1)
	
	
	## Plot triangular trapezoid distribution
	curve(dtrapezoid(x, min=0, mode1=1/2, mode2=1/2, max=1, n1=2, n3=2, alpha=1),
		from=0, to=1)
	
	
	## Explore effects of n1, n3, and alpha parameters
	# plyr and ggplot2 are required for this example
	require(plyr)
	require(ggplot2)
	
	x <- seq(from = 0, to = 1, by = 0.01)	
	
# Create a lists of arguments, varying n1, n3, and alpha
	arguments <- list()
	arguments[['A']] <- list(x = x, n1 = 2, n3 = 2, alpha = 0.8)
	arguments[['B']] <- list(x = x, n1 = 1.5, n3 = 1.5, alpha = 1)
	arguments[['C']] <- list(x = x, n1 = 2.5, n3 = 2.5, alpha = 1.5)
	arguments[['D']] <- list(x = x, n1 = 1.5, n3 = 2.5, alpha = 0.5)
	arguments[['E']] <- list(x = x, n1 = 2.5, n3 = 1.5, alpha = 1)
	arguments[['F']] <- list(x = x, n1 = 0.5, n3 = 0.5, alpha = 1.5)
	arguments[['G']] <- list(x = x, n1 = 1.5, n3 = 0.5, alpha = 0.5)
	arguments[['H']] <- list(x = x, n1 = 2.5, n3 = 0.5, alpha = 1)
	arguments[['I']] <- list(x = x, n1 = 0.5, n3 = 1.5, alpha = 1.5)
	arguments[['J']] <- list(x = x, n1 = 0.5, n3 = 2.5, alpha = 0.5)
	
# Calculate the distributions
	plot.data <- ldply(arguments, function(z)
		{
			x <- z$x
			density <- dtrapezoid(x = z$x, min = 0, mode1 = 0.2, mode2 = 0.8,
				max = 1, n1 = z$n1, n3 = z$n3, alpha = z$alpha)
			args <- paste("n1 = ", z$n1, ", n3 = ", z$n3, ", alpha = ", z$alpha, 
				sep="", collapse="")
			out <- data.frame(x, density, args)
		})
	
# Create labels for later use in displaying the arguments on the plots 
	plot.data$label <- paste(plot.data$.id, ": ", plot.data$args, sep="")
	
# Create plots
	generalizedTrapezoids <- ggplot(data = plot.data, aes(x = x, y = density)) + 
		geom_line() + theme_bw() + 
		facet_wrap(~label, ncol = 2, scales = "free_y")
	print(generalizedTrapezoids)
}
\keyword{distribution}
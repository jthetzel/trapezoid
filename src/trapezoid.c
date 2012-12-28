#include <R.h>
#include <Rinternals.h>

/* Function to calculate quantile given an vector of probabilities*/
SEXP qtrapezoid(SEXP p_, SEXP a_, SEXP b_, SEXP c_, SEXP d_, SEXP m_, SEXP n_, SEXP alpha_)
{
	double *p_array = REAL(p_);
	double a = asReal(a_);
	double b = asReal(b_);
	double c = asReal(c_);
	double d = asReal(d_);
	double m = asReal(m_);
	double n = asReal(n_);
	double alpha = asReal(alpha_);
	int len = length(p_);

	SEXP result;
	PROTECT(result = allocVector(REALSXP, len));
	double *q_array = REAL(result);

	/* Calculate pi1, pi2, and pi3 */
	double pi1 = (2 * alpha * (b - a) * n) /((2 * alpha * (b - a) * n) + ((alpha + 1) * (c - b) * m * n) +	(2 * (d - c) * m));
	double pi2 = ((alpha + 1) * (c - b) * m * n) /	((2 * alpha * (b - a) * n) + ((alpha + 1) * (c - b) * m * n) + (2 * (d - c) * m));
	double pi3 = (2 * (d - c) * m) / ((2 * alpha * (b - a) * n) + ((alpha + 1) * (c - b) * m * n) + (2 * (d - c) * m));

	/* Calculate quantile conditional on p */
	for (int i = 0; i < len; i++)
	{
		double p = p_array[i];

		if ((0 <= p) && (p <= pi1))
		{
			q_array[i] = pow(((p*(2*alpha*(b-a)*n+(alpha+1)*(c-b)*m*n+2*(d-c)*m))/(2*alpha*(b-a)*n)),(1/m))*(b-a)+a;
		}
		else if ((pi1 < p) && (p <= (1 - pi3)) && (alpha != 1))
		{
			q_array[i] = ((sqrt(pow((((-2)*b*m*n*(1-alpha))/(2*(c-b))+2*m*n*(((2*c-b)*(alpha-1))/(2*(c-b))+1)),2)/pow((2*alpha*(b-a)*n+(alpha+1)*(c-b)*m*n+2*(d-c)*m),2)-(((2*alpha*(b-a)*n+(-2)*b*m*n*(((2*c-b)*(alpha-1))/(2*(c-b))+1))/(2*alpha*(b-a)*n+(alpha+1)*(c-b)*m*n+2*(d-c)*m)-p)*4*2*m*n*(1-alpha))/(2*(c-b)*(2*alpha*(b-a)*n+(alpha+1)*(c-b)*m*n+2*(d-c)*m)))-(((-2)*b*m*n*(1-alpha))/(2*(c-b))+2*m*n*(((2*c-b)*(alpha-1))/(2*(c-b))+1))/(2*alpha*(b-a)*n+(alpha+1)*(c-b)*m*n+2*(d-c)*m))*2*(c-b)*(2*alpha*(b-a)*n+(alpha+1)*(c-b)*m*n+2*(d-c)*m))/(2*2*m*n*(1-alpha));
		}
		else if ((pi1 < p) && (p <= (1 - pi3)) && (alpha == 1))
		{
			q_array[i] = b + (((p - pi1) / (pi2)) * (c - b));
		}
		else if (((1 - pi3) < p) && (p <= 1))
		{
			q_array[i] = d-pow((((1-p)*(2*alpha*(b-a)*n+(alpha+1)*(c-b)*m*n+2*(d-c)*m))/(2*(d-c)*m)),(1/n))*(d-c);
		}
		else
			q_array[i] = NAN;
	}

	UNPROTECT(1);
	return(result);
}

/* Function to calculate probability given an vector of quantiles*/
SEXP ptrapezoid(SEXP q_, SEXP a_, SEXP b_, SEXP c_, SEXP d_, SEXP m_, SEXP n_, SEXP alpha_)
{
	double *q_array = REAL(q_);
	double a = asReal(a_);
	double b = asReal(b_);
	double c = asReal(c_);
	double d = asReal(d_);
	double m = asReal(m_);
	double n = asReal(n_);
	double alpha = asReal(alpha_);
	int len = length(q_);

	SEXP result;
	PROTECT(result = allocVector(REALSXP, len));
	double *p_array = REAL(result);

	/* Calculate probability conditional on q */
	for (int i = 0; i < len; i++)
	{
		double q = q_array[i];

		if (q < a)
		{
			p_array[i] = 0;
		}
		else if ((a <= q) && (q < b))
		{
			p_array[i] = (2 * alpha * (b - a) * n) / ((2 * alpha * (b - a) * n) + ((alpha + 1) * (c - b) * m * n) + (2 * (d - c) * m)) * pow((q - a) / (b - a), m);
		}
		else if ((b <= q) && (q < c))
		{
			p_array[i] = ((2 * alpha * (b - a) * n) + (2 * (q - b) * m * n * (1 + ((alpha - 1) / 2) * ((2*c - b - q) / (c - b))))) / ((2 * alpha * (b - a) * n) + ((alpha + 1) * (c - b) * m * n) + (2 * (d - c) * m));
		}
		else if ((c <= q) && (q < d))
		{
			p_array[i] = 1 - ((2 * (d - c) * m) / ((2 * alpha * (b - a) * n) + ((alpha + 1) * (c - b) * m * n) + (2 * (d - c) * m))) * pow((d - q) / (d - c), n);
		}
		else if (q >= d)
		{
			p_array[i] = 1;
		}
		else
		{
			p_array[i] = NAN;
		}
	}

	UNPROTECT(1);
	return(result);
}

/* Function to calculate cumulative density given an vector of quantiles*/
SEXP dtrapezoid(SEXP x_, SEXP a_, SEXP b_, SEXP c_, SEXP d_, SEXP m_, SEXP n_, SEXP alpha_)
{
	double *x_array = REAL(x_);
	double a = asReal(a_);
	double b = asReal(b_);
	double c = asReal(c_);
	double d = asReal(d_);
	double m = asReal(m_);
	double n = asReal(n_);
	double alpha = asReal(alpha_);
	int len = length(x_);

	SEXP result;
	PROTECT(result = allocVector(REALSXP, len));
	double *d_array = REAL(result);

	/* Calculate normalizing_constant */
	double normalizing_constant = (2 * m * n) / ((2 * alpha * (b - a) * n) + ((alpha + 1) * (c - b) * m * n) + (2 * (d - c) * m));

	/* Calculate cumulative density conditional on x */
	for (int i = 0; i < len; i++)
	{
		double x = x_array[i];

		if ((a <= x) && (x < b))
		{
			d_array[i] = normalizing_constant * alpha * pow((x - a) / (b - a),(m - 1));
		}
		else if ((b <= x) && (x < c))
		{
			d_array[i] = normalizing_constant * (((1 - alpha) * ((x - b) / (c - b))) + alpha);
		}
		else if ((c <= x) && (x <= d))
		{
			d_array[i] = normalizing_constant * pow((d - x) / (d - c), (n - 1));
		}
		else
		{
			d_array[i] = 0;
		}
	}

	UNPROTECT(1);
	return(result);
}

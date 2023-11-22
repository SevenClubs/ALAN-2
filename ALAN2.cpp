#include <iostream>
#include <math.h>

double exp_Taylor (double x, int N) 
{
	double ans = 0;
	for (int i = 0; i <= N; ++i) 
	{
		double num = pow(x, i);
		for (int j = 2; j <= i; ++j)
			num /= j;
		ans += num;
	}
	return ans;
}

double exp_Taylor2(double x, int N)
{
	if (x > 0)
		return exp_Taylor(x, N);
	else
		return 1 / exp_Taylor(-x, N);
}

int main()
{
	std::cout << "\nEhi, mi daresti la x e la N, bellezza? ;)\n";
	double x;
	int N;
	std::cin >> x >> N;
	
	double r0 = exp(x);
	double r1 = exp_Taylor(x, N);
	double r2 = exp_Taylor2(x, N);
	
	std::cout << "\nRisultati esperimento:\nr0 = " << r0 << "\nr1 = " << r1 << "\nr2 = " << r2;
	
	double abs_err1 = abs(r0 - r1);
	double rel_err1 = abs_err1 / r0;
	
	double abs_err2 = abs(r0 - r2);
	double rel_err2 = abs_err2 / r0;
	
	std::cout << "\n\n - Errori sul primo algoritmo:\n\nErrore assoluto = " << abs_err1 << "\nErrore relativo = " << rel_err1 * 100 << "%";
	std::cout << "\n\n - Errori sul secondo algoritmo:\n\nErrore assoluto = " << abs_err2 << "\nErrore relativo = " << rel_err2 * 100 << "%\n\n---------------------------\n\n";

}

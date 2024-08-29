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
	double x_val[2] = {0.5, 30};
	int N_val[5] = {3, 10, 50, 100, 150}; 
	for (int i = 0; i < 2; ++i)
	{
		double x = x_val[i];
		for (int j = 0; j < 5; ++j)
		{	
			int N = N_val[j];
			
			double expected_res = exp(x);
			double res = exp_Taylor(x, N);
			
			std::cout << "\nRisultati esperimento(x = " << x << ", N = " << N << "):\n\nValore atteso = " << expected_res << "\nValore calcolato = " << res;
			
			double abs_err = abs(expected_res - res);
			double rel_err = abs_err / expected_res;
			
			std::cout << "\n\n - Errori sull'algoritmo:\n\nErrore assoluto = " << abs_err << "\nErrore relativo = " << rel_err << "\n\n---------------------------\n\n";
		}
	}
	for (int i = 0; i < 2; ++i)
	{
		double x = -x_val[i];
		for (int j = 0; j < 5; ++j)
		{	
			int N = N_val[j];
			
			double expected_res = exp(x);
			double res = exp_Taylor(x, N);
			
			std::cout << "\nRisultati esperimento(x = " << x << ", N = " << N << ", algoritmo 1):\n\nValore atteso = " << expected_res 
				<< "\nValore calcolato = " << res;
			
			double abs_err = abs(expected_res - res);
			double rel_err = abs_err / expected_res;
			
			std::cout << "\n\n - Errori sull'algoritmo 1:\n\nErrore assoluto = " << abs_err << "\nErrore relativo = " << rel_err << "\n\n---------------------------\n\n";
		}
	}

	for (int i = 0; i < 2; ++i)
	{
		double x = -x_val[i];
		for (int j = 0; j < 5; ++j)
		{	
			int N = N_val[j];
			
			double expected_res = exp(x);
			double res = exp_Taylor2(x, N);
			
			std::cout << "\nRisultati esperimento(x = " << x << ", N = " << N << ", algoritmo 2):\n\nValore atteso = " << expected_res 
				<< "\nValore calcolato = " << res;
			
			double abs_err = abs(expected_res - res);
			double rel_err = abs_err / expected_res;
			
			std::cout << "\n\n - Errori sull'algoritmo 2:\n\nErrore assoluto = " << abs_err << "\nErrore relativo = " << rel_err << "\n\n---------------------------\n\n";
		}
	}
}

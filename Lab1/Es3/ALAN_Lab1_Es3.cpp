#include <iostream>
#include <math.h>

int main()
{
	long int d = 0;
	long int f = 0;
	
	{
		float one = 1;
		float power_two = pow(2, 0);
		while (one + power_two > one)
		{
			++f;
			power_two = pow(2, -f);
		}
		--f;
		std::cout << "\nf = " << f << "\neps = " << pow(2, -f) << "\n\n";
	}
	
	{
		double one = 1;
		double power_two = pow(2, 0);
		while (one + power_two > one)
		{
			++d;
			power_two = pow(2, -d);
		}
		--d;
		std::cout << "d = " << d << "\neps = " << pow(2, -d) << "\n\n";
	}
}

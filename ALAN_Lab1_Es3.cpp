#include <iostream>
#include <math.h>

int main()
{
	long int d = 0;
	long int o = 0;
	
	{
		double one = 1;
		double lollo = pow(2, 0);
		while (one + lollo > one)
		{
			++o;
			std::cout << lollo << "\n";
			lollo = pow(2, -o);
			if (o > 100)
				return 0;
		}
	}
	
	{
		float one = 1;
		float lollo = pow(2, 0);
		while (one + lollo > one)
		{
			++d;
			std::cout << d << "\n";
			lollo = pow(2, -d);
		}
	}
	
	--d;
	--o;
	
	std::cout << "d = " << d << ", o = " << o << "\n";
}

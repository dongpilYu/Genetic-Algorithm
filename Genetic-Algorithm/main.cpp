#include "ga.h"
//#include "simulation.h"

int main()
{
	double avg_best_fiteness[GENERATIONS] = { 0, };
	double std_best_fiteness[GENERATIONS] = { 0, };
	double var_best_fiteness[GENERATIONS] = { 0, };

	for (int i = 0; i < ITERATION; i++)
	{
		srand(time(0));
		set_input_parameters();
		Population pop;
		pop.run();
		pop.statistics_info_calc();
		cout << "<------------" << i << "th iteration------------>" << endl;
		//pop.print_population(MODE_DETAIL);
	}
	cout << "<------------Result------------>" << endl;
	

	/*
	for (int i = 0; i < GENERATIONS ; i++)
	{
		avg_best_fiteness[i] = sum_best_fiteness[i] / ITERATION;
		if (i % 10000 == 0)
		{
			cout << avg_best_fiteness[i] << " " << sum_avg_fiteness[i] / ITERATION << endl;
		}
	}
	*/
/*
	for (int i = 0; i < GENERATIONS; i++)
	{
		avg_best_fiteness[i] = sum_best_fiteness[i] / ITERATION;
		for (int j = 0; j < ITERATION; j++)
		{
			var_best_fiteness[i] += pow(gene_best_fitness[j][i] - avg_best_fiteness[i], 2);
		}
		std_best_fiteness[i] = sqrt(var_best_fiteness[i]);

		if ((i == 0) || (i % 10000 == 0) || (i == GENERATIONS - 1))
		{
			cout << avg_best_fiteness[i] << " " << std_best_fiteness[i] << endl;
		}
	}

	int temp;
	cin >> temp;
	*/

	//Random
	cout.setf(ios::fixed);
	cout.precision(6);
	int random_avg_best_fitness = sum_best_fiteness[0] / ITERATION;
	/*
	int random_var_best_fitness = 0;
	int random_std_best_fitness = 0;

	for (int i = 0; i < ITERATION; i++)
	{
		random_var_best_fitness += pow(random_gene_best_fitness[i] - random_avg_best_fitness, 2);
	}
	random_std_best_fitness = sqrt(random_var_best_fitness);

	cout << random_avg_best_fitness << " " << random_std_best_fitness << endl;
*/
	for (int i = 0; i < ITERATION; i++)
	{
		cout << random_gene_best_fitness[i] << endl;
	}
	int temp;
	cin >> temp;
}

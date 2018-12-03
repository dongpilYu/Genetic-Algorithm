#include "ga.h"
#include "simulation.h"

int main()
{
	/*요기*/
	clock_t begin, end;
	begin = clock();

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
		/*cout << "assignment: " << endl;
		
		for (int j = 0; j < K; j++)
		{
			cout << assignment[j] << endl;
		}
		cout << endl;*/
		//pop.print_population(MODE_DETAIL);
		pop.print_best_assignment();
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
	/*요기*/

	cout << "----------적합도-----------" << endl;
	for (int i = 0; i < GENERATIONS; i++)
	{
		avg_best_fiteness[i] = sum_best_fiteness[i] / ITERATION;
		for (int j = 0; j < ITERATION; j++)
		{
			var_best_fiteness[i] += pow(gene_best_fitness[j][i] - avg_best_fiteness[i], 2);
		}
		std_best_fiteness[i] = sqrt(var_best_fiteness[i]);

		if ((i == 0) || (i % 100 == 0) || (i == GENERATIONS - 1))
		{
			cout << avg_best_fiteness[i] << " " << std_best_fiteness[i] << endl;
		}
	}

	//Random
	/*
	cout.setf(ios::fixed);
	cout.precision(6);
	int random_avg_best_fitness = sum_best_fiteness[0] / ITERATION;
	
	int random_var_best_fitness = 0;
	int random_std_best_fitness = 0;

	for (int i = 0; i < ITERATION; i++)
	{
		random_var_best_fitness += pow(random_gene_best_fitness[i] - random_avg_best_fitness, 2);
	}
	random_std_best_fitness = sqrt(random_var_best_fitness);

	cout << random_avg_best_fitness << " " << random_std_best_fitness << endl;

	for (int i = 0; i < ITERATION; i++)
	{
		cout << random_gene_best_fitness[i] << endl;
	}//요기 밑에 풀기*/
	end = clock();
	cout << "수행시간 : " << (end - begin) << endl;
	
	
}



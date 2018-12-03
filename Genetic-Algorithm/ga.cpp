#include "ga.h"
#include "simulation.h"

int W[K] = { 0, }, C[K] = { 0, }, R[M] = { 0, };
double P[M] = { 0, }, Q[M] = { 0, };
double sum_best_fiteness[GENERATIONS] = { 0, };
double sum_avg_fiteness[GENERATIONS] = { 0, };
double gene_best_fitness[ITERATION][GENERATIONS] = { 0, };
double random_gene_best_fitness[ITERATION] = { 0, };
int assignment[K] = { 0, };
int iteration_num = 0;

int SCENARIO[S] = {
	8500,  1200, 45000, 3800,
	8500,   600, 45000,  800,
	500,   1700,  1200, 2500,
	45000,  800,   500,   50
};

/*
Post: Randomly generates data within range
*/
Chromosome::Chromosome()
{
	int rand_binary;

	for (int row = 0; row < M; row++) {
		for (int col = 0; col < K; col++)
		{
			//data[row][col] = (rand() % ((int)(2 * (double)R[row] / K) - LB + 1)) + LB;
			//data[row][col] = C[col] + rand() % (int)(C[col] * RAND_PERCENTAGE);
			rand_binary = rand() % 2;
			if (rand_binary == 0)
			{
				data[row][col] = C[col] + rand() % (int)(C[col] * RAND_PERCENTAGE);
			}
			else
			{
				data[row][col] = C[col] - rand() % (int)(C[col] * RAND_PERCENTAGE);
			}
		}
	}
	cumulative_sum_calc();
	
	repair();
	//print_chromo();
	//cumulative_sum_calc();
	evaluate();
	//cout << "Chromosome Constructor" << endl;
}
/*
Pre: Arrays that store integer types are declared
Post: Calculates the cumulative sum of rows and columns and stores them in the array row_sum and column_sum
*/
void Chromosome::cumulative_sum_calc()
{
	for (int row = 0; row < M; row++) {
		int sum = 0;
		for (int col = 0; col < K; col++)
		{
			sum += data[row][col];
		}
		row_sum[row] = sum;
	}

	for (int col = 0; col < K; col++) {
		int sum = 0;
		for (int row = 0; row < M; row++)
		{
			sum += data[row][col];
		}
		column_sum[col] = sum;
	}
	/*
	cout << "Column Sum" << endl;
	for (int col = 0; col < K; col++)
	{
	cout << column_sum[col] << " ";
	}
	cout << "\n";

	cout << "Row Sum" << endl;
	for (int row = 0; row < M; row++)
	{
	cout << row_sum[row] << " ";
	}
	cout << "\n";*/
}

/*
Post: Prints a chromosome
*/
void Chromosome::print_chromo()
{
	for (int row = 0; row < M; row++) {
		for (int col = 0; col < K; col++)
		{
			//cout.width(K);
			cout << setw(3) << data[row][col] << " "; 
		}
		cout << "\n";
	}
	cout << "\n";
}

/*
Post: Generates a random real number from 0.0 to 1.0 for each data,
and randomly changes the data if the generated real number is less than the defined mutation rate
*/
void Chromosome::mutate()
{
	// one point flip
	for (int row = 0; row < M; row++) {
		for (int col = 0; col < K; col++)
		{
			if ((double)rand() / RAND_MAX < MUTATION_RATE)
			{
				data[row][col] = rand() % R[row] - LB + 1;
			}
		}
	}
	//cout << "Mutation" << endl;
	return;
}

/*
Pre: Chromosome created by crossing and mutation is not suitable for constraints
Post: The chromosome is modified to fit the constraints entered from the file
*/
void Chromosome::repair()
{
	//cout << "Repair start" << endl;
	cumulative_sum_calc();

	/* Constraint R */
	for (int row = 0; row < M; row++)
	{
		while (row_sum[row] != R[row])
		{
			if (row_sum[row] < R[row])
			{
				data[row][rand() % K]++;
			}
			else if (row_sum[row] > R[row])
			{
				int rand_row = rand() % K;

				if (data[row][rand_row] > MIN_ASN)
				{
					data[row][rand_row]--;
				}
			}
			cumulative_sum_calc();				//for debug
		}
	}
	assert_feasible(MODE_R_ERROR);

	//print_chromo();
	//cout << "Constraint R is done" << endl;

	/*max constraint*/
	for (int col = 0; col < K; col++)
	{
		if (data[0][col] > MAX_ASN)
		{
			int excess_value = 0;
			excess_value = (data[0][col] - MAX_ASN) / (K - 1);
			
			for (int i = 0; i < col; i++)
			{
				data[0][i] += excess_value;
			}
			data[0][col] = MAX_ASN;
			for (int i = col+1; i < K; i++)
			{
				data[0][i] += excess_value;
			}
		}
	}

	/* Constraint C */
	/*
	int flag_arr[K] = { 0, };
	int count = 0;
	while (1)
	{
		count++;
		if (count > 1000000)  //임시코드
		{
			cout << "Repair count over" << endl;
			print_chromo();
			cout << "Column Sum" << endl;
			for (int col = 0; col < K; col++)
			{
				cout << column_sum[col] << " ";
			}
			cout << "\n";
			cout << "Row Sum" << endl;
			for (int row = 0; row < M; row++)
			{
				cout << row_sum[row] << " ";
			}
			cout << "\n";
			break;
		}
		for (int col = 0; col < K; col++) //flag check, 원래는 MAX_ASN이 C[col]
		{
			if (column_sum[col] > MAX_ASN) flag_arr[col] = FLAG_EXCESS;
			else if (column_sum[col] < MAX_ASN) flag_arr[col] = FLAG_SHORTAGE;
			else flag_arr[col] = FLAG_SAME;
		}

		//Escape examination
		bool flag_escape = ESCAPABLE;
		for (int col = 0; col < K; col++)
		{
			if (flag_arr[col] == FLAG_EXCESS)
			{
				flag_escape = NON_ESCAPABLE;
			}
		}

		if (flag_escape == ESCAPABLE)
		{
			assert_feasible(MODE_C_ERROR);
			assert_feasible(MODE_R_ERROR);
			break;
		}

		int rand_row = rand() % M;
		int rand_col = 0;
		int excess_col, shortage_col;
		int old_excess_data = 0, old_shortage_data = 0;

		while (flag_arr[rand_col] != FLAG_EXCESS)
		{
			rand_col = rand() % K;
		}
		excess_col = rand_col;

		if (count > 999900) //디버깅용 임시코드
			cout << "excess_col is " << data[rand_row][excess_col] << ", " << rand_col << "th column" << endl;

		while (flag_arr[rand_col] != FLAG_SHORTAGE)
		{
			rand_col = rand() % K;
		}
		shortage_col = rand_col;

		if (count > 999900) //디버깅용 임시코드
			cout << "shortage_col is " << data[rand_row][shortage_col] << ", " << rand_col << "th column" << endl;

		if (data[rand_row][shortage_col] < data[rand_row][excess_col])
		{
			old_excess_data = data[rand_row][excess_col];
			old_shortage_data = data[rand_row][shortage_col];

			if ((data[rand_row][shortage_col] < KEEP + LB) && (data[rand_row][excess_col] >= KEEP + LB))
			{
				//cout << "old shortage col is " << data[rand_row][shortage_col] << " and excess col is " << data[rand_row][excess_col] << endl;

				int temp = 0;
				temp = data[rand_row][shortage_col] + KEEP;
				data[rand_row][shortage_col] = data[rand_row][excess_col] - KEEP;
				data[rand_row][excess_col] = temp;
				//cout << "swapping with keep: shortage col is "<< data[rand_row][shortage_col] << " and excess col is " << data[rand_row][excess_col] << endl;
			}
			else
			{
				int temp = 0;
				temp = data[rand_row][shortage_col];
				data[rand_row][shortage_col] = data[rand_row][excess_col];
				data[rand_row][excess_col] = temp;
			}

			//If the if statement does not swap properly
			if ((old_excess_data == data[rand_row][excess_col]) && (old_shortage_data == data[rand_row][shortage_col]))
			{
				int temp = 0;
				temp = data[rand_row][shortage_col];
				data[rand_row][shortage_col] = data[rand_row][excess_col];
				data[rand_row][excess_col] = temp;
			}
		}
		cumulative_sum_calc();
	}*/
	//cout << "Repair done" << endl;
}

/*
Pre: column_sum, row_sum must contain value
Post: Decide whether to run Repair again
*/
void Chromosome::assert_feasible(int mode)
{
	for (int row = 0; row < M; row++) {
		for (int col = 0; col < K; col++)
		{
			if (data[row][col] < MIN_ASN) //원래는 LB였음
			{
				//cout << "Repaired less than LB" << endl;
				data[row][col]++;
				//print_chromo();
				repair();
			}
		}
	}
	if (mode == MODE_C_ERROR) //col, MAX_ASN이 원래는 C[col]
	{
		for (int col = 0; col < K; col++)
		{
			if (column_sum[col] > MAX_ASN)
			{
				cout << "Repair Error in constraint C" << endl;
				//repair();
			}

		}
	}
	else if (mode == MODE_R_ERROR) //row
	{
		for (int row = 0; row < M; row++)
		{
			if (row_sum[row] != R[row])
			{
				cout << "Repair Error in constraint R" << endl;
				repair();
			}
		}
	}
}


/*
Post: Evaluate the created GA
*/
void Chromosome::evaluate() //It's a temporary code, and we'll have to fix it later
{
	double eval_sum = 0;
	double eval_basic_sum = 0;
	//double damage[K] = { 0, };

	cumulative_sum_calc();
	/*
	for (int col = 0; col < K; col++)//기존 Evaluation
	{
	eval_sum += column_sum[col] * W[col];
	}

	for (int col = 0; col < K; col++) {
		double prod = 1;
		for (int row = 0; row < M; row++)
		{
			prod *= pow(Q[row], data[row][col]); //row마다 파괴확률을 곱
		}
		damage[col] = prod;
		assert(prod >= 0 && prod < 1, "Wrong range");
	}
	for (int col = 0; col < K; col++)
	{
		eval_sum += damage[col] * W[col]; // damage 마다 가중치를 곱
	}

	fitness = eval_sum;
	*/
	
	for (int i = 0; i < K; i++)
	{
		assignment[i] = data[0][i];
	}
	
	int weight = 1;
	double time[S];

	sim::simulation(data, time);
	
	
	for (int i = 0; i<K; ++i) {
		//std::cout <<"time: " << time[i] << ',';
		double prod = 1;
		double basic_prod = 1;
		if (time[i] < 2) //2일 내에 끝나는 경우
		{
			weight = PENALTY_1DAY;
		}
		else if (time[i] >= 2 && time[i] < 3) //2~3일
		{
			weight = PENALTY_2DAY;
		}
		else // 3일 넘게 걸리는 경우
		{
			weight = PENALTY_EXCESS;
		}
		prod = SCENARIO[i] * time[i] * weight;
		basic_prod = SCENARIO[i] * time[i];
		eval_sum += prod;
		eval_basic_sum += basic_prod;
		//cout << i << " : " << time[i] << endl;
	}
	fitness = eval_sum / 165650;
	
	//std::cout << std::endl; 
	

	/*
	for (int i = 0; i<M; ++i) {
		data[i][0] = cur_mar_skimmer[i][0] + cur_gro_skimmer[i][0];
		std::cout << data[i][0] << ',';
	}
	std::cout << std::endl;
	*/
	
}

void Chromosome::simple_simulation(const int data[M][K], double(&time)[S])
{
	sim::simple_simulation(data, time);
}

/////////////////////////////////////////////////////////////////////////////////////
/*
Pre: Population declared
Post: Individuals are created for a defined population size. Compute the fitness statistically and put the value into the variable
*/
Population::Population(int pop_size)
{
	individual = new Chromosome[POP_SIZE];


	statistics_info_calc();
	cout << "Population Constructor" << endl;
}

/*
Pre: The main function declares a population object and calls run()
Post: Evolutionary operations up to the number of generations
*/
void Population::run()
{
	// steady state GA
	cout << "run() Start" << endl;
	for (int i = 0; i < GENERATIONS; i++)
	{
		int parent1 = select();
		int parent2 = select();
		while (parent1 == parent2)
		{
			parent2 = select();
			//cout << parent1 << "," << parent2 << endl;
		}
		crossover(parent1, parent2);
		assert(parent1 != parent2, "Wrong selection");
		offspring.mutate();
		//offspring.print_chromo();
		offspring.repair();
		//offspring.print_chromo();
		offspring.evaluate();
		replacement(parent1, parent2);
		statistics_info_calc();

		sum_best_fiteness[i] += get_best_fitness();
		sum_avg_fiteness[i] += get_avg_fitness();
		gene_best_fitness[iteration_num][i] = get_best_fitness();
		//random
		random_gene_best_fitness[iteration_num] = get_best_fitness();
	}
	iteration_num++;
	//offspring.print_chromo();
	cout << "run() End" << endl;
}

/*
Pre: Enter one of the defined modes
Post: Print informations according to mode
*/
void Population::print_population(int mode)
{
	Chromosome individual;
	if (mode == MODE_BRIEF)
	{
		cout << "********Brief Mode********" << endl;
		cout << "(1) Best fitness: " << best_fitness << endl;
		cout << "(2) Average fitness: " << avg_fitness << endl;
	}
	else if (mode == MODE_DETAIL)
	{
		cout << "********Detail Mode********" << endl;
		cout << "(1) Best fitness: " << best_fitness << endl;
		cout << "(2) Average fitness: " << avg_fitness << endl;
		cout << "(3) Worst fitness: " << worst_fitness << endl;
		cout << "(4) Best index" << best_index << endl;
		cout << "(5) Worst index" << worst_index << endl;

		for (int i = 0; i < POP_SIZE; i++)
		{
			cout << "(6) Chromosome: " << endl;
			cout << individual.get_fitness() << endl;
		}
	}
	else
	{
		cout << "********Test Mode********" << endl;
		cout << "(1) Best fitness: " << best_fitness << endl;
		cout << "(2) Average fitness: " << avg_fitness << endl;
		cout << "(3) Worst fitness: " << worst_fitness << endl;
		cout << "(4) Best index" << best_index << endl;
		cout << "(5) Worst index" << worst_index << endl;
	}
}

/*
Pre:
sum of fitness must be greater than 0
Post: Enter a value of Parents
*/
int Population::select()
{
	double point;
	double sum_of_fitness, sum;
	int idx;

	idx = 0;
	sum_of_fitness = 0;
	sum = 0;

	for (int i = 0; i < POP_SIZE; i++)
	{
		sum_of_fitness += individual[i].get_fitness();
		//cout << "in for"<< idx << "," << sum_of_fitness << "," << sum << endl;
	}

	point = sum_of_fitness * (double)rand() / RAND_MAX;
	//cout << "after make point" << idx << "," << sum_of_fitness << "," << sum << endl;
	//cout << "point : " << point << endl;

	for (int i = 0; i < POP_SIZE; i++)
	{
		sum += individual[i].get_fitness();
		if (point <= sum)
		{
			idx = i;
			break;
		}
	}
	//cout << "Selection : "<< idx << endl;
	return idx;
}

/*
Pre: Two different parents must be selected
Post: Two chromosomes are mixed based on one point
*/
void Population::crossover(int parent1, int parent2)
{/*
	//block uniform 
	int m, k;

	m = rand() % M + 1;
	k = rand() % K + 1;

	for (int row = 0; row < M; row++) {
		for (int col = 0; col < K; col++)
		{
			if (row < m && col < k)
			{
				offspring.put_data(row, col, individual[parent1].get_data(row, col));

			}
			else if (row < m && col >= k)
			{
				offspring.put_data(row, col, individual[parent2].get_data(row, col));
			}
			else if (row >= m && col < k)
			{
				offspring.put_data(row, col, individual[parent1].get_data(row, col));
			}
			else
			{
				offspring.put_data(row, col, individual[parent2].get_data(row, col));
			}
		}
	}

*/
//uniform
	int rand_binary;

	for (int row = 0; row < M; row++) {
		for (int col = 0; col < K; col++)
		{
			rand_binary = rand() % 2;
			if (rand_binary == 0)
			{
				offspring.put_data(row, col, individual[parent1].get_data(row, col));
			}
			else
			{
				offspring.put_data(row, col, individual[parent2].get_data(row, col));
			}
		}
	}

	//cout << "Crossover" << endl;
	//offspring.print_chromo();
}

/*
Pre: Enter two parents
Post: Replace parent offspring with parent
*/
void Population::replacement(int parent1, int parent2)
{
	double distance1 = 0, distance2 = 0;

	distance1 = (double)fabs(individual[parent1].get_fitness() - offspring.get_fitness());
	distance2 = (double)fabs(individual[parent2].get_fitness() - offspring.get_fitness());
	/*
	//  < Maximizing code >
	if (distance1 <= distance2) {
	if( (parent1 != best_index) && (individual[parent1].get_fitness() < offspring.get_fitness()))	//elitism
	individual[parent1] = offspring;
	}
	else {
	if ((parent2 != best_index) && (individual[parent2].get_fitness() < offspring.get_fitness()))	//elitism
	individual[parent2] = offspring;
	}
	*/

	//  < Minimizing code >
	if (distance1 <= distance2) {
		if ((parent1 != best_index) && (individual[parent1].get_fitness() > offspring.get_fitness()))	//elitism
			individual[parent1] = offspring;
	}
	else {
		if ((parent2 != best_index) && (individual[parent2].get_fitness() > offspring.get_fitness()))	//elitism
			individual[parent2] = offspring;
	}

	//cout << "Replacement" << endl;
}

void Population::print_best_assignment()
{
	int assignment[K] = { 0, };
	cout << "best assinment" << endl;
		for (int row = 0; row < M; row++) {
			for (int col = 0; col < K; col++)
			{
				//cout.width(K);
				cout << setw(3) << individual[best_index].get_data(row, col) << endl;
				assignment[col] = individual[best_index].get_data(row, col);
			}
			//cout << "\n";
		}
	//cout << "\n";
		assignment_eval(assignment);
}

/*
Pre: Enter data used as a constraint in file format
Post: Parameters used as constraints are set
*/
void set_input_parameters()
{
	vector<int> numbers;
	ifstream inputFile1("input_test.txt");

	if (inputFile1.good()) {
		int current_number = 0;
		while (inputFile1 >> current_number) {
			numbers.push_back(current_number);
		}
		inputFile1.close();

		for (int i = 0; i < K; i++)
		{
			W[i] = numbers[i];
			//cout << W[i] << endl;
		}
		for (int i = K; i < K + K; i++)
		{
			C[i - K] = numbers[i];
			//cout << C[i - K] << endl;
		}
		for (int i = K + K; i < K + K + M; i++)
		{
			R[i - (K + K)] = numbers[i];
			//cout << R[i - (K + K)] << endl;
		}
	}

	vector<double> probabilities;
	ifstream inputFile2("damage_input.txt");

	if (inputFile2.good()) {
		double current_number = 0.0;
		while (inputFile2 >> current_number) {
			probabilities.push_back(current_number);
		}
		inputFile2.close();

		for (int i = 0; i < M; i++)
		{
			P[i] = probabilities[i];
			Q[i] = 1 - P[i];
			//cout << P[i]<<", " <<Q[i] << endl;
		}
	}
	cout << "Parameter entered" << endl;
}

/*
Pre: Member variable must be declared
Post: Variables that store statistical information are calculated
*/
void Population::statistics_info_calc()
{
	//cout << "statistics_info_calc" << endl;
	int max_idx = 0, min_idx = 0;
	double max = 0;
	double min = 1000;
	double sum_of_fitness = 0;

	for (int i = 0; i < POP_SIZE; i++)
	{
		//cout << i << "th fitness" << individual[i].get_fitness() << endl;
		if (max < individual[i].get_fitness())
		{
			max = individual[i].get_fitness();
			max_idx = i;
		}
		if (min > individual[i].get_fitness())
		{
			min = individual[i].get_fitness();
			min_idx = i;
		}
		sum_of_fitness += individual[i].get_fitness();
	}
	/*
	// < Maximizing code >
	best_fitness = max;
	best_index = max_idx;
	worst_fitness = min;
	worst_index = min_idx;
	avg_fitness = sum_of_fitness / POP_SIZE;
	*/

	// < Minimizing code >
	best_fitness = min;
	best_index = min_idx;
	worst_fitness = max;
	worst_index = max_idx;
	avg_fitness = sum_of_fitness / POP_SIZE;


	//cout << "Calculation" << endl;
}


//배치의 평가 결과
void assignment_eval(int assignment[K])
{
	int test_data[M][K] = { 0, };
	double time[S];
	double weight = 1;
	double basic_fitness = 1.0;
	double time_weight_fitness = 1.0;
	double eval_basic_sum = 0;
	double eval_weight_sum = 0;

	for (int i = 0; i < K; i++)
	{
		test_data[0][i] = assignment[i];
		//cout << C[i] << endl;
	}
	sim::simulation(test_data, time);

	for (int i = 0; i<K; ++i) {
		//std::cout <<"time: " << time[i] << ',';
		double basic_prod = 1;
		double weight_prod = 1;
		if (time[i] < 2) //2일 내에 끝나는 경우
		{
			weight = PENALTY_1DAY;
		}
		else if (time[i] >= 2 && time[i] < 3) //2~3일
		{
			weight = PENALTY_2DAY;
		}
		else // 3일 넘게 걸리는 경우
		{
			weight = PENALTY_EXCESS;
		}
		basic_prod = sim::SCENARIO[i] * time[i];
		weight_prod = sim::SCENARIO[i] * time[i] * weight;

		eval_basic_sum += basic_prod;
		eval_weight_sum += weight_prod;
		cout << time[i] << endl;
	}
	basic_fitness = eval_basic_sum / 165650;
	time_weight_fitness = eval_weight_sum / 165650;
	cout << "assignment's fitness is: " << basic_fitness << " " << time_weight_fitness << endl;
}
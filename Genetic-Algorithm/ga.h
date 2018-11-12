#include <iostream>
#include <cstdlib>
#include <string>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cassert>
#include <cmath>
#ifndef GA_H
#define GA_H

using namespace std;

const int M = 1; // 무기종류, 예비실험 10
const int K = 16; // 타겟종류,repair 예비실험 12
const int LB = 1;
const int ITERATION = 10;
const int KEEP = 3;
const int TIME_LIMIT = 72;
const int MODE_BRIEF = 0, MODE_DETAIL = 1, MODE_TEST = 2;
const int MODE_C_ERROR = 0, MODE_R_ERROR = 1, MODE_RANGE_ERROR = 2;
const int FLAG_SHORTAGE = 0, FLAG_SAME = 1, FLAG_EXCESS = 2;
const int NON_ESCAPABLE = 0, ESCAPABLE = 1;
const double CROSSOVER_RATE = 0.7;
const double MUTATION_RATE = 0.001;
static const int GENERATIONS = 200;
extern int W[K], C[K], R[M];
extern double P[M], Q[M];			//destroy probability, expected damage 
extern double sum_best_fiteness[GENERATIONS];
extern double sum_avg_fiteness[GENERATIONS];
extern double gene_best_fitness[ITERATION][GENERATIONS];
extern int iteration_num;
extern double random_gene_best_fitness[ITERATION];
//extern double sum_avg_fiteness[GENERATIONS];
void set_input_parameters();

/*
*  0:인천,    1:평택,  2:대산,  3:군산,
*  4:목포,    5:완도,  6:여수,  7:제주,
*  8:서귀포,  9:통영, 10:창원, 11:부산,
* 12:울산,   13:포항, 14:동해, 15:속초
*/
constexpr int S = 16; // # of simulation
extern int SCENARIO[S];

/*
* Maritime Equipment(해상 장비)
* Ground Equipment(육상 장비?)
* Oil-skimmer(유회수기)
*/
const int cur_mar_skimmer[S][1] = {
	569, 814, 594, 600,
	739, 220, 1167, 305,
	309, 287, 474, 1662,
	1617, 450, 505, 200
};
const int cur_gro_skimmer[S][1] = {
	699, 268, 769, 760,
	128, 95, 955, 209,
	414, 630, 314, 472,
	1341, 286, 243, 45
};

const double alpha = 0.2;  // 효율계수
const double beta = 1.0 / 3; // 동원률


class Chromosome
{
public:
	Chromosome();
	void print_chromo();
	void mutate();			// one point flip
	void repair();			// greedy repair
	void evaluate();
	void assert_feasible(int mode);
	int get_data(int row, int col) { return data[row][col]; }
	void put_data(int row, int col, int val) { data[row][col] = val; }
	double get_fitness() { return fitness; }

	void simple_simulation(const int data[M][K], double(&time)[S]);
private:
	int data[M][K];
	double column_sum[K], row_sum[M];
	double fitness;
	void cumulative_sum_calc();
};

class Population
{
public:
	static const int POP_SIZE = 100;

	Population(int pop_size = POP_SIZE);
	void run();
	void statistics_info_calc();
	void print_population(int mode);
	double get_best_fitness() { return best_fitness; }
	double get_avg_fitness() { return avg_fitness; }
	Chromosome& get_best_chromosome() { return individual[best_index]; }

private:
	int best_index, worst_index;
	Chromosome *individual, offspring;
	double best_fitness, avg_fitness, worst_fitness;
	//double gene_best_fitness[6];

	int select();									// roulette wheel
	void crossover(int parent1, int parent2);		// block uniform 
	void replacement(int parent1, int parent2);		// similar survival 
};
#endif // !GA_H

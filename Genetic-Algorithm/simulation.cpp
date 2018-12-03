//
//  simulation.cpp
//  simulation
//
//  Created by jazz4rabbit on 11/9/18.
//

#include <iostream>
#include <algorithm>
#include <cmath>
//#include "ga.h"

#include "simulation.h"


int util::h(int t)
{
	if (t >= constant::T_w_0 && t <= constant::T_w_0 + 24) {
		return std::min(constant::T_w, t - constant::T_w_0);
	}
	else if (t < constant::T_w_0) {
		return h(t + 24) - 10;
	}
	else
		return h(t - 24) + 10;
}

int util::g(int t1, int t2)
{
	return util::h(t1 + t2) - util::h(t2);
}


util::matrix<double, K, K> util::D2T(util::matrix<int, K, K> D, int speed)
{
	util::matrix<double, K, K> T;
	
	for (size_t i = 0; i < K; i++)
		for (size_t j = 0; j < K; j++)
			T[i][j] = D[i][j] / static_cast<double>(speed);

	return T;
}



void sim::simple_simulation(const int data[M][K], double (&time)[S])
{
    for (int i=0; i< S; ++i) {
        auto target = static_cast<double>(SCENARIO[i]) / 3;
        auto oil_per_h = static_cast<double>(data[0][i])*sim::alpha*sim::beta;

		//std::cout << oil_per_h << std::endl;
        if (oil_per_h == 0) {
            time[i] = std::numeric_limits<double>::infinity();
            continue;
        }
        int t = 0;
        while (target > 0) {
            target -= oil_per_h;
            t++;
        }

        time[i] = t/8.0;
    }
    return ;
}

void sim::simulation(const int data[M][K], double(&time)[S])
{
	using util::g;
	using util::h;
	using util::D2T;

	int tick_count;
	auto T = D2T(constant::D, constant::ship_speed);

	for (size_t i = 0; i < S; i++)
	{
		tick_count = 0;
		double Q_i = SCENARIO[i];
		auto target = Q_i / (constant::alpha * constant::beta * constant::gamma) / 3.0;

		std::array<int, K> work_times;
		// get work_times
		for (size_t j = 0; j < K; j++) {
			work_times[j] = -g(T[i][j], constant::T_0 + T[i][i]);
		}

		while (true) {
			tick_count++;
			int work_acc = 0;
			for (size_t j = 0; j < K; j++) {
				work_times[j]++;
				work_acc += std::max(work_times[j], 0) * data[0][j];
			}
			if (target <= work_acc) {
				time[i] = std::ceil(tick_count * constant::gamma) / 8.0;
				break;
			}
		}
	}

	return;
}


/*
int main() {
    int data[M][1];
    double time[S];

    for (int i=0; i<M; ++i) {
        data[i][0] = cur_mar_skimmer[i][0] + cur_gro_skimmer[i][0];
        std::cout << data[i][0] << ',';
    }
    std::cout << std::endl;

    simple_simulation(data, time);
    
    for (int i=0; i<M; ++i) {
        std::cout << time[i] << ',';
    }
    std::cout << std::endl;
}
*/

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

        int q = t / 8;
        int r = t % 8;

        // 하루의 일할 수 있는 시간이 8시간이기 때문에
        // 일한 시간이 9시간인 경우, (9/8)*24 + 9%8 로 계산
        time[i] = r != 0 ? q*24 + r : q*24 - 16;
    }
    return ;
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

//
//  simulation.h
//  simulation
//
//  Created by jazz4rabbit on 11/9/18.
//

#include "ga_parameters.h"

#ifndef SIMULATION_H
#define SIMULATION_H

namespace sim {

	/*
	 *  0:인천,    1:평택,  2:대산,  3:군산,
	 *  4:목포,    5:완도,  6:여수,  7:제주,
	 *  8:서귀포,  9:통영, 10:창원, 11:부산,
	 * 12:울산,   13:포항, 14:동해, 15:속초
	 */
	constexpr int S = 16; // # of simulation
	const int SCENARIO[S] = {
		8500,  1200, 45000, 3800,
		8500,   600, 45000,  800,
		500,   1700,  1200, 2500,
		45000,  800,   500,   50
	};

	/*
	 * Maritime Equipment(해상 장비)
	 * Ground Equipment(육상 장비?)
	 * Oil-skimmer(유회수기)
	 */
	constexpr int cur_mar_skimmer[S][1] = {
		569, 814, 594, 600,
		739, 220, 1167, 305,
		309, 287, 474, 1662,
		1617, 450, 505, 200
	};
	constexpr int cur_gro_skimmer[S][1] = {
		699, 268, 769, 760,
		128, 95, 955, 209,
		414, 630, 314, 472,
		1341, 286, 243, 45
	};

	constexpr double alpha = 0.2;  // 효율계수
	constexpr double beta = 1.0 / 3; // 동원률


	/*
	 * void simple_simulation(int data[M][K], doulbe (&time)[K])
	 *
	 * S개의 시나리오와 data로 부터 time을 산출
	 *
	 * (제약조건:
	 *   K가 무기의 종류지만, 우리는 0번째 무기만을
	 *   이용함. 오전 7시에 사건 발생, 하루의 일할 수
	 *   있는 시간은 8시간이고, 사건 발생 지점에 할당된
	 *   유회수기만을 사용하여 1시간 단위로 시뮬레이션
	 * )
	 * -----------------------------------------------
	 * Global Varaible
	 *
	 * int SCENARIO[S] : 유출유 사고 시나리오
	 *                   (index:지역, value: 유출유(톤))
	 * double alpha : 0.2 효율 계수
	 * double beta  : 1/3 동원률
	 *
	 * -----------------------------------------------
	 * Function Parameter
	 *
	 * int data[M][K]
	 *  : M개의 지점에 첫 번째 무기(유회수기)가 배치된
	 *    상태.
	 *
	 * (&time)[S]
	 *  : 시나리오로 부터 유회수기를 이용하여 유출유를
	 *    회수하는데 걸리는 시간을 산출
	 * -----------------------------------------------
	 *
	 * 예제 :
	 * // 첫 번째 지역에 사고 1톤,
	 * // 두 번째 지역에 사고 1톤,
	 * // ...
	 * // 네 번째 지역에 사고 4톤
	 * int SCENARIO[S] = { 1, 1, 2, 1 }
	 *
	 * // 첫 번째 지역에 유회수기 3개 배치,
	 * // ...
	 * // 네 번째 지역에 유회수기 1개 배치
	 * int data[4][1] = { {3},
	 *                    {6},
	 *                    {6},
	 *                    {1}
	 *                  }
	 *
	 * double time[4];
	 * simple_simulation(data, time);
	 *
	 * time == { 5, 3, 5, 31 }
	 *
	 */
	void simple_simulation(const int data[M][K], double(&time)[S]);

}
#endif /* SIMULATION_H */


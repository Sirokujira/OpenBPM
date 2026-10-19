/*
outputChars.c

calculate and output to oth.log
*/

#include "obpm.h"
#include "user_define_complex.h"
#include "obpm_prototype.h"


// farfield = 1 : 遠方界用の閉曲面 (Surface*) を準備する。
//   obpm.out (writeout) が Surface 配列を書き出すため、FDTD ソルバーと
//   BPM の obpm.out 出力時は 1 を渡す。BPM で -no-fdtd-out のときは
//   Surface を誰も参照しないので 0 を渡して準備を省く
//   (全境界セルでの節点材料計算 NodeE_c/NodeH_c が不要になる)。
void outputChars(FILE *fp, int farfield)
{
	// setup far field
	if (farfield) {
		alloc_farfield();
		setup_farfield();
	}

	// input imepedanece
	if (NFeed && NFreq1) {
		calcZin();
		outputZin(fp);
	}

	// input power (writeout は NFeed * NFreq2 個を書き出すため NFreq1 とは独立に確保する)
	if (NFeed && NFreq2) {
		calcPin();  // for post
	}

	// S-parameters
	if (NPoint && NFreq1) {
		calcSpara();
		outputSpara(fp);
	}

	// coupling
	if (NFeed && NPoint && NFreq1) {
		outputCoupling(fp);
	}

	// cross section
	if (IPlanewave && NFreq2) {
		outputCross(fp);
	}
}

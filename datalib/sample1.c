/*
sample1.c

OpenBPMデータ作成ライブラリ、サンプルプログラム No.1

コンパイル+実行:
Windows + VC++:
> cl.exe /O2 sample1.c obpm_datalib.c
> sample1.exe
Linux + gcc:
$ gcc -O sample1.c obpm_datalib.c -o sample1
$ ./sample1
*/

#include "obpm_datalib.h"

int main(void)
{
	// initialize

	obpm_init();

	// title

	obpm_title("sample1");

	// mesh

	obpm_xsection(2, -75e-3, +75e-3);
	obpm_xdivision(1, 30);

	obpm_ysection(2, -75e-3, +75e-3);
	obpm_ydivision(1, 30);

	obpm_zsection(4, -75e-3, -25e-3, +25e-3, +75e-3);
	obpm_zdivision(3, 10, 11, 10);

	// material

	obpm_material(2.0, 0.0, 1.0, 0.0, "");

	// geometry

	obpm_geometry(1, 1, 0e-3, 0e-3, 0e-3, 0e-3, -25e-3, +25e-3);

	// feed

	obpm_feed('Z', 0e-3, 0e-3, 0e-3, 1, 0, 50);
	//obpm_rfeed(10);

	// ABC

	//obpm_pml(5, 2, 1e-5);

	// frequency

	obpm_frequency1(2e9, 3e9, 10);
	obpm_frequency2(3e9, 3e9, 0);

	// solver

	obpm_solver(1000, 100, 1e-3);

	// iteration

	obpm_plotiter(1);

	// waveform and spectrum

	//obpm_plotfeed(1);
	//obpm_plotpoint(1);

	// frequency

	obpm_plotsmith(1);
	obpm_plotzin(1, 0, 0, 0);
	obpm_plotyin(1, 0, 0, 0);
	obpm_plotref(1, 0, 0, 0);
	//obpm_plotfar0d(90, 0, 1, 0, 0, 0);

	// far-1d

	obpm_plotfar1d('X', 72, 0);

	// far-2d

	obpm_plotfar2d(18, 36);

	// near-1d

	obpm_plotnear1d("E", 'Z', 30e-3, 0e-3);

	// near-2d

	obpm_plotnear2d("E", 'X', 30e-3);

	// output options

	//obpm_far1dstyle(1);
	//obpm_far1dcomponent(1, 0, 0);
	//obpm_far1ddb(1);
	//obpm_far1dnorm();
	//obpm_far1dscale(-30, +10, 4);
	//obpm_far2dcomponent(1, 0, 0, 0, 0, 0, 0);
	//obpm_far2ddb(1);
	//obpm_far2dscale(-30, +10);
	//obpm_far2dobj(0.5);
	//obpm_near1ddb(1);
	//obpm_near1dscale(-30, +10, 4);
	//obpm_near2ddim(1, 1);
	//obpm_near2ddb(1);
	//obpm_near2dscale(-30, +10);
	//obpm_near2dcontour(0);
	//obpm_near2dobj(1);
	//obpm_near2dzoom(-50e-3, 50e-3, -50e-3, 50e-3);
	//obpm_window2d(750, 500, 15, 0);
	//obpm_window3d(600, 600, 12, 60, 30);

	// output

	obpm_outdata("sample1.obpm");

	return 0;
}

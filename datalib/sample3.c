/*
sample3.c

OpenBPMデータ作成ライブラリ、サンプルプログラム No.3

コンパイル+実行:
Windows + VC++:
> cl.exe /O2 sample3.c obpm_datalib.c
> sample3.exe
Linux + gcc:
$ gcc -O sample3.c obpm_datalib.c -o sample3
$ ./sample3
*/

#include "obpm_datalib.h"

int main(void)
{
	// initialize

	obpm_init();

	// title

	obpm_title("sample3");

	// mesh

	obpm_xsection(2, -150e-3, +150e-3);
	obpm_xdivision(1, 30);

	obpm_ysection(2, -150e-3, +150e-3);
	obpm_ydivision(1, 30);

	obpm_zsection(2, -150e-3, +150e-3);
	obpm_zdivision(1, 30);

	// material

	obpm_material(2.0, 0.0, 1.0, 0.0, "");

	// geometry

	obpm_geometry(1, 1, -50e-3, +50e-3, -50e-3, +50e-3, -50e-3, +50e-3);

	// planewave

	obpm_planewave(90, 0, 1);

	// ABC

	//obpm_pml(5, 2, 1e-5);

	// frequency

	obpm_frequency1(3e9, 3e9, 0);
	obpm_frequency2(3e9, 3e9, 0);

	// solver

	obpm_solver(1000, 50, 1e-3);

	// iteration

	obpm_plotiter(1);

	// far-0d

	//obpm_plotfar0d(90, 0);

	// far-1d

	obpm_plotfar1d('Z', 72, 0);

	// far-2d

	obpm_plotfar2d(18, 36);

	// near-1d

	obpm_plotnear1d("E", 'X', 0e-3, 0e-3);

	// near-2d

	obpm_plotnear2d("E", 'Z', 0e-3);

	// output options

	//obpm_far1dstyle(1);
	//obpm_far1dcomponent(1, 0, 0);
	//obpm_far1ddb(1);
	//obpm_far1dscale(-30, +10, 4);
	//obpm_far2dcomponent(1, 0, 0, 0, 0, 0,0);
	//obpm_far2ddb(1);
	//obpm_far2dscale(-30, +10);
	//obpm_far2dobj(0.5);
	//obpm_near1ddb(1);
	//obpm_near1dscale(-30, +10, 4);
	//obpm_near1dnoinc();
	//obpm_near2ddim(1, 1);
	//obpm_near2ddb(1);
	//obpm_near2dscale(-30, +10);
	//obpm_near2dcontour(0);
	//obpm_near2dobj(1);
	//obpm_near2dnoinc();

	// output

	obpm_outdata("sample3.obpm");

	return 0;
}

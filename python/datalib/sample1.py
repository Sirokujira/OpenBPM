# -*- coding: utf-8 -*-
"""
sample1.py
OpenBPM データ作成ライブラリ (Python) サンプルプログラム (1)

使い方:
(1) 本ファイルを編集し、obpm_datalib.py と同じフォルダにおく
(2) > python sample1.py
(3) OpenBPM入力ファイル(sample1.obpm)が出力される
"""

import obpm_datalib as obpm

obpm.title('dipole antenna')

obpm.xsection(-0.075, 30, +0.075)
obpm.ysection(-0.075, 30, +0.075)
obpm.zsection(-0.075, 10, -0.025, 11, 0.025, 10, 0.075)

#obpm.material(2.0, 0.0, 1.0, 0.0)

obpm.geometry(1, 1, 0.0, 0.0, 0.0, 0.0, -0.025, +0.025)

obpm.feed('Z', 0.0, 0.0, 0.0, 1.0, 0.0, 50)
#obpm.planewave(90, 0, 1)
#obpm.pml(5, 2, 1e-5)

obpm.frequency1(2e9, 3e9, 10)
obpm.frequency2(3e9, 3e9, 0)

obpm.solver(1000, 50, 1e-3)

obpm.plotiter(1)
obpm.plotzin(1)
obpm.plotyin(1)
obpm.plotref(1)

obpm.plotfar1d('X', 90)
#obpm.far1dstyle(2)
#obpm.far1dcomponent(1, 1, 1)
#obpm.far1ddb(1)
#obpm.far1dnorm(1)
#obpm.far1dscale(-30, 10, 4)

obpm.plotfar2d(18, 36)
#obpm.far2dcomponent(1, 0, 0, 0, 0, 0, 0)
#obpm.far2ddb(1)
#obpm.far2dscale(-20, 10, 6)
#obpm.far2dobj(0.5)

obpm.plotnear1d('E', 'Z', 0.03, 0.0)
#obpm.near1ddb(1)
#obpm.near1dnoinc(1)
#obpm.near1dscale(-30, 20, 5)

obpm.plotnear2d('E', 'X', 0.03)
#obpm.near2ddim(1, 1)
#obpm.near2dframe(20)
#obpm.near2ddb(1)
#obpm.near2dscale(-40, +20)
#obpm.near2dcontour(1)
obpm.near2dobj(1)
#obpm.near2dnoinc(1)
#obpm.near2dzoom(-0.1, +0.1, -0.1, +0.1)

obpm.output('sample1.obpm')

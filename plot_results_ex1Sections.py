import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt


#  0: disp_x_1
#  1: disp_y_1
#  2: disp_z_1
#  3:  rot_x_1
#  4:  rot_y_1
#  5:  rot_z_1
#  6: disp_x_2
#  7: disp_y_2
#  8: disp_z_2
#  9:  rot_x_2
# 10:  rot_y_2
# 11:  rot_z_2
# 12: forc_x_1
# 13: forc_y_1
# 14: forc_z_1
# 15:  mom_x_1
# 16:  mom_y_1
# 17:  mom_z_1
# 18: forc_x_2
# 19: forc_y_2
# 20: forc_z_2
# 21:  mom_x_2
# 22:  mom_y_2
# 23:  mom_z_2


data03 = np.loadtxt("results_ex1_3sections_0_4.dat", usecols=(10, 16), max_rows=9)
data06 = np.loadtxt("results_ex1_6sections_0_01.dat", usecols=(10, 16), max_rows=360)
data10 = np.loadtxt("results_ex1_10sections_0_01.dat", usecols=(10, 16), max_rows=360)
data20 = np.loadtxt("results_ex1_20sections_0_004.dat", usecols=(10, 16), max_rows=900)

colors = {
	"TUM_blue1" : (0/255., 82/255., 147/255.),
	"TUM_blue2" : (0/255., 51/255., 89/255.),
	"TUM_grey1" : (88/255., 88/255., 90/255.),
	"TUM_grey2" : (156/255., 157/255., 159/255.),
	"TUM_grey3" : (217/255., 218/255., 219/255.),
	"TUM_ivory" : (218/255., 215/255., 203/255.),
	"TUM_orange" : (227/255., 114/255., 34/255.),
	"TUM_green" : (162/255., 173/255., 0/255.),
	"TUM_bluel1" : (152/255., 198/255., 234/255.),
	"TUM_bluel2" : (100/255., 160/255., 200/255.),
}

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(-data03[:,0]*1e6/40, data03[:,1], "--", color=colors["TUM_grey2"], label="3 sections")
ax.plot(-data06[:,0]*1e6/40, data06[:,1], ":",  color=colors["TUM_orange"], label="6 sections")
ax.plot(-data10[:,0]*1e6/40, data10[:,1], "-*",  color=colors["TUM_green"], label="10 sections", markevery=20)
ax.plot(-data20[:,0]*1e6/40, data20[:,1], "-",  color=colors["TUM_blue1"], label="20 sections")
ax.set(
	title=r"Moment - Curvature",
	xlabel=r"$\phi_y [\mu~rad/in]$",
	ylabel=r"$M_y~[kip\cdot in]$"
)
ax.grid()
ax.legend()
plt.show()

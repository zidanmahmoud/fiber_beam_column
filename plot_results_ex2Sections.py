
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import numpy as np

# data2 = np.loadtxt("results_ex2_2sections_0_1.dat", usecols=(8, 20), max_rows=10)
# data3 = np.loadtxt("results_ex2_3sections_0_1.dat", usecols=(8, 20), max_rows=10)
data4 = np.loadtxt("results_ex2_4sections0_01.dat")
# data6 = np.loadtxt("results_ex2_6sections0_005.dat")

colors = {
	"TUM_blue1" : (0/255, 82/255, 147/255),
	"TUM_blue2" : (0/255, 51/255, 89/255),
	"TUM_grey1" : (88/255, 88/255, 90/255),
	"TUM_grey2" : (156/255, 157/255, 159/255),
	"TUM_grey3" : (217/255, 218/255, 219/255),
	"TUM_ivory" : (218/255, 215/255, 203/255),
	"TUM_orange" : (227/255, 114/255, 34/255),
	"TUM_green" : (162/255, 173/255, 0/255),
	"TUM_bluel1" : (152/255, 198/255, 234/255),
	"TUM_bluel2" : (100/255, 160/255, 200/255),
}

fig = plt.figure()
ax = fig.add_subplot(111)
# ax.plot(data2[:,0], data2[:,1], "--", color=colors["TUM_grey1"], label="2 sections")
# ax.plot(data3[:,0], data3[:,1], "-.", color=colors["TUM_orange"], label="3 sections")
# ax.plot(data4[:,0], data4[:,1], "-*", color=colors["TUM_grey2"], label="4 sections", markevery=10)
ax.plot(data4[:,0], data4[:,1], "-" , color=colors["TUM_blue1"], label="6 sections")
ax.set(
	title=r"Load - Displacement",
	xlabel=r"Displacement in Z direction [in]",
	ylabel=r"$\lambda$"
)
ax.grid()
# ax.legend()
plt.show()
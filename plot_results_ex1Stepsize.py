import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import numpy as np

data04 = np.loadtxt("results_ex1_4sections_0_4.dat", usecols=(10, 16))
data01 = np.loadtxt("results_ex1_4sections_0_1.dat", usecols=(10, 16))
data005 = np.loadtxt("results_ex1_4sections_0_05.dat", usecols=(10, 16))

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
ax.plot(-data04[:,0]*1e6/40, data04[:,1], "-*", color=colors["TUM_orange"], label="0.4")
ax.plot(-data01[:,0]*1e6/40, data01[:,1], "--", color=colors["TUM_grey1"], label="0.1")
ax.plot(-data005[:,0]*1e6/40, data005[:,1], "-", color=colors["TUM_blue1"], label="0.05")
ax.set(
	title=r"Load - Displacement",
	xlabel=r"Displacement in Z direction [in]",
	ylabel=r"$\lambda$"
)
ax.grid()
ax.legend()
plt.show()


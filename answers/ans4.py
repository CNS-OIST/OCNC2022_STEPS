# Uncomment the line below and run the cell to reload the empty exercise
#%load ./answers/blank4.py

%matplotlib inline
import matplotlib.pyplot as plt

plt.plot(counts_wm.time[0], counts_wm.data[0])
plt.ylabel('Number of molecules')
plt.xlabel('Time [s]')
plt.legend(counts_wm.labels)
plt.show()

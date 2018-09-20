# Numpy is the numerical library dadi is built upon
from numpy import array

import dadi


data = dadi.Spectrum.from_file('w_c_all.fs')

import pylab
pylab.figure(1)
dadi.Plotting.plot_single_2d_sfs(data, vmin=1,pop_ids =('P. monodii','P. glaucum'))
# This ensures that the figure pops up. It may be unecessary if you are using
# ipython.
pylab.show()
# Save the figure
pylab.savefig('Wild_Cultivated_all.png', dpi=50)



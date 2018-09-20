"""
Adapdated by Abdou Rahmane WADE
"""
# Numpy is the numerical library dadi is built upon
from numpy import array

import dadi

# In demographic_models.py, we've defined a custom model for this problem
import Model2

# Load the data
data = dadi.Spectrum.from_file('w_c_centre.fs')
ns = data.sample_sizes

# These are the grid point settings will use for extrapolation.
pts_l = [40,50,60]

# The Demographics1D and Demographics2D modules contain a few simple models,
# mostly as examples. We could use one of those.
func = dadi.Demographics2D.split_mig
# Instead, we'll work with our custom model
func = Model2.prior_onegrow_mig_ancestralfixed

# Now let's optimize parameters for this model.

# The upper_bound and lower_bound lists are for use in optimization.
# Occasionally the optimizer will try wacky parameter values. We in particular
# want to exclude values with very long times, very small population sizes, or
# very high migration rates, as they will take a long time to evaluate.
# Parameters are: (nu1F, nu2B, nu2F, m12, m21, TB, TE,e)
upper_bound = [100, 100, 100, 10, 10, 3, 3, 0.1]
lower_bound = [1e-2, 1e-2, 1e-2, 0, 0, 0, 0, 0]

# This is our initial guess for the parameters, which is somewhat arbitrary.
p0 = [2,0.1,2,1,1,0.2,0.2,0.001]
# Make the extrapolating version of our demographic model function.
func_ex = dadi.Numerics.make_extrap_log_func(func)

# Perturb our parameters before optimization. This does so by taking each
# parameter a up to a factor of two up or down.
p0 = dadi.Misc.perturb_params(p0, fold=1, upper_bound=upper_bound,
                              lower_bound=lower_bound)
 
# Do the optimization. By default we assume that theta is a free parameter,
# since it's trivial to find given the other parameters. If you want to fix
# theta, add a multinom=False to the call.
# The maxiter argument restricts how long the optimizer will run. For real 
# runs, you will want to set this value higher (at least 10), to encourage
# better convergence. You will also want to run optimization several times
# using multiple sets of intial parameters, to be confident you've actually

popt = [0.172149   ,  0.0302491  ,  0.320733   ,  6.24472    ,  2.46278    ,  0.0087688  ,  0.603221   ,  0.0709079]
print('Best-fit parameters: {0}'.format(popt))

# Calculate the best-fit model AFS.
model = func_ex(popt, ns, pts_l)
# Likelihood of the data given the model AFS.
ll_model = dadi.Inference.ll_multinom(model, data)
print('Maximum log composite likelihood: {0}'.format(ll_model))
# The optimal value of theta given the model.
theta = dadi.Inference.optimal_sfs_scaling(model, data)
print('Optimal value of theta: {0}'.format(theta))

# Plot a comparison of the resulting fs with the data.
import pylab
pylab.figure(1)
dadi.Plotting.plot_2d_comp_multinom(model, data, vmin=1, resid_range=50,
                                    pop_ids =('P. monodii','P. glaucum'))
# This ensures that the figure pops up. It may be unecessary if you are using
# ipython.
pylab.show()
# Save the figure
pylab.savefig('ModelMigUnifBotTimeNoAncError_centre.png', dpi=50)

# Let's generate some data using ms, if you have it installed.
mscore = ModelMigUnifBotTimeNoAncError.prior_onegrow_mig_mscore(popt)
# I find that it's most efficient to simulate with theta=1, average over many
# iterations, and then scale up.
mscommand = dadi.Misc.ms_command(1., ns, mscore, int(1e5))
# If you have ms installed, uncomment these lines to see the results.

# We use Python's os module to call this command from within the script.
import os
return_code = os.system('{0} > test.msout'.format(mscommand))
# We check the return code, so the script doesn't crash if you don't have ms
# installed
if return_code == 0:
    msdata = dadi.Spectrum.from_ms_file('test.msout')
    pylab.figure(2)
    dadi.Plotting.plot_2d_comp_multinom(model, theta*msdata, vmin=1,
                                        pop_ids=('P. monodii','P. glaucum'))
    pylab.show()

# Estimate parameter uncertainties using the Godambe Information Matrix, to
# account for linkage in the data. To use the GIM approach, we need to have
# spectra from bootstrapping our data.  Let's load the ones we've provided for
# the example.  
# (We're using Python list comprehension syntax to do this in one line.)
all_boot = [dadi.Spectrum.from_file('bootstraps/{0:02d}.fs'.format(ii)) 
            for ii in range(100)]
uncerts = dadi.Godambe.GIM_uncert(func_ex, pts_l, all_boot, popt, data, 
                                  multinom=True)
# uncert contains the estimated standard deviations of each parameter, with
# theta as the final entry in the list.
print('Estimated parameter standard deviations from GIM: {0}'.format(uncerts))

# For comparison, we can estimate uncertainties with the Fisher Information
# Matrix, which doesn't account for linkage in the data and thus underestimates
# uncertainty. (Although it's a fine approach if you think your data is truly
# unlinked.)
uncerts_fim = dadi.Godambe.FIM_uncert(func_ex, pts_l, popt, data, multinom=True)
print('Estimated parameter standard deviations from FIM: {0}'.format(uncerts_fim))

print('Factors by which FIM underestimates parameter uncertainties: {0}'.format(uncerts/uncerts_fim))

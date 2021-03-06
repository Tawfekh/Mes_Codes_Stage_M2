"""
Adapdated by Abdou Rahmane WADE
"""
import numpy
import dadi

def prior_onegrow_mig_ancestralfixed((nu1F, nu2B, nu2F, m12, m21, TB, TE,e), (n1,n2), pts):
    """
    Model with growth, split, bottleneck in pop2 , exp recovery, asymmetric migration rate

    nu1F: The final size for pop1
    nu2B: The bottleneck size for pop2
    nu2F: The final size for pop2
    TB: The time between the split and end of bottleneck
    TE: The end of bottleneck and present
    m12: rate of migration from pop 2 into pop 1
    m12: rate of migration from pop 1 into pop 2
    e: polarization error frac

    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = yy = dadi.Numerics.default_grid(pts)

    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # The divergence
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    phi = dadi.Integration.two_pops(phi, xx, TB, nu1=nu1F, nu2=nu2B)
    # We need to define a function to describe the non-constant population 2
    # size. lambda is a convenient way to do so.
    nu2_func = lambda t: nu2B*(nu2F/nu2B)**(t/TE)
    phis = dadi.Integration.two_pops(phi, xx, TE, nu1=nu1F, nu2=nu2_func,
                                    m12=m12, m21=m21)

    phie= dadi.Numerics.reverse_array(phis)
    phi=(1-e)*phis+e*phie
    sfs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,yy))
    return sfs

import dadi
import numpy as np
import Optimize_Functions
from Find_Best_Projections import dadi_test_projections

## import data
## example for N-Coastal
dd = dadi.Misc.make_data_dict_vcf("./0miss_vcfs/CEN_str_full_06_thin.recode.vcf", "popfile_cen_str.txt")
print(len(dd))

## best projections
dadi_test_projections.run_projections(dd, pop_ids=["CEN"], polarized=False, maxproj=[40], min_frac=0.1)
## best projection: S-Coastal: 16, N-Coastal: 40, S-Oceanic: 16

fs = dadi.Spectrum.from_data_dict(dd, ["CEN"], projections=[40], polarized=False)
print(fs)

## single population summary statistics

thetaW = fs.Watterson_theta()  ## Watterson's theta
print(thetaW)

pi = fs.pi()  ## expected heterozygosity!!
print(pi)

D = fs.Tajima_D()
print(D)


## try to fold the spectrum, but spectrum is already folded

# folded = fs.fold()
# print(folded)

## first model: standard neutral model

def snm(notused, ns, pts):
    """
    Standard neutral model.

    ns = (n1,)

    n1: Number of samples in resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx,))
    return fs


## second model: one sudden change (two epoch model):

def two_epoch(params, ns, pts):
    """
    Instantaneous size change some time ago.

    params = (nu,T)
    ns = (n1,)

    nu: Ratio of contemporary to ancient population size
    T: Time in the past at which size change happened (in units of 2*Na
       generations)
    n1: Number of samples in resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu, T = params

    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)

    phi = dadi.Integration.one_pop(phi, xx, T, nu)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx,))
    return fs

## third model: two sudden changes (three epoch model):

def three_epoch(params, ns, pts):
    """
    params = (nuB,nuF,TB,TF)
    ns = (n1,)

    nuB: Ratio of bottleneck population size to ancient pop size
    nuF: Ratio of contemporary to ancient pop size
    TB: Length of bottleneck (in units of 2*Na generations)
    TF: Time since bottleneck recovery (in units of 2*Na generations)

    n1: Number of samples in resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nuB, nuF, TB, TF = params

    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)

    phi = dadi.Integration.one_pop(phi, xx, TB, nuB)
    phi = dadi.Integration.one_pop(phi, xx, TF, nuF)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx,))
    return fs

## fouth model: one sudden change followed by exponential growth/size reduction (bottlegrowth_model):

def bottlegrowth(params, ns, pts):
    """
    Instantanous size change followed by exponential growth.

    params = (nuB,nuF,T)
    ns = (n1,)

    nuB: Ratio of population size after instantanous change to ancient
         population size
    nuF: Ratio of contemporary to ancient population size
    T: Time in the past at which instantaneous change happened and growth began
       (in units of 2*Na generations)
    n1: Number of samples in resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nuB,nuF,T = params

    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)

    nu_func = lambda t: nuB*np.exp(np.log(nuF/nuB) * t/T)
    phi = dadi.Integration.one_pop(phi, xx, T, nu_func)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx,))
    return fs


## save functions in variable names
my_demo_func = snm
my_two_epochs_func = two_epoch
my_three_epochs_func = three_epoch
my_bottle_func = bottlegrowth

## setting up initial values

pts = np.array([50, 60, 70])
ns = [40]
neutral_params = np.array([])

## extrapolate
func_ex = dadi.Numerics.make_extrap_log_func(my_demo_func)
func_two_epoch_ex = dadi.Numerics.make_extrap_log_func(my_two_epochs_func)
func_three_epoch_ex = dadi.Numerics.make_extrap_log_func(my_three_epochs_func)
func_bottle_ex = dadi.Numerics.make_extrap_log_func(my_bottle_func)


## create models
model0 = func_ex(neutral_params, ns, pts)

## Calculate log-likelihood for neutral model
ll_snm = dadi.Inference.ll_multinom(model0, fs)
print('Log-likelihood: %s' % (ll_snm))

## Calculate AIC for neutral model

aic_snm = 0 - (ll_snm * 2)
print('AIC: %s' % (aic_snm))


## Optimize Function

reps = [10,20,30]
maxiters = [5,10,50]
folds = [3,2,1] ## more conservative


## write a loop that runs the command many times

for i in range(0,5):
    var_name = f"TM_filter_{i}_opt"
    print(var_name)
    Optimize_Functions.Optimize_Routine(fs, pts, var_name, "two_epoch", my_two_epochs_func, 3, 2, fs_folded=True,
                                        reps=reps, maxiters=maxiters, folds=folds)
    Optimize_Functions.Optimize_Routine(fs, pts, var_name, "three_epoch", my_three_epochs_func, 3, 4, fs_folded=True,
                                        reps=reps, maxiters=maxiters, folds=folds)
    Optimize_Functions.Optimize_Routine(fs, pts, var_name, "bottle", my_bottle_func, 3, 3, fs_folded=True,
                                        reps=reps, maxiters=maxiters, folds=folds)
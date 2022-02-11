import matplotlib.pyplot as pyplot
import pylab
import dadi
import numpy as np

## create a pyplot figure
pyplot.figure()

## import empirical SFS
## example for S-Oceanic
dd = dadi.Misc.make_data_dict_vcf("./0miss_vcfs/TM_full_06_thin.recode.vcf", "popfile_tm.txt")
print(len(dd))

## best projection: S-Coastal: 16, N-Coastal: 40, S-Oceanic: 16
fs = dadi.Spectrum.from_data_dict(dd, ["TM"], projections=[16], polarized=False)
print(fs)

## create a simulated SFS

## create functions as specified by dadi

my_demo_func = dadi.Demographics1D.snm
my_two_epochs_func = dadi.Demographics1D.two_epoch
my_three_epochs_func = dadi.Demographics1D.three_epoch
my_bottle_func = dadi.Demographics1D.bottlegrowth

## extrapolate

func_ex = dadi.Numerics.make_extrap_log_func(my_demo_func)
func_two_epoch_ex = dadi.Numerics.make_extrap_log_func(my_two_epochs_func)
func_three_epoch_ex = dadi.Numerics.make_extrap_log_func(my_three_epochs_func)
func_bottle_ex = dadi.Numerics.make_extrap_log_func(my_bottle_func)

## set up parameters (best parameters calculated by the dadi likelihood optimizations runs in single_pop_filtered.py)

neutral_params = []
two_epoch_params = [1.387,0.4142] ## N-Coastal: 7.3044,0.5709, S-Coastal: 4.167,0.7319
three_epoch_params = [2.209,3.2139,3.0023,0.0532] ## N-Coastal: 3.1664,11.8915,0.8239,0.3645, S-Coastal: 4.7206,1.3716,0.6469,0.0128
bottlegrowth_params = [0.8562,1.4597,1.3754] ## N-Coastal: 0.9933,12.0659,1.1209, S-Coastal: 9.8342,3.1371,0.6226

## grid size and number of seqs

pts = np.array([40, 50, 60])
ns = [16]

## simulate SFS based on model parameters

model0 = func_ex(neutral_params, ns, pts)
model_two_epoch = func_two_epoch_ex(two_epoch_params, ns, pts)
model_three_epoch = func_three_epoch_ex(three_epoch_params, ns, pts)
model_bottlegrowth = func_bottle_ex(bottlegrowth_params, ns, pts)


## scale optimally

model0 = dadi.Inference.optimally_scaled_sfs(model0, fs)
model_bottlegrowth = dadi.Inference.optimally_scaled_sfs(model_bottlegrowth, fs)
model_three_epoch = dadi.Inference.optimally_scaled_sfs(model_three_epoch, fs)
model_two_epoch = dadi.Inference.optimally_scaled_sfs(model_two_epoch, fs)


## plots
def plot_1d_comp_Poisson(model, data, fig_num=None, residual='Anscombe',
                         plot_masked=False, show=True):
    """
    Poisson comparison between 1d model and data.


    model: 1-dimensional model SFS
    data: 1-dimensional data SFS
    fig_num: Clear and use figure fig_num for display. If None, an new figure
             window is created.
    residual: 'Anscombe' for Anscombe residuals, which are more normally
              distributed for Poisson sampling. 'linear' for the linear
              residuals, which can be less biased.
    plot_masked: Additionally plots (in open circles) results for points in the
                 model or data that were masked.
    show: If True, execute pylab.show command to make sure plot displays.
    """
    if fig_num is None:
        f = pylab.gcf()
    else:
        f = pylab.figure(fig_num, figsize=(7,7))
    pylab.clf()

    if data.folded and not model.folded:
        model = model.fold()

    masked_model, masked_data = dadi.Numerics.intersect_masks(model, data)

    ax = pylab.subplot(2,1,1)
    pylab.semilogy(masked_data, '-ob', label='data')
    pylab.semilogy(masked_model, '-or', label='model')

    if plot_masked:
        pylab.semilogy(masked_data.data, '--ob', mfc='w', zorder=-100)
        pylab.semilogy(masked_model.data, '--or', mfc='w', zorder=-100)

    ax.legend(loc='upper right')

    pylab.subplot(2,1,2, sharex = ax)
    if residual == 'Anscombe':
        resid = dadi.Inference.Anscombe_Poisson_residual(masked_model, masked_data)
    elif residual == 'linear':
        resid = dadi.Inference.linear_Poisson_residual(masked_model, masked_data)
    else:
        raise ValueError("Unknown class of residual '%s'." % residual)
    pylab.plot(resid, '-og')
    if plot_masked:
        pylab.plot(resid.data, '--og', mfc='w', zorder=-100)

    ax.set_xlim(0, data.shape[0]-1)
    if show:
        pylab.show()


def plot_1d_comp_multinom(model, data, fig_num=None, residual='Anscombe',
                          plot_masked=False):
    """
    Mulitnomial comparison between 1d model and data.


    model: 1-dimensional model SFS
    data: 1-dimensional data SFS
    fig_num: Clear and use figure fig_num for display. If None, an new figure
             window is created.
    residual: 'Anscombe' for Anscombe residuals, which are more normally
              distributed for Poisson sampling. 'linear' for the linear
              residuals, which can be less biased.
    plot_masked: Additionally plots (in open circles) results for points in the
                 model or data that were masked.

    This comparison is multinomial in that it rescales the model to optimally
    fit the data.
    """
    model = dadi.Inference.optimally_scaled_sfs(model, data)

    plot_1d_comp_Poisson(model, data, fig_num, residual,
                         plot_masked)



#dadi.Plotting.plot_1d_fs(fs)
#plot_1d_comp_Poisson(model0, fs)
#plot_1d_comp_Poisson(model_bottlegrowth, fs)
#plot_1d_comp_Poisson(model_three_epoch, fs)
#plot_1d_comp_Poisson(model_two_epoch, fs)


pyplot.show()
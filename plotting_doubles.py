import matplotlib.pyplot as pyplot
import pylab
import dadi
import numpy as np
from Two_Population_Pipeline import Models_2D

## create a pyplot figure
pyplot.figure()

## import empirical SFS
dd = dadi.Misc.make_data_dict_vcf("./0miss_vcfs/ALBACEN_str_full_06_thin.recode.vcf", "popfile_albacen.txt")

# create/compute JSFS
# N-Coastal - S-Coastal
fs = dadi.Spectrum.from_data_dict(dd, ["ALBA", "CEN"], projections=[12, 38], polarized=False)

params = [0.0124,2.3312,1.3078,8.4609,2.6401,0.0628,0.5251,0.4042] ## sec_con_asym_mig_size
pts = np.array([40, 50, 60])
ns = [12,38]

## Model for plotting

sec_con_asym_mig_size_func = Models_2D.sec_contact_asym_mig_size
sec_con_asym_mig_size_ex = dadi.Numerics.make_extrap_log_func(sec_con_asym_mig_size_func)
model = sec_con_asym_mig_size_ex(params, ns, pts)

pylab.figure(figsize=(8,6))
dadi.Plotting.plot_2d_comp_multinom(model, fs, vmin=1, resid_range=5,
                                    pop_ids =('S-Coastal','N-Coastal'), show=True)


## S-Oceanic - S-Coastal
dd = dadi.Misc.make_data_dict_vcf("./0miss_vcfs/TMALBA_str_full_06_thin.recode.vcf", "popfile_albatm.txt")
fs = dadi.Spectrum.from_data_dict(dd, ["ALBA", "TM"], projections=[10,16], polarized=False)

params = [14.0558,1.1865,1.0188,1.521,6.1548]## asym_mig
pts = np.array([40, 50, 60])
ns = [10,16]

## Model for plotting


asym_mig_func = Models_2D.asym_mig
asym_mig_func_ex = dadi.Numerics.make_extrap_log_func(asym_mig_func)
model = asym_mig_func_ex(params, ns, pts)

# asym_mig_size_func = Models_2D.asym_mig_size
# asym_mig_size_ex = dadi.Numerics.make_extrap_log_func(asym_mig_size_func)
# model = asym_mig_size_ex(params, ns, pts)

pylab.figure(figsize=(8,6))
dadi.Plotting.plot_2d_comp_multinom(model, fs, vmin=1, resid_range=5,
                                    pop_ids =('S-Coastal','S-Oceanic'), show=True)




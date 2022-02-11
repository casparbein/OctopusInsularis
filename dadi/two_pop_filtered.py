import dadi
import itertools
import matplotlib.pyplot as pyplot
from Two_Population_Pipeline import Models_2D
import numpy as np
import Optimize_Functions
from Find_Best_Projections import dadi_test_projections

## import data
## example for N-Coastal - S-Coastal

dd = dadi.Misc.make_data_dict_vcf("~/ALBACEN_str_full_06_thin.recode.vcf", "~/popfile_albacen.txt")

run_projections(dd_at, pop_ids=["pop1", "pop3"], polarized=False, maxproj=[26, 24], min_frac=0.25)
## best projection is [12,38]

# create/compute JSFS

fs = dadi.Spectrum.from_data_dict(dd, ["ALBA", "CEN"], projections=[12, 38], polarized=False)
#print(fs_at)

## add 2pop divergence models to test
my_base_func = dadi.Demographics2D.snm

## basic models
no_mig_func = Models_2D.no_mig
syn_mig_func = Models_2D.sym_mig
asym_mig_func = Models_2D.asym_mig

## basic models with instantaneous size change
no_mig_size_func = Models_2D.no_mig_size
sym_mig_size_func = Models_2D.sym_mig_size
asym_mig_size_func = Models_2D.asym_mig_size

# Ancient migration or secondary contact with instantaneous size change
sec_contact_sym_mig_size_func = Models_2D.sec_contact_sym_mig_size
sec_contact_asym_mig_size_func = Models_2D.sec_contact_asym_mig_size

## Island Models: Vicariance and continous migration
vic_no_mig_func = Models_2D.vic_no_mig
vic_anc_asym_mig_func = Models_2D.vic_anc_asym_mig
vic_sec_asym_mig_func = Models_2D.vic_sec_contact_asym_mig

## Island Models: Founder events and continuous migration
founder_nomig_func = Models_2D.founder_nomig
founder_sym_mig_func = Models_2D.founder_sym
founder_asym_mig_func = Models_2D.founder_asym

## model fitting
## setting up initial values, not done for all models except the neutral one and the grid size
neutral_params = np.array([])
pts = np.array([50, 60, 70])

## Optimize Function
pts = np.array([60, 70, 80])
reps = [15,10,5,5]
maxiters = [15,15,25,50]
folds = [3,2,1,1]

## write a loop that runs the command many times
upper = [1000000,20,0.99]
lower = [0.001,0.001,0.001]
for i in range(0,5):
    var_name = f"run_CENALBA_{i}"
    print(var_name)
    Optimize_Functions.Optimize_Routine(fs, pts, var_name, "no_mig_ca", no_mig_func, 4, 3, fs_folded=True,
                                        reps=reps, maxiters=maxiters, folds=folds)
    Optimize_Functions.Optimize_Routine(fs, pts, var_name, "sym_mig_ca", syn_mig_func, 4, 4, fs_folded=True,
                                        reps=reps, maxiters=maxiters, folds=folds)
    Optimize_Functions.Optimize_Routine(fs, pts, var_name, "asym_mig_ca", asym_mig_func, 4, 5, fs_folded=True,
                                        reps=reps, maxiters=maxiters, folds=folds)
    Optimize_Functions.Optimize_Routine(fs, pts, var_name, "no_mig_size_ca", no_mig_size_func, 4, 6, fs_folded=True,
                                        reps=reps, maxiters=maxiters, folds=folds)
    Optimize_Functions.Optimize_Routine(fs, pts, var_name, "sym_mig_size_ca", sym_mig_size_func, 4, 7, fs_folded=True,
                                        reps=reps, maxiters=maxiters, folds=folds)
    Optimize_Functions.Optimize_Routine(fs, pts, var_name, "asym_mig_size_ca", asym_mig_size_func, 4, 8, fs_folded=True,
                                        reps=reps, maxiters=maxiters, folds=folds)
    Optimize_Functions.Optimize_Routine(fs, pts, var_name, "vic_no_mig_ca",vic_no_mig_func , 4, 2, fs_folded=True,
                                        reps=reps, maxiters=maxiters, folds=folds, in_upper=[20,0.99], in_lower=[0.001,0.001])
    Optimize_Functions.Optimize_Routine(fs, pts, var_name, "vic_anc_asym_mig_ca", vic_anc_asym_mig_func, 4, 5, fs_folded=True,
                                        reps=reps, maxiters=maxiters, folds=folds, in_upper=[20,20,20,20,0.99], in_lower=[0.001,0.001,0.001,0.001,0.001])
    Optimize_Functions.Optimize_Routine(fs, pts, var_name, "vic_sec_asym_mig_ca", vic_sec_asym_mig_func, 4, 5, fs_folded=True,
                                        reps=reps, maxiters=maxiters, folds=folds, in_upper=[20,20,20,20,0.99], in_lower=[0.001,0.001,0.001,0.001,0.001])
    Optimize_Functions.Optimize_Routine(fs, pts, var_name, "founder_nomig_ca", founder_nomig_func, 4, 3, fs_folded=True,
                                        reps=reps, maxiters=maxiters, folds=folds, in_lower=lower, in_upper=upper)
    Optimize_Functions.Optimize_Routine(fs, pts, var_name, "founder_sym_ca", founder_sym_mig_func, 4, 4, fs_folded=True,
                                        reps=reps, maxiters=maxiters, folds=folds, in_upper=[1000000,20,20,0.99], in_lower=[0.001,0.001,0.001,0.001])
    Optimize_Functions.Optimize_Routine(fs, pts, var_name, "founder_asym_ca", founder_asym_mig_func, 4, 5, fs_folded=True,
                                        reps=reps, maxiters=maxiters, folds=folds, in_upper=[1000000,20,20,20,0.99], in_lower=[0.001,0.001,0.001,0.001,0.001])
    Optimize_Functions.Optimize_Routine(fs, pts, var_name, "sec_contact_sym_mig_size_ca", sec_contact_sym_mig_size_func, 4, 7, fs_folded=True,
                                        reps=reps, maxiters=maxiters, folds=folds)
    Optimize_Functions.Optimize_Routine(fs, pts, var_name, "sec_contact_asym_mig_size_ca", sec_contact_asym_mig_size_func, 4, 8, fs_folded=True,
                                        reps=reps, maxiters=maxiters, folds=folds)

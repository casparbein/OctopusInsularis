import dadi
import matplotlib.pyplot as pyplot
from Two_Population_Pipeline import Models_2D
import numpy as np
import Optimize_Functions
from Find_Best_Projections import dadi_test_projections
from Goodness_of_Fit import Optimize_Functions_GOF

## N-Coastal - S-Coastal
dd = dadi.Misc.make_data_dict_vcf("~/ALBACEN_str_full_06_thin.recode.vcf", "~/popfile_albacen.txt")

## S-Oceanic - S-Coastal
#dd = dadi.Misc.make_data_dict_vcf("~/TMALBA_str_full_06_thin.recode.vcf", "~/popfile_albatm.txt")

print(dd)
len(dd)

# create/compute JSFS
# N-Coastal - S-Coastal
fs = dadi.Spectrum.from_data_dict(dd, ["ALBA", "CEN"], projections=[12, 38], polarized=False)

# S-Oceanic - S-Coastal
#fs = dadi.Spectrum.from_data_dict(dd, ["ALBA", "TM"], projections=[10, 16], polarized=False)

## max projections
max_proj = [24,66]
#max_proj = [24,24]

# Convert this dictionary into folded AFS object based on
# MAXIMUM projection sizes
# [polarized = False] creates folded spectrum object
max_fs = dadi.Spectrum.from_data_dict(dd, pop_ids=["ALBA", "CEN"], projections = max_proj, polarized = False)
#max_fs = dadi.Spectrum.from_data_dict(dd, pop_ids=["ALBA", "TM"], projections = max_proj, polarized = False)

## grid size
pts = [80,90,100]

## model
## N-Coastal - S-Coastal
sec_con_asym_mig_size_func = Models_2D.sec_contact_asym_mig_size

## S-Oceanic - S-Coastal
#asym_mig_func = Models_2D.asym_mig

##best fitting parameters
emp_params = [0.0124,2.3312,1.3078,8.4609,2.6401,0.0628,0.5251,0.4042] ## sec_con_asym_mig_size, N-Coastal - S-Coastal 
#emp_params = [14.0558,1.1865,1.0188,1.521,6.1548] ## asym_mig, S-Oceanic - S-Coastal

## folded
fs_folded = True

## account for different projection sizes
fs_for_sims = Optimize_Functions_GOF.Get_Empirical(fs, max_fs, pts, "Empirical", "sec_con_asym_mig_size", sec_con_asym_mig_size_func, emp_params, [12,38], fs_folded=fs_folded)

## simulations

#Set the number of simulations to perform here. This should be ~100 or more.
sims = 100

#**************
#Enter the number of parameters found in the model to test.
pnum = 8
#pnum = 5

#**************
#Set the number of rounds here.
rounds = 4
reps = [15,5,5,5]
maxiters = [5,5,5,10]
folds = [3,2,1,1]

## labels
p_labels = "nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2"
#p_labels = "nu1, nu2, m12, m21, T"

Optimize_Functions_GOF.Perform_Sims(sims, fs_for_sims, pts, "sec_con_asym_mig_size", 
                                    sec_con_asym_mig_size_func, rounds, pnum, projections=[12,38], 
                                    fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels=p_labels)

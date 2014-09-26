## this is the script for running reactive vac similations for three
## populations
import scipy.integrate as spi
import numpy as np
import pylab as pl
import warnings
from matplotlib.font_manager import FontProperties as fmp
import sys
import pdb
import vacfunctions as vf
import re
import os
import time


start=time.time()

## this time we are going to expect a Hot Cold configuration given
## in the form spotconfig="HHC"
# check that spotconfig is of length 3

# if len(sys.argv) == 2:
#     spotconfig = re.split("=",sys.argv[1])[1]
#     if spotconfig[0] != "H" or spotconfig[1] != "C":
#         print "Warning: this only works for configs starting with HC at the moment"
# else:
#    print "No hotspot config specified, assuming Hot-Cold-Hot"
spotconfig = "HCC"

pvacs = np.linspace(0,499999.0/500000,41)
times = range(0,300,1)
Rhots = np.linspace(1.5,2.5,5)
Rnots = np.linspace(0.75,1.5,5)
Rs = vf.expandgrid(Rnots,Rhots)

Rs0 = Rs['Var2'] if spotconfig[0] == "H" else  Rs['Var1']
Rs1 = Rs['Var2'] if spotconfig[1] == "H" else  Rs['Var1']
Rs2 = Rs['Var2'] if spotconfig[2] == "H" else  Rs['Var1']

ptiles = np.linspace(0,1.0,41)
my_N = 500000
my_nlocs = 3
len_Rs = len(Rs["Var1"])
my_hmax = 0.5 # max step size for numerical integrator

ci_null = np.zeros(len_Rs)
fs = np.zeros([len_Rs,len(ptiles),len(pvacs),4])
null_run = [0]*len_Rs
timings_save = np.zeros(len_Rs*len(ptiles)).reshape(len_Rs,len(ptiles))
#threshold for which if there are no cases at a step we say the epimidc is over
vac_end_threshold = .05


my_alphas=np.array([0.0,0.01,0.2])

for a in range(len(my_alphas)):
    my_alpha = my_alphas[a]
    time_stamp = re.sub("\.","",str(time.time()))
    run_stamp = spotconfig + "_alpha" + re.sub("\.","p",str(my_alpha)[0:len(str(my_alpha))]) + "_" + time_stamp
    path = "Generated_Data/npy_" + run_stamp  + "/"
    os.makedirs(path)
    #sys.stdout = open("/Users/aazman/Documents/ReactiveVac/Run_Logs/log_" + run_stamp  + ".log",'w')
    #temp = sys.stdout # store orignla stout for later

    for r in range(len_Rs):
        print ("\n R0 = %.2f, R1 = %.2f, R2 = %.2f, scenario %s of %s" %
               (Rs0[r],Rs1[r],Rs2[r],r,len_Rs))
        ins = vf.make_inputs(my_N,my_nlocs,
                             np.array([Rs0[r],Rs1[r],Rs2[r]]),my_alpha)

        # do a null run of the model
        null_run[r] = spi.odeint(vf.SIR_dx_dt,
                                 ins[0],hmax=my_hmax,
                                 t=times,
                                 args=(ins[1],ins[2],ins[3],0,0,
                                       np.zeros(my_nlocs)))
        # get cumulative incidence
        ci = vf.get_ci_columns(null_run[r])
        #print ci

        cases_per_step = np.array([np.diff(ci[:,x]) for x in range(ci.shape[1])])
        cases_per_step = np.array([sum(locs) for locs in cases_per_step.T])

        # figure out when the epidemic ended
        epidemic_end = ((cases_per_step < vac_end_threshold).nonzero())[0]

        #if we didnt go long enough
        while len(epidemic_end) == 0:
            times = range(0,max(times)+30,1)
            null_run[r] = spi.odeint(vf.SIR_dx_dt,
                                     ins[0],
                                     t=times,
                                     args=(ins[1],ins[2],ins[3],0,0,
                                           np.zeros(my_nlocs)))
            # array with cis in all locations per column
            ci = vf.get_ci_columns(null_run[r])

            cases_per_step = np.array(
                [np.diff(ci[:,x]) for x in range(ci.shape[1])]
                )
            cases_per_step = np.array(
                [sum(locs) for locs in cases_per_step.T]
                )
            # figure out when the epidemic ended
            epidemic_end = (
                (cases_per_step < vac_end_threshold).nonzero()
                )[0]

            #        print "epidemic end day = %s" % (epidemic_end[0])

        # get the final size for the this uncontrolled epidemic
        ci_null[r] = max(sum(x) for x in ci[0:epidemic_end[0],:])
        global_ci = np.array([sum(x) for x in ci[0:epidemic_end[0],:]])
        timings_save[r,:] = np.interp(ptiles,global_ci/ci_null[r],
                                      range(epidemic_end[0]))

        # using exact time of percentiles. could consider casting as an
        # integer later
        timings_ref = [x for x in timings_save[r,:]]
        # loop over percent vaccinated
        for pv in range(len(pvacs)):
            pct_vac = pvacs[pv]
            total_vac = pct_vac*(my_N)
            # get NaN by dividing by zero in S/(N-V)
            #print "vaccinating %d pct of the total population" % (pct_vac*100)
            sys.stdout.write('|')

            ## vaccination distribution scenarios
            vac_daily_target_hs = np.concatenate(([total_vac],
                                                  [0]*(my_nlocs-1)))
            vac_daily_target_nhs = np.concatenate(([0],
                                                   [total_vac],[0]))
            vac_daily_diffuse = np.array([total_vac/my_nlocs]*my_nlocs)

            ## now for targeting both hotspots or both coldsptos
            if spotconfig[2] == "H": # share between loc 0 and 2
                vac_daily_target_double = np.concatenate(([total_vac/(my_nlocs-1)],
                                                          [0],
                                                          [total_vac/(my_nlocs-1)]))
            else:
                vac_daily_target_double = np.concatenate(([0],
                                                          [total_vac/(my_nlocs-1)],
                                                          [total_vac/(my_nlocs-1)]))

            for timing in range(len(timings_ref)):
                vac_start = timings_ref[timing]
                vac_end = timings_ref[timing]+1 # TODO: check that this +1 is required in python

                sys.stdout.write('.')

                # print "vac started at %s (%0.2f percentile)" % (vac_start,ptiles[timing])

                #print "targeting the hotspot"
                #pdb.set_trace()

                tmp_run = spi.odeint(vf.SIR_dx_dt,
                                     ins[0],
                                     t=times,hmax=my_hmax,
                                     args=(ins[1],ins[2],ins[3],
                                           vac_start,vac_end,
                                           vac_daily_target_hs))

                fs[r,timing,pv,0] = sum(vf.get_final_attack_size(tmp_run))
                # print fs[r,timing,pv,0]
                #print "targeting the coldspot"
                tmp_run = spi.odeint(vf.SIR_dx_dt,
                                     ins[0],
                                     t=times,hmax=my_hmax,
                                     args=(ins[1],ins[2],ins[3],
                                           vac_start,vac_end,
                                           vac_daily_target_nhs))

                fs[r,timing,pv,1] = sum(vf.get_final_attack_size(tmp_run))
                #print fs[r,timing,pv,1]

                #print "pro-rata vaccination"
                tmp_run = spi.odeint(vf.SIR_dx_dt,
                                     ins[0],
                                     t=times,hmax=my_hmax,
                                     args=(ins[1],ins[2],ins[3],
                                           vac_start,vac_end,
                                           vac_daily_diffuse))

                fs[r,timing,pv,2] = sum(vf.get_final_attack_size(tmp_run))
                #print fs[r,timing,pv,2]

                #print "double spot vaccination"

                tmp_run = spi.odeint(vf.SIR_dx_dt,
                                     ins[0],
                                     t=times,hmax=my_hmax,
                                     args=(ins[1],ins[2],ins[3],
                                           vac_start,vac_end,
                                           vac_daily_target_double))

                fs[r,timing,pv,3] = sum(vf.get_final_attack_size(tmp_run))
                #print fs[r,timing,pv,3]

                # we can only easily import 2-d arrays into R so we will go ahead and make this
        np.save(path + "fs_r" + str(r) + "_type1_" + run_stamp + ".npy",fs[r,:,:,0])
        np.save(path + "fs_r" + str(r) + "_type2_" + run_stamp + ".npy",fs[r,:,:,1])
        np.save(path + "fs_r" + str(r) + "_type3_" + run_stamp + ".npy",fs[r,:,:,2])
        np.save(path + "fs_r" + str(r) + "_type4_" + run_stamp + ".npy",fs[r,:,:,3])

        # saving the full output from every state in the model for each set of R
        np.save(path + "nullrun_r" + str(r) + "_" + run_stamp + ".npy",null_run[r])

        np.save(path + "timings_save_" + run_stamp + ".npy",timings_save)

# sys.stdout.close()
# sys.stdout = temp
elapsed = (time.time() - start)
print "Elapsed time is %s" % elapsed

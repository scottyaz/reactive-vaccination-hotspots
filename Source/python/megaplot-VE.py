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

pvacs = np.linspace(0.0,499999.0/500000,41) # note this is % of single patch not
# full population
times = range(0,300,1)
Rhots = np.linspace(1.5,2.5,5)
Rnots = np.linspace(0.75,1.5,5)
#Rnots = np.linspace(0.9,1.5,5)

my_theta = 0.6 # vaccine efficiacy
Rs = vf.expandgrid(Rnots,Rhots)
ptiles = np.linspace(0.0,1.0,41)
my_N = 500000
my_nlocs = 2
my_hmax= 0.5
len_Rs = len(Rs["Var1"])

Rs0 = Rs['Var2'] # hot
Rs1 = Rs['Var1'] # cold

ci_null = np.zeros(len_Rs)
fs = np.zeros([len_Rs,len(ptiles),len(pvacs),3])
null_run = [0]*len_Rs
timings_save = np.zeros(len_Rs*len(ptiles)).reshape(len_Rs,len(ptiles))
#threshold for which if there are no cases at a step we say the epimidc is over
vac_end_threshold = .05

alphas = np.array([0.0,0.01,0.2])

for a in range(len(alphas)):
    my_alpha = alphas[a]
    print "running simulations for alpha = %s" % my_alpha
    start=time.time()
    time_stamp = re.sub("\.","",str(time.time()))
    run_stamp = "alpha" + re.sub("\.","p",str(my_alpha)[0:len(str(my_alpha))]) + "_" + time_stamp

    # make directory for output
    path = "Generated_Data/npy_" + run_stamp  + "/"
    os.makedirs(path)
    #    sys.stdout = open("/Users/aazman/Documents/ReactiveVac/Run_Logs/log_" + run_stamp  + ".log",'w')
    #    temp = sys.stdout # store orignla stout for later

    for r in range(len(Rs0)):
        print ("\n R1 = %.2f, R2 = %.2f, scenario %s of %s" %
        (Rs0[r],Rs1[r],r,len_Rs))
        ## ins[0] - initial state
        ## ins[1] - contact mat
        ## ins[2] - Rs
        ## ins[3] - theta
        ## ins[4] - gamma
        ins = vf.make_inputs_VE(my_N,my_nlocs,
                                np.array([Rs0[r],Rs1[r]]),my_alpha,
                                my_theta)
        ins_orig = vf.make_inputs(my_N,my_nlocs,
                                np.array([Rs0[r],Rs1[r]]),my_alpha
                                )

        # do a null run of the model
        null_run[r] = spi.odeint(vf.SIR_VE_dx_dt,
                                 ins[0],hmax=0.9, # state = ins[0]
                                 t=times,
                                 args=(ins[1],ins[2],ins[4],ins[3],0,0,
                                               np.zeros(my_nlocs)))


        test = spi.odeint(vf.SIR_dx_dt,
                          ins_orig[0],hmax=0.9, # state = ins[0]
                          t=times,
                          args=(ins_orig[1],ins_orig[2],ins_orig[3],0,0,
                                np.zeros(my_nlocs)))

        # get cumulative incidence
        ci = vf.get_ci_columns(null_run[r],6)
        #print ci

        cases_per_step = np.array([np.diff(ci[:,x]) for x in range(ci.shape[1])])
        cases_per_step = np.array([sum(locs) for locs in cases_per_step.T])

        # figure out when the epidemic ended
        epidemic_end = ((cases_per_step < vac_end_threshold).nonzero())[0]

        #if we didnt go long enough
        while len(epidemic_end) == 0:
            times = range(0,max(times)+30,1)
            null_run[r] = spi.odeint(vf.SIR_VE_dx_dt,
                                     ins[0],
                                     t=times,
                                     hmax=my_hmax,args=(ins[1],ins[2],ins[4],
                                                        ins[3],0,0,
                                                        np.zeros(my_nlocs)))
            # array with cis in all locations per column
            ci = vf.get_ci_columns(null_run[r],6)
            cases_per_step = np.array([np.diff(ci[:,x]) for x in
                                       range(ci.shape[1])])
            cases_per_step = np.array(
            [sum(locs) for locs in cases_per_step.T]
            )
            # figure out when the epidemic ended
            epidemic_end = (
            (cases_per_step < vac_end_threshold).nonzero()
            )[0]

            # print "epidemic end day = %s" % (epidemic_end[0])

            # get the final size for the this uncontrolled epidemic
        ci_null[r] = max(sum(x) for x in ci[0:epidemic_end[0],:])
        global_ci = np.array([sum(x) for x in ci[0:epidemic_end[0],:]])
        # determine the vaccination timeings based on percentiles of uncon epidemic
        timings_save[r,:] = np.interp(ptiles,global_ci/ci_null[r],
                                      range(epidemic_end[0]))

        timings_ref = [x for x in timings_save[r,:]]
        # loop over percent vaccinated
        for pv in range(len(pvacs)):

            #print "pv: %s" %  pv
            pct_vac = pvacs[pv]
            #print "vaccinating %d pct of the total population" % (pct_vac*100)
            sys.stdout.write('|')

            vac_daily_diffuse = np.array([pct_vac*my_N/my_nlocs]*my_nlocs)
            vac_daily_target_hs = np.concatenate(([pct_vac*my_N],
                                                  [0]*(my_nlocs-1)))
            vac_daily_target_nhs = np.concatenate(([0],[pct_vac*my_N/(my_nlocs-1)]*
                                                   (my_nlocs-1)))


            for timing in range(len(timings_ref)):
                #print timing

                vac_start = timings_ref[timing]
                vac_end = timings_ref[timing]+1

                sys.stdout.write('.')
                #print "vac started at %s (%0.2f percentile)" % (vac_start,ptiles[timing])

                #print "targeting the hotspot"
                tmp_run1 = spi.odeint(vf.SIR_VE_dx_dt,
                                     ins[0],
                                     t=times,hmax=my_hmax,
                                     args=(ins[1],ins[2],ins[4],ins[3],
                                           vac_start,vac_end,
                                           vac_daily_target_hs))
                #print sum(vf.get_final_attack_size(tmp_run1))
                fs[r,timing,pv,0] = sum(vf.get_final_attack_size(tmp_run1,6))

                #print "targeting the coldspot"
                tmp_run2 =spi.odeint(vf.SIR_VE_dx_dt,
                                    ins[0],
                                    t=times,hmax=my_hmax,
                                    args=(ins[1],ins[2],ins[4],ins[3],
                                          vac_start,vac_end,
                                          vac_daily_target_nhs))

                #print sum(vf.get_final_attack_size(tmp_run2))
                fs[r,timing,pv,1] = sum(vf.get_final_attack_size(tmp_run2,6))

                #print "pro-rata vaccination"
                tmp_run3 = spi.odeint(vf.SIR_VE_dx_dt,
                                     ins[0],
                                     t=times,hmax=my_hmax,
                                     args=(ins[1],ins[2],ins[4],ins[3],
                                           vac_start,vac_end,
                                           vac_daily_diffuse))

                # print sum(vf.get_final_attack_size(tmp_run3))
                fs[r,timing,pv,2] = sum(vf.get_final_attack_size(tmp_run3,6))

    elapsed = (time.time() - start)
                # we can only easily import 2-d arrays into R so we will go ahead and make this
    for rs in range(len_Rs):
        np.save(path + "fs_r" + str(rs) + "_type1_VE" + str(int(my_theta*100)) + "_" + run_stamp + ".npy",fs[rs,:,:,0])
        np.save(path + "fs_r" + str(rs) + "_type2_VE" + str(int(my_theta*100)) + "_" + run_stamp + ".npy",fs[rs,:,:,1])
        np.save(path + "fs_r" + str(rs) + "_type3_VE" + str(int(my_theta*100)) + "_" + run_stamp + ".npy",fs[rs,:,:,2])
        # saving the full output from every state in the model for each set of Rs
        np.save(path + "nullrun_r" + str(rs) + "_VE" + str(int(my_theta*100)) + "_" + run_stamp + ".npy",null_run[rs])
        np.save(path + "timings_save_VE" + str(int(my_theta*100)) + "_" + run_stamp + ".npy",timings_save)
        #        sys.stdout.close()
        # sys.stdout = temp

        print "Elapsed time is %s" % elapsed

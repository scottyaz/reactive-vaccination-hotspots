###############################################################
## Script for finding the optimal allocation to a hotspot in ##
## a simple 2-patch metapopulation                           ##
###############################################################
import numpy as np
import scipy.integrate as spi
import scipy.optimize
import pdb
import time
import os
import sys
import warnings
import re
import vacfunctions as vf

##########################################
## some functions used for optimization ##
##########################################

def two_pop_finalsize(individual,**options):
    """evalutates final size for a vaccine scenario:
    input:
    1. individual = vaccines in hotspot (will impute the other loc)
    2. options = all the deets for the simulation
    returns:
    tuple of the total funal size for the controlled epidemic
    [NOTE: this is used for optimization purposes]
    """

    my_alpha = options.get("alpha")
    my_N = options.get("N")
    my_nlocs = options.get("nlocs")
    R0 = options.get("R0")
    R1 = options.get("R1")
    # R2 = options.get("R2")
    my_vac_start = options.get("vac_start")
    my_times = options.get("times")
    my_hmax = options.get("hmax")
    vac_max = options.get("vac_max")

    individual = np.array([individual,vac_max-individual])

    # pdb.set_trace()

    ins = vf.make_inputs(my_N,
                         my_nlocs,
                         np.array([R0,R1]),
                         my_alpha,1.0/3,[0.0])

    # make vaccine allocation vector
    vac_daily = np.array(individual)

    tmp_run = spi.odeint(vf.SIR_dx_dt,
                         ins[0],
                         t=my_times,
                         hmax=my_hmax,
                         args=(ins[1],
                               ins[2],
                               ins[3],
                               my_vac_start,
                               my_vac_start+1,
                               vac_daily))

    final_size = sum(vf.get_final_attack_size(tmp_run))[0]

    return final_size

def two_pop_finalsize_helper(individual,alpha,N,nlocs,
             R0,R1,
             vac_start,
             times,
             hmax,
             vac_max):

    """helper evaluate function with all the defaults loaded """

    rc = two_pop_finalsize(individual,
             alpha=alpha,
             N=N,
             nlocs=nlocs,
             R0=R0,
             R1=R1,
             vac_start=vac_start,
             times=times,
             hmax=hmax,
             vac_max=vac_max)

    return rc


####################################
## Some default parameters to set ##
####################################
pvacs = np.linspace(0.0,499999.0/500000,21)
times = range(0,300,1)
Rhots = np.linspace(1.5,2.5,5)
Rnots = np.linspace(0.75,1.5,5)

Rs = vf.expandgrid(Rnots,Rhots)
ptiles = np.linspace(0.0,1,21)
my_N = 500000
my_nlocs = 2
my_hmax= 0.5
len_Rs = len(Rs["Var1"])

Rs0 = Rs['Var2'] # hot
Rs1 = Rs['Var1'] #cold

ci_null = np.zeros(len_Rs)
opt = np.zeros([len_Rs,len(ptiles),len(pvacs)])
null_run = [0]*len_Rs
timings_save = np.zeros(len_Rs*len(ptiles)).reshape(len_Rs,len(ptiles))
#threshold for which if there are no cases at a step we say the epimidc is over
vac_end_threshold = .05

alphas = np.array([0.2])

##################################
## Start of the evaluation loop ##
##################################

for a in range(len(alphas)):
    my_alpha = alphas[a]
    print "running optimization simulations for alpha = %s" % my_alpha
    start=time.time()
    time_stamp = re.sub("\.","",str(time.time()))
    run_stamp = "opt_alpha" + re.sub("\.","p",str(my_alpha)[0:len(str(my_alpha))]) + "_" + time_stamp

    # make directory for output
    path = "Generated_Data/npy_" + run_stamp  + "/"
    os.makedirs(path)
    #    sys.stdout = open("/Users/aazman/Documents/ReactiveVac/Run_Logs/log_" + run_stamp  + ".log",'w')
    #    temp = sys.stdout # store orignla stout for later

    for r in range(len(Rs0)):
        print ("\n R1 = %.2f, R2 = %.2f, scenario %s of %s" %
        (Rs0[r],Rs1[r],r,len_Rs))
        ins = vf.make_inputs(my_N,my_nlocs,
                             np.array([Rs0[r],Rs1[r]]),my_alpha)

        # do a null run of the model
        null_run[r] = spi.odeint(vf.SIR_dx_dt,
                                 ins[0],hmax=0.9,
                                 t=times,args=(ins[1],ins[2],ins[3],0,0,
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
                             hmax=my_hmax,args=(ins[1],ins[2],ins[3],0,0,
                                           np.zeros(my_nlocs)))
    # array with cis in all locations per column
            ci = vf.get_ci_columns(null_run[r])
            cases_per_step = np.array([np.diff(ci[:,x]) for x in
                                       range(ci.shape[1])])
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

        timings_ref = [x for x in timings_save[r,:]]
        # loop over percent vaccinated
        for pv in range(len(pvacs)):

            #print "pv: %s" %  pv
            pct_vac = pvacs[pv]
            num_vac = pct_vac*my_N
            #print "vaccinating %d pct of the total population" % (pct_vac*100)
            sys.stdout.write('|')

            for timing in range(len(timings_ref)):
                #print timing

                vac_start = timings_ref[timing]
                vac_end = timings_ref[timing]+1

                sys.stdout.write('.')

                #print "vac started at %s (%0.2f percentile)" % (vac_start,ptiles[timing])

                tmp = scipy.optimize.minimize_scalar(fun=two_pop_finalsize_helper,
                                               args=(my_alpha,my_N,
                                                     my_nlocs,Rs0[r],Rs1[r],
                                                     vac_start,times,
                                                     my_hmax,num_vac),
                                                     method="bounded",
                                                     bounds=(0,num_vac))

                #                print tmp


                opt[r,timing,pv] = tmp['x']


    elapsed = (time.time() - start)
                # we can only easily import 2-d arrays into R so we will go ahead and make this
    for rs in range(len_Rs):
        np.save(path + "opt_r" + str(rs) + "_type1_" + run_stamp + ".npy",opt[rs,:,:])
        # saving the full output from every state in the model for each set of Rs
        np.save(path + "nullrun_r" + str(rs) + "_" + run_stamp + ".npy",null_run[rs])
        np.save(path + "timings_save_" + run_stamp + ".npy",timings_save)

        print "Elapsed time is %s" % elapsed

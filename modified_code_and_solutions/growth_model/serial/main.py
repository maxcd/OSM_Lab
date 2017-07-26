#======================================================================
#
#     This routine solves an infinite horizon growth model
#     with dynamic programming and sparse grids
#
#     The model is described in Scheidegger & Bilionis (2017)
#     https://papers.ssrn.com/sol3/papers.cfm?abstract_id=2927400
#
#     external libraries needed:
#     - IPOPT (https://projects.coin-or.org/Ipopt)
#     - PYIPOPT (https://github.com/xuy/pyipopt)
#     - TASMANIAN (http://tasmanian.ornl.gov/)
#
#     Simon Scheidegger, 11/16 ; 07/17
#======================================================================

import nonlinear_solver_initial as solver     #solves opt. problems for terminal VF
import nonlinear_solver_iterate as solviter   #solves opt. problems during VFI
from parameters import *                      #parameters of model
import interpolation as interpol              #interface to sparse grid library/terminal VF
import interpolation_iter as interpol_iter    #interface to sparse grid library/iteration
import postprocessing as post                 #computes the L2 and Linfinity error of the model

import TasmanianSG                            #sparse grid library
import numpy as np


#======================================================================
# Start with Value Function Iteration
# terminal value function
def run_all(n_agents):
    valnew=TasmanianSG.TasmanianSparseGrid()

    if (numstart==0):
        gridlist = []
        for tT in range(ntheta):
            valnew=TasmanianSG.TasmanianSparseGrid()
            valnew=interpol.sparse_grid(n_agents, iDepth, theta[tT])
            gridlist.append(valnew)
            valnew.write("valnew_1." + str(numstart) + "theta_" + str(tT) +  ".txt") #write file to disk for restart

    # value function during iteration
    else:
        gridlist = []
        for tT in range(ntheta):
            valnew.read("valnew_1." + str(numstart) + "theta_" + str(tT) +  ".txt")  #write file to disk for restart
            gridlist.append(valnew)

    #valold=TasmanianSG.TasmanianSparseGrid()
    #valold=valnew

    # avals_list = []

    for i in range(numstart, numits):
        print " ================================================= "
        print "             Iteration", i
        print " ================================================= " 
        for tT in range(ntheta):
            thet = theta[tT]
            valnew = TasmanianSG.TasmanianSparseGrid()
            valnew = interpol_iter.sparse_grid_iter(n_agents, iDepth, gridlist, thet)

            valnew.write("valnew_1." + str(i+1) + "theta_" + str(tT) + ".txt")
            gridlist[tT].copyGrid(valnew)
        
        # valnew0=TasmanianSG.TasmanianSparseGrid()
        # valnew1=TasmanianSG.TasmanianSparseGrid()
        # valnew2=TasmanianSG.TasmanianSparseGrid()
        # valnew3=TasmanianSG.TasmanianSparseGrid()
        # valnew4=TasmanianSG.TasmanianSparseGrid()
        #
        # valnew0=interpol_iter.sparse_grid_iter(n_agents, iDepth, valold, theta[0])
        # valnew1=interpol_iter.sparse_grid_iter(n_agents, iDepth, valold, theta[1])
        # valnew2=interpol_iter.sparse_grid_iter(n_agents, iDepth, valold, theta[2])
        # valnew3=interpol_iter.sparse_grid_iter(n_agents, iDepth, valold, theta[3])
        # valnew4=interpol_iter.sparse_grid_iter(n_agents, iDepth, valold, theta[4])

        # evaluate all grids at the same points
        # chosen arbitrarily to be the points of the third grid where theta=1
        # eval_points = valnew2.getPoints()
        #
        # aVals0 = valnew0.evaluateBatch(eval_points)[:,0]
        # aVals1 = valnew1.evaluateBatch(eval_points)[:,0]
        # aVals2 = valnew2.evaluateBatch(eval_points)[:,0]
        # aVals3 = valnew3.evaluateBatch(eval_points)[:,0]
        # aVals4 = valnew4.evaluateBatch(eval_points)[:,0]

        # print aVals0, aVals4

        # aVals_new = 0.2 *(aVals0 + aVals1 + aVals2 + aVals3 + aVals4)
        # aVals_new = np.reshape(aVals_new, (eval_points.shape[0], 1))

        #f=open("aVals_new.txt", 'a')
        #np.savetxt(f, aVals_new, fmt='% 2.16f')
        #f.close()

        # print "==================================================================="
        # print " print avals shape here :"
        # print aVals_new.shape
        #
        # valold=TasmanianSG.TasmanianSparseGrid()
        # valold.copyGrid(valnew2)
        # valold.loadNeededPoints(aVals_new)
        #valold=valnew

        # valold.write("valnew_1." + str(i+1) + "theta_" + str(tT) + ".txt")

    #======================================================================
    print "==============================================================="
    print " "
    print " Computation of a growth model of dimension ", n_agents ," finished after ", numits, " steps"
    print " "
    print "==============================================================="
    #======================================================================

    # compute errors
    post.ls_error(n_agents, numstart, numits, No_samples)

    #======================================================================
    print "==============================================================="
    print " "
    print " Errors are computed -- see errors_theta_NUM.txt"
    print " "
    print " Groupwork of Max, Clint and Ben at OSM Lab 2017"
    print "==============================================================="
    #======================================================================

run_all(n_agents)

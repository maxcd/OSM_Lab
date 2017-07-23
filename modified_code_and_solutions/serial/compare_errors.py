from main import run_all
import numpy as np

dims = np.array([1, 2, 3], dtype=int)
final_errs = np.zeros((len(dims), 2))

for i in range(len(dims)):
    agents = dims[i]
    print "\n\n=========================================="
    print "   Solving Model with", agents, " agents. "
    print "=========================================="  
    avg, max_err = run_all(agents)
    final_errs[i] = [avg, max_err]
       
print "Comparing errors for different numbers of dimensions and  refinement 2:"
print "dimensioins    average error     maximum error"  
for i in range(len(final_errs)):
    print(i, "\t", final_errs[i,0], "\t",final_errs[i,1]) 
    

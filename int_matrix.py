import numpy as np
import time
#import sys
#sys.path.append('C:/Users/tevid/Desktop/DNA_Origami/Code/InteractionMatrix')
#import pyximport; pyximport.install()
#import cy_int_matrix_func as imf
import int_matrix_func as imf
#import matplotlib.pyplot as plt
#import multiprocessing as mp
#import sys
import gc
#import os, psutil

#for each group of three columns we can only choose three total elements. If any of those choices are not in location (i,i)
#then it restricts one of the choices in the next group of three columns. Maybe the easiest way to do this is to brute force
#the calculation.

start = time.time()

#how many species do we want to have
k=4
save_src = './k'+str(k)+'-valid-tiling/'
save_csv = 'k'+str(k)+'-intmatrices.csv'
save_to_csv = False
save_img = False
element_run = False

if element_run:
    interaction_matricies = imf.generate_all_matrices_elements(k)
else:
    interaction_matricies = imf.generate_all_matrices(k)

#we should now be done with the loop
print('all* int matrices', len(interaction_matricies))

#process = psutil.Process(os.getpid())
#print('memory used', process.memory_info().rss)

valid_loops = []

if element_run:
    for e_M in interaction_matricies:
        #print(np.array(e_M).T)
        #M = imf.generate_m_from_elements(e_M, k)
        M = imf.generate_m_from_indices(e_M,k)
        #print(M)
        if imf.valid_coloring(M,k): valid_loops.append(M)
        
else:
    for M in interaction_matricies:
        M = np.unpackbits(M, axis=1, count=3*k)
        if imf.valid_coloring(M,k): valid_loops.append(M)

print('removed bad loops', len(valid_loops))

non_degen_list = imf.remove_degenerate_matrices(valid_loops, k)
      
print('removed degenerate matrices', len(non_degen_list))

final_list = non_degen_list

for i in range(len(final_list)):
    imf.draw_int_matrix(final_list[i], save_src, i,save=save_img, show=True)

    fname = save_src+'k'+str(k)+'-intmatrix-'+str(i)+'.csv'
    #print(fname)
    if save_to_csv:
        np.savetxt(fname, final_list[i].astype(int),delimiter=',')
        
print(time.time() - start)
        
#run_times = np.array([0.122,0.376,1.506,15.6,847.07,15601])

gc.collect()
        
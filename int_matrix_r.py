import numpy as np
import time
#import sys
#sys.path.append('C:/Users/tevid/Desktop/DNA_Origami/Code/InteractionMatrix')
#import pyximport; pyximport.install()
#import cy_int_matrix_func as imf
import int_matrix_r_func as imf
import matplotlib.pyplot as plt
#import multiprocessing as mp
#import sys
import gc
#import os, psutil

#for each group of three columns we can only choose three total elements. If any of those choices are not in location (i,i)
#then it restricts one of the choices in the next group of three columns. Maybe the easiest way to do this is to brute force
#the calculation.

start = time.time()

#how many species do we want to have
k=10
save_src = './k'+str(k)+'-valid-tiling/'
save_csv = 'k'+str(k)+'-intmatrices.csv'
save_to_csv = False
save_to_img = False
save_params = [save_src, save_csv, save_to_csv, save_to_img]

interaction_matrices = imf.generate_all_recursion(k,save_params)

print('generated matrices:', len(interaction_matrices))

np.save('valid_mat_list_k'+str(k)+'.npy', interaction_matrices)

'''
#now need to see which of the elements in this list are degenerate
final_list = imf.remove_degenerate_matrices(interaction_matrices, k)
print('final list length:', len(final_list))

#final_list = interaction_matrices

for i in range(len(final_list)):
    imf.draw_int_matrix(final_list[i], save_src, i,save=save_to_img, show=True)

    fname = save_src+'k'+str(k)+'-intmatrix-'+str(i)+'.csv'
    #print(fname)
    if save_to_csv:
        np.savetxt(fname, final_list[i].astype(int),delimiter=',')

np.save('final_mat_list.npy', final_list)
'''
print(time.time() - start)
        
#plt.show()

times_k_run = [0.2075, 0.258137, 1.8657, 15.36, 571.9]
times_k_run = [0.0209, 0.0827, 0.9295, 9.1839, 88.66, 829, 14527, 36457]
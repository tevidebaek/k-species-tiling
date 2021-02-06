import numpy as np
import time
#import sys
#sys.path.append('C:/Users/tevid/Desktop/DNA_Origami/Code/InteractionMatrix')
#import pyximport; pyximport.install()
#import cy_int_matrix_func as imf
import int_matrix_func as imf
import matplotlib.pyplot as plt
#import multiprocessing as mp

#for each group of three columns we can only choose three total elements. If any of those choices are not in location (i,i)
#then it restricts one of the choices in the next group of three columns. Maybe the easiest way to do this is to brute force
#the calculation.

start = time.time()

#how many species do we want to have
k=5
save_src = './k'+str(k)+'-valid-tiling/'
save_csv = 'k'+str(k)+'-intmatrices.csv'
save_to_csv = False
save_img = False

interaction_matricies = imf.generate_all_matrices(k)

#we should now be done with the loop
print('all* int matrices', len(interaction_matricies))
    
#decimated_matricies = imf.remove_bad_line_interactions(interaction_matricies, k)

#print('removed bad lines', len(decimated_matricies))

#decimated_matricies = interaction_matricies

valid_loops = []

#color_pool = mp.Pool()

#def check_color(M):
#    return imf.valid_coloring(M,k)

#color_bool = color_pool.map(check_color, interaction_matricies)

#valid_loops = interaction_matricies[color_bool]
#for M in decimated_matricies:
for M in interaction_matricies:    
    if imf.valid_coloring(M,k): valid_loops.append(M)

print('removed bad loops', len(valid_loops))

non_degen_list = imf.remove_degenerate_matrices(valid_loops, k)
#non_degen_list = imf.remove_degenerate_matrices(decimated_matricies, k)
      
print('removed degenerate matrices', len(non_degen_list))

#final_list = []

#for M in non_degen_list:
    
#    if imf.valid_coloring(M,k): final_list.append(M)

#print('removed bad loops', len(final_list))

final_list = non_degen_list

for i in range(len(final_list)):
    imf.draw_int_matrix(final_list[i], save_src, i,save=save_img, show=True)
    fname = save_src+'k'+str(k)+'-intmatrix-'+str(i)+'.csv'
    print(fname)
    if save_to_csv:
        np.savetxt(fname, final_list[i].astype(int),delimiter=',')
        
print(time.time() - start)
        
run_times = np.array([0.122,0.376,1.506,15.6,847.07,15601])

#plt.figure()
#plt.scatter([2.,3.,4.,5.,6.,7.], run_times, c='k')
#plt.semilogy()
#plt.ylabel('run time [s]', fontsize=16)
#plt.xlabel('k-species', fontsize=16)
#plt.show()

        
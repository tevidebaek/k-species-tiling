import numpy as np
import itertools
import matplotlib.pyplot as plt
import time
import os, psutil

#want to make drawings of the interaction matricies
def draw_int_matrix(M, img_dir='./', num_m=0, save=False, show=False):
    #print(M)
    k = int(len(M[0])/3)
    #print(k)
    
    fig, ax = plt.subplots(figsize=(k,k))
    
    #draw a dark line at each of the column-set/row-set boundaries
    for i in range(k+1):
        ax.plot([0,3*k],[3*i,3*i],'k', lw=2)
        ax.plot([3*i,3*i],[0,3*k],'k', lw=2)
        
        if i<k:
            #draw thin lines
            ax.plot([0,3*k],[3*i+1, 3*i+1],'k', lw=0.5)
            ax.plot([0,3*k],[3*i+2, 3*i+2],'k', lw=0.5)
            
            ax.plot([3*i+1, 3*i+1],[0,3*k],'k', lw=0.5)
            ax.plot([3*i+2, 3*i+2],[0,3*k],'k', lw=0.5)
    
    loc_elements = np.where(M>0)

    for i in range(np.shape(loc_elements)[1]):
        ax.scatter(loc_elements[0][i]+.5,3*k-loc_elements[1][i]-.5, c='w', edgecolor='k', s=80)
    
    ax.axis('off')

    if save==True:
        fig.savefig(img_dir+'intmatrix_k'+str(k)+'_'+str(num_m)+'.png',dpi=100)
    
    if show==False: plt.close(fig)
    
def swap_colors(M,a,b):
    #given a matrix we want to swap the column-set and row-set of color a with color b
    
    #make a copy of the matrix
    M_swap = np.copy(M)
    
    for i in range(3):
        M_swap[:,[3*a+i, 3*b+i]] = M_swap[:,[3*b+i,3*a+i]]
        M_swap[[3*a+i, 3*b+i],:] = M_swap[[3*b+i,3*a+i],:]
    
    return M_swap

def shift_index(M,shift):
    #shift can be 0,1,2
    if shift<0: return
    if shift>2: return
    
    #make a copy of the matrix
    M_swap = np.copy(M)
    
    #find k of matrix
    k = int(len(M_swap[0])/3)
    
    #print(k)
    
    if shift==1:
        for i in range(k):   #loops through each column-set and row-set
            M_swap[:,[3*i, 3*i + 1]] = M_swap[:, [3*i + 1, 3*i]]
            M_swap[[3*i, 3*i + 1],:] = M_swap[[3*i + 1, 3*i],:]

            M_swap[:,[3*i, 3*i + 2]] = M_swap[:, [3*i + 2, 3*i]]
            M_swap[[3*i, 3*i + 2],:] = M_swap[[3*i + 2, 3*i],:]
    elif shift==2:
        for i in range(k):   #loops through each column-set and row-set
            M_swap[:,[3*i, 3*i + 2]] = M_swap[:, [3*i + 2, 3*i]]
            M_swap[[3*i, 3*i + 2],:] = M_swap[[3*i + 2, 3*i],:]

            M_swap[:,[3*i, 3*i + 1]] = M_swap[:, [3*i + 1, 3*i]]
            M_swap[[3*i, 3*i + 1],:] = M_swap[[3*i + 1, 3*i],:]
            
            
    return M_swap

def get_color_permutations(M, k):
    M_ps = [M]

    for cnt in range(k-1):
        new_Mps = []
        
        for M_p in M_ps:
            for i in range(cnt,k):
                #mps = M_ps.copy()
                mp = swap_colors(M_p, cnt, i)
                new_Mps.append(mp)
                
        M_ps = new_Mps.copy()
        
    return M_ps

def check_same(M1,M2):
    if (M1==M2).all()==True: return True
    else: return False
    
def get_element_info(M):
    
    element_locations = np.where(M==1)
    
    element_info = []
    
    for i,j in zip(element_locations[0], element_locations[1]):
        #print((i,j))

        #given a pair of points find c-r, r-r
        c_s = np.floor(i/3.)
        r_s = np.floor(j/3.)
        side = (i-3*c_s)%3

        element_info.append([c_s, r_s, side])


    return np.array(element_info)
    
    
def find_loop(element_info, element_indx):
    #need to keep track of the starting element
    start_element = element_info[element_indx]
    
    next_element = np.array([-1, start_element[0], (start_element[2]-1)%3])
    
    current_element = next_element
    
    loop_elements = []
    #loop_elements.append(start_element)
    
    while not (current_element==start_element).all():
        #print(next_element)
        loc_next = np.where((element_info.T[1]==next_element[1])&((element_info.T[2]==next_element[2])))[0][0]
        current_element = element_info[loc_next]
        next_element = np.array([-1, current_element[0], (current_element[2]-1)%3])
        loop_elements.append(loc_next)
        
    return loop_elements
    
    
def find_all_loops(M):
    #M = non_degen_list[1]

    #draw_int_matrix(M,show=True)

    element_info = get_element_info(M)

    loops = []

    for i in range(len(element_info)):
        loop = find_loop(element_info,i)
        loops.append(np.sort(loop))
        if len(loop) not in [3,6]: return False, False

    unique_loops = []
    #go over the loops and keep only the unique ones

    unique_loops.append(loops[0])
    for l in loops[1:]:
        #check each following loop to see if it is already in the list
        new = True
        for ul in unique_loops:
            if len(ul)==len(l):
                if check_same(l,ul): new=False

        if new: unique_loops.append(l)

    return unique_loops, element_info

def get_loop_props(loop, element_info):
    #need to get both the colors involved in the loop as well as the length
    loop_len = len(loop)
    
    colors = []
    for element in loop:
        colors.append(element_info[element][0])
        colors.append(element_info[element][1])
    colors = np.unique(colors)
    
    return loop_len, colors

def valid_coloring(M, k):

    unique_loops, element_info = find_all_loops(M)
    
    if unique_loops==False: return False
    
    loop_props = []

    for ul in unique_loops:
        loop_len, colors = get_loop_props(ul, element_info)
        #print(loop_len, colors)

        if loop_len not in [3,6]: 
            print('bad loop length', loop_len)
            return False
        loop_props.append([loop_len, colors])

    #an important thing to keep in mind is that each loop corresponds to the coloring of a vertex within the tiling.
    #to have a valid tiling verticies must be allowed to be placed next to each other. What this means is that if there
    #is ever a loop in the interaction matrix that has a color that does not appear in any other loop, then that vertex 
    #cannot be placed next to any other and thus does not participate in the tiling. This means that the interaction
    #matrix supports two distinict tiling patterns that do not use all k colors. We will throw these out.

    #make a list for each color of the lists it is a part of
    start_loop = loop_props[0]

    overlapping_colors = start_loop[1]

    if len(loop_props)>1:
        for loop in loop_props[1:]:
            #check if any colors in loop match overlapping_colors
            if len(np.intersect1d(loop[1],overlapping_colors))!=0:
                overlapping_colors = np.union1d(overlapping_colors,loop[1])

    #check to see if all the colors are there
    if len(overlapping_colors)==k: connected_int_matrix = True
    else: connected_int_matrix = False

    return connected_int_matrix

def remove_degenerate_matrices(interaction_matricies,k):

    non_degen_list = []
    degen_list = interaction_matricies.copy()
    
    degen_elements = []
    for degen in degen_list:
        degen_elements.append(get_element_info(degen))
        
    while len(degen_list)>0:
        M = degen_list[0]
        
        degenerate_M = []
        
        def add_degen_M(M_add):
            if len(degenerate_M)!=0:
                for dM in degenerate_M:
                    if check_same(M_add,dM): 
                        return None
                degenerate_M.append(M_add)
            else: 
                degenerate_M.append(M_add)
                
        #we need to know how many possible color swaps there are (order does not matter)
        
        M_shift_1 = shift_index(M,1)
        M_shift_2 = shift_index(M,2)
        add_degen_M(M_shift_1)
        add_degen_M(M_shift_2)
        
        M_perms = get_color_permutations(M,k)
        
        for M_p in M_perms:
            M_shift_1 = shift_index(M_p,1)
            M_shift_2 = shift_index(M_p,2)
            
            add_degen_M(M_p)
            add_degen_M(M_shift_1)
            add_degen_M(M_shift_2)
    
        #should now have list of all degenerate matrices of M
        #check to see if any of them are already in non-degen-list
        
        new=True
        
        #if len(np.intersect1d(degenerate_M, non_degen_list))>0: new=False
        for dM in degenerate_M:
            for non_degen in non_degen_list:
                if check_same(dM, non_degen): new=False
                    
        if new:
            non_degen_list.append(M)

        #at this point try to do a more reasonable removal of the degenerate matrices from the original list
        idx_to_remove = []
    
        degen_copy = degen_list.copy()
        
        #degen_list = np.setdiff1d(degen_copy, degenerate_M)
        
        for dM in degenerate_M:
            for i in range(len(degen_list)):
                if check_same(dM, degen_list[i]):
                    idx_to_remove.append(i)      
                    continue

        for i in sorted(idx_to_remove, key=int, reverse=True):
            del(degen_copy[i])
        
        degen_list = degen_copy
    
        print(len(degen_list))
        
    print('non-degen length', len(non_degen_list))
    
    return non_degen_list

def remove_bad_line_interactions(M_list, k):
    #here we will remove elements that have lines forming along different lattice vectors.
    #can do this by checking each (col_set,row_set) if they have two elements
    
    valid_M = []
    
    #loop through all (row_set,col_set) and make list of which ones are two elements and which one is missing on the diagonal
    for M in M_list:

        M = np.unpackbits(M, axis=1, count=3*k)

        if check_bad_lines:
            valid_M.append(M)
            
    return valid_M
            
def check_bad_lines(M, k):
    
    set_props = []
    
    for i in range(k):
        for j in range(i,k):
            set_indx = ([0+3*i,1+3*i,2+3*i],[0+3*j,1+3*j,2+3*j])
            set_diag = M[set_indx]
            #print(set_diag)
            if np.sum(set_diag)==2:
                set_props.append(np.where(set_diag==0)[0][0])
            #if np.sum(set_diag)==3:
            #    if k>2: return False #removes two-line tilings from list
                
    if len(np.unique(set_props))<2:  
        return True
    
    else: 
        #print(set_props)
        #print(M)
        return False

'''
Based off of comments on the memory structure of trees it is more efficient to search a tree from a node all the way down to its root. Then we
will have a full matrix and can decide if we want to keep it or not. Then we can go up one level and follow the next branch down to the root.

To do this the general flow should be something like:
Generate the top matrix and keep track of the col_indx we are on.
Generate the list of combinations that we can do at this level.
    Loop through the combinations:
    For the first one genearte a new matrix and increment the col_indx
        If col_indx == k: check quality of matrix
        Else: generate combs for next level and check the elements
'''

def get_color_permutations_less(M, k):
    #this breaks down the permutations into sets of four swaps since we know that we will never fill some rows
    M_ps = [M]

    #how many iterations of swaps do we need to do
    num_full_swap = np.floor(k/4)

    for set_num in range(num_full_swap):
        
        for cnt in range(3):
            new_Mps = []
            
            for M_p in M_ps:
                for i in range(cnt,4):
                    #mps = M_ps.copy()
                    mp = swap_colors(M_p, cnt+set_num, i+set_num)
                    new_Mps.append(mp)
                    
            M_ps = new_Mps.copy()
        
    return M_ps


def check_degenerate(M, degen_list):
    
        ele_to_pop = []
        for i in range(len(degen_list)):
            if check_same(M, degen_list[i]): 
                ele_to_pop.append(i)
                #degen_list.pop(i)
                #return False
        if len(ele_to_pop)>0:
            for i in ele_to_pop[::-1]:
                degen_list.pop(i)
            return False

        return True

def add_degenerates(M,k, degen_list):
    #given the new matrix M, get all the degenerate copies of it
        
    temp_list = []

    def add_degen_M(M_add):
        if len(temp_list)!=0:
            #for dM in temp_list:
            #    if check_same(M_add,dM): 
            #        return None
            temp_list.append(M_add)
        else: 
            temp_list.append(M_add)
            
    #we need to know how many possible color swaps there are (order does not matter)
    
    #add_degen_M(M)

    M_shift_1 = shift_index(M,1)
    M_shift_2 = shift_index(M,2)
    #degen_list.append(M_shift_1)
    #degen_list.append(M_shift_2)
    add_degen_M(M_shift_1)
    add_degen_M(M_shift_2)
    
    M_perms = get_color_permutations(M,k)
    
    cnt=0

    for M_p in M_perms:
        M_shift_1 = shift_index(M_p,1)
        M_shift_2 = shift_index(M_p,2)
        
        #degen_list.append(M_p)
        #degen_list.append(M_shift_1)
        #degen_list.append(M_shift_2)

        add_degen_M(M_p)
        add_degen_M(M_shift_1)
        add_degen_M(M_shift_2)


        if cnt%100==0: print('degen', cnt, end='\r')
        cnt+=1

    for t in temp_list: degen_list.append(t)
    #degen_list = list(np.unique(degen_list))

def recursive_branch(k, col_set_indx,lowest_row_set, M, valid_list, degenerate_list, node_struc, save_params):
    '''
    This takes a generated matrix and sees if it is done or needs to look at the next level of the  tree

    After getting to the final level of the tree we need to store matrix if it passes our checks
    '''

    #things to print, what is the depth, and what is the 

    #print the current tree strucutre we are searching
    if col_set_indx<k/2.+1:
        process = psutil.Process(os.getpid())
        print("Tree strucutre", node_struc[:int(k/2)], 'mem (MB)', process.memory_info().rss/1000000.,'v:', len(valid_list), 'd:', len(degenerate_list) , end='\r')

    #print(col_set_indx, end='\r', flush=True)

    if col_set_indx==k:
        #check the quality of the matrix
        if valid_coloring(M,k):
            #we also need to see if this Matrix is already degenerate with another one already added
            #if check_degenerate(M, degenerate_list):
            i=len(valid_list)
            valid_list.append(M)
            #add_degenerates(M, k, degenerate_list)
            #save_src, save_csv, save_to_csv, save_to_img = save_params
            #if save_to_img: draw_int_matrix(M, save_src, i,save=save_img, show=True)
            #if save_to_csv: 
            #    fname = save_src+'k'+str(k)+'-intmatrix-'+str(i)+'.csv'
            #    np.savetxt(fname, M.astype(int),delimiter=',')

    else:
        #available columns in set
        col_indx = np.array([0,1,2])+3*col_set_indx

        #we want to see which columns in the column set of M do not have elements already.
        col_to_fill = []
        for col in col_indx: #looping through available columns
            if any(M[col])>0: continue
            else: col_to_fill.append(col)
                
        #if all columns of set are full then readd the matrix and go to next one, no additional elements to add
        if len(col_to_fill)==0:
            recursive_branch(k,col_set_indx+1, lowest_row_set,M, valid_list, degenerate_list, node_struc, save_params)
            return
            
        #if there are elements to add, do that next
        #for the first column we will never need to look past the fourth row since we can always
        #relabel things. This will also have consequences for the permutations we do
        if col_set_indx==0: combs = list(itertools.product(range(min(k,4)), repeat=len(col_to_fill)))    
        #for each column we have k choices, and in the set these choices do not interfer with each other
        #how many rows down do we need to fill? If nothing in the column, then 4 below the lowest row set filled
        #if there are elements in the col set, then fill len(col_to_fill)+1 below lowest row set
        else: combs = list(itertools.product(range(col_set_indx,min(k,len(col_to_fill)+1+lowest_row_set)), repeat=len(col_to_fill)))  #get all combinations of [0,...,k-1]
        #print(combs)
    
        cnt = 0

        #for each combination we need to create a new matrix
        for i in combs:
            lowest_row_set = max(i)

            cnt+=1  
            node_struc[col_set_indx] = cnt/len(combs)

            ######## add the extra elements to the matrix #########
            new_M = np.copy(M)

            #fill as many rows as the length of col_to_fill
            row_to_fill = []
            
            for j in range(len(col_to_fill)):
                #need to choose the correct row to fill. Recall that we will only be filling the diagonal elements of the row sets
                row = col_to_fill[j]%3 + 3*i[j]
                
                #check if row is already filled
                if any(new_M[row])>0: continue
                else: row_to_fill.append(row)
            
            #if no rows to fill from this combination, skip it
            if len(row_to_fill)!=len(col_to_fill): continue
            
            #print('row to fill', row_to_fill)
            
            for col, row in zip(col_to_fill, row_to_fill):
                new_M[col, row]=1
                new_M[row, col]=1

            ######## end of adding elements #########

            if np.sum(new_M[col_indx, col_indx])<3:
                if check_bad_lines(new_M,k):
                    #at this point the generated matrix passes our checks, go to the next branch
                    recursive_branch(k,col_set_indx+1, lowest_row_set,new_M, valid_list, degenerate_list, node_struc, save_params)

def generate_all_recursion(k, save_params):

    #start with a list of all the matrices we will eventually make
    interaction_matricies = []
    degenerate_list = []

    M = np.zeros((3*k, 3*k)).astype(int)  #initialize our first matrix
    ### end initilization of variables ###
        
    #need a counter for the column set that we are looking at. A column set are the group of three columns that correspond to
    #the interactions of a single particle species
    col_set_indx = 0   
    lowest_row_set = -1

    node_struc = np.zeros(k)

    #active the recursive look up
    recursive_branch(k,col_set_indx, lowest_row_set, M, interaction_matricies, degenerate_list, node_struc, save_params)

    return interaction_matricies
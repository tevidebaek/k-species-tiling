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
    #elements_to_check = list(range(len(element_info)))
    
    #while len(elements_to_check)>0:
    #    loop = find_loop(element_info,elements_to_check[0])
    #    loops.append(np.sort(loop))
    #    
    #    #remove all elements that are in the loop
    #    elements_to_check = np.setdiff1d(loop, elements_to_check)
    
    for i in range(len(element_info)):
        loop = find_loop(element_info,i)
        loops.append(np.sort(loop))
        if len(loop) not in [1,2,3,6]: return False, False

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

    #print('Number of loops:', np.shape(loops)[0])
    #print('Number of unique loops:', np.shape(unique_loops)[0])
    #print(unique_loops)
    
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
        if loop_len not in [3,6]: return False
        loop_props.append([loop_len, colors])
    
    #print(loop_props)

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

    #print(connected_int_matrix)
    return connected_int_matrix

def remove_degenerate_matrices(interaction_matricies,k):

    non_degen_list = []
    degen_list = interaction_matricies.copy()
    
    degen_elements = []
    for degen in degen_list:
        degen_elements.append(get_element_info(degen))
        
    
    #print(check_same(degen_list[2], degen_list[4]))
    
    #for dege_m in degen_list:
    #    draw_int_matrix(dege_m, show=True)
    
    #for M in degen_list:
    while len(degen_list)>0:
        M = degen_list[0]
    
        degenerate_M = []
        
        def add_degen_M(M_add):
            if len(degenerate_M)!=0:
                for dM in degenerate_M:
                    if check_same(M_add,dM): 
                        #print('did not add')
                        return None
                #print('added matrix')
                degenerate_M.append(M_add)
            else: 
                #print('added matrix')
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
                if check_same(dM, degen_list[i]): idx_to_remove.append(i)      
        
        for i in sorted(idx_to_remove, key=int, reverse=True):
            del(degen_copy[i])
        
        degen_list = degen_copy
        
        #print(idx_to_remove)
    
        print(len(degen_list))
        
    #print('degen_length', len(degen_list))
    print('non-degen length', len(non_degen_list))
    
    return non_degen_list
        
def generate_all_matrices(k):
    #start with a list of all the matrices we will eventually make
    interaction_matricies = []    

    M = np.zeros((3*k, 3*k))  #initialize our first matrix
    interaction_matricies.append(M)
    
    ### end initilization of variables ###
        
    #need a counter for the column set that we are looking at. A column set are the group of three columns that correspond to
    #the interactions of a single particle species
    col_set_indx = 0    
        
    while col_set_indx<k:
        
        print('current col', col_set_indx, ' of ', k, ', list length', len(interaction_matricies))
        #print(col_set_indx)
        
        #available columns in set
        col_indx = np.array([0,1,2])+3*col_set_indx
        #print(col_indx)
        
        new_int_matricies = [] 
        
        #look at each existing matrix to see what new matricies we could make from these
        for M in interaction_matricies:
            
            #we want to see which columns in the column set of M do not have elements already.
            col_to_fill = []
            for col in col_indx: #looping through available columns
                if any(M[col])>0: continue
                else: col_to_fill.append(col)
                   
            #if all columns of set are full then readd the matrix and go to next one, no additional elements to add
            if len(col_to_fill)==0:
                new_int_matricies.append(np.copy(M))
                continue
                
            #if there are elements to add, do that next
            #print("columns to fill", col_to_fill)
            
            #for each column we have k choices, and in the set these choices do not interfer with each other
            combs = list(itertools.product(range(k), repeat=len(col_to_fill)))  #get all combinations of [0,...,k-1]
            #print(combs)
            
            #for each combination we need to create a new matrix
            for i in combs:
                M_p = np.copy(M)
    
                #print(M_p)
                
                #fill as many rows as the length of col_to_fill
                row_to_fill = []
                
                for j in range(len(col_to_fill)):
                    #need to choose the correct row to fill. Recall that we will only be filling the diagonal elements of the row sets
                    row = col_to_fill[j]%3 + 3*i[j]
                    
                    #check if row is already filled
                    if any(M_p[row])>0: continue
                    else: row_to_fill.append(row)
                
                #if no rows to fill from this combination, skip it
                if len(row_to_fill)!=len(col_to_fill): continue
                
                #print('row to fill', row_to_fill)
                
                for col, row in zip(col_to_fill, row_to_fill):
                    M_p[col, row]=1
                    M_p[row, col]=1
                    
                if np.sum(M_p[col_indx, col_indx])<3:
                    if check_bad_lines(M_p,k):
                        new_int_matricies.append(M_p)
                        
                
    
                #for new_m in new_int_matrices:
                #    print(new_m)
    
        #after looping through all existing Ms we replace the list with the new list
        interaction_matricies = new_int_matricies.copy()
        
        process = psutil.Process(os.getpid())
        print('memory used (MB)', process.memory_info().rss/1000000.)
        
        
        #interate col_indx
        col_set_indx += 1
        
        #for m in interaction_matricies:
        #    print(interaction_matricies)
    return interaction_matricies

def remove_bad_line_interactions(M_list, k):
    #here we will remove elements that have lines forming along different lattice vectors.
    #can do this by checking each (col_set,row_set) if they have two elements
    
    valid_M = []
    
    #loop through all (row_set,col_set) and make list of which ones are two elements and which one is missing on the diagonal
    for M in M_list:

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

def check_bad_lines_elements(M, k):
    #here M is assumed to be the list of indices for non-zero elements
    
    set_props = []
    
    for i in range(k):
        for j in range(i,k):
            set_indx = [[0+3*i, 0+3*j],[1+3*i, 1+3*j],[2+3*i, 2+3*j]]
            #print(set_indx)
            diag_count = 0
            #pos = []
            
            for l in range(3):
                if set_indx[l] in M: diag_count +=1
                else: #pos.append(i) 
                    pos=l
                
            #print(set_diag)
            if diag_count==2:
                set_props.append(pos)

    if len(np.unique(set_props))<2:
        #print(set_props)
        return True
    
    else: 
        #print(set_props)
        #print(M)
        return False

def generate_m_from_elements(element_info, k):
    
    M = np.zeros((3*k, 3*k))
    
    if len(element_info)==0: return M
    
    for element in element_info:
        M[int(3*element[0]+element[2]), int(3*element[1]+element[2])]=1
        
    return M

def get_element_index(element_info,k):
    indices = []
    
    if len(element_info)==0: return [None, None]
    
    for element in element_info:
        indices.append([3*element[0]+element[2], 3*element[1] + element[2]])
    return indices

def generate_m_from_indices(indices,k):
    M = np.zeros((3*k, 3*k))
    for ind in indices: M[ind[0], ind[1]]=1
    return M

def generate_all_matrices_elements(k):
    #start with a list of all the matrices we will eventually make
    interaction_matricies = []    

    #M = np.zeros((3*k, 3*k))  #initialize our first matrix
    #interaction_matricies.append(get_element_info(M))
    
    interaction_matricies.append([])
    
    ### end initilization of variables ###
        
    #need a counter for the column set that we are looking at. A column set are the group of three columns that correspond to
    #the interactions of a single particle species
    col_set_indx = 0    
        
    while col_set_indx<k:
        
        print('current col', col_set_indx, ' of ', k, ', list length', len(interaction_matricies))
        #print(col_set_indx)
        
        #available columns in set
        col_indx = np.array([0,1,2])+3*col_set_indx
        #print(col_indx)
        
        new_int_matricies = [] 
        
        #look at each existing matrix to see what new matricies we could make from these
        for M in interaction_matricies:
            
            #generate a matrix to manipulate
            #M = generate_m_from_elements(e_M, k)  #want to rewrite this section to remove this step of generating a matrix from the elements
            
            indices = np.array(M)
            #print(indices.T)
            
            #we want to see which columns in the column set of M do not have elements already.
            col_to_fill = []

            if len(indices)==0: col_to_fill = list(col_indx)
            else:
                for col in col_indx:
                    if col in indices.T[0]: continue
                    else: col_to_fill.append(col)
            
            #col_to_fill = []
            #for col in col_indx: #looping through available columns
            #    if len(M[col])>0: continue
            #    else: col_to_fill.append(col)
                   
            #if all columns of set are full then readd the matrix and go to next one, no additional elements to add
            if len(col_to_fill)==0:
                new_int_matricies.append(M)
                continue
                
            #if there are elements to add, do that next
            #print("columns to fill", col_to_fill)
            
            #for each column we have k choices, and in the set these choices do not interfer with each other
            combs = list(itertools.product(range(k), repeat=len(col_to_fill)))  #get all combinations of [0,...,k-1]
            #print(combs)
            
            #for each combination we need to create a new matrix
            for i in combs:
                i_copy = M.copy()
    
                #print(type(i_copy))
                #print(M_p)
                
                #fill as many rows as the length of col_to_fill
                row_to_fill = []
                
                if len(indices)==0:
                    for j in range(len(col_to_fill)):
                        row_to_fill.append(col_to_fill[j]%3 + 3*i[j])
                
                else:
                    for j in range(len(col_to_fill)):
                        #need to choose the correct row to fill. Recall that we will only be filling the diagonal elements of the row sets
                        row = col_to_fill[j]%3 + 3*i[j]
                    
                        if row in indices.T[1]: continue
                        else: row_to_fill.append(row)
                        
                        #check if row is already filled
                        #if any(M_p[row])>0: continue
                        #else: row_to_fill.append(row)
                
                #if no rows to fill from this combination, skip it
                if len(row_to_fill)!=len(col_to_fill): continue
                
                #print('row to fill', row_to_fill)
                
                for col, row in zip(col_to_fill, row_to_fill):
                    i_copy.append([col,row])
                    if row!=col: i_copy.append([row,col])
                
                #for col, row in zip(col_to_fill, row_to_fill):
                #    M_p[col, row]=1
                #    M_p[row, col]=1
                
                diag_count=0
                for col in col_indx:
                    if [col,col] in i_copy: diag_count+=1
            
                if diag_count<3:
                    if check_bad_lines_elements(i_copy, k):
                        new_int_matricies.append(i_copy)
            
                        #test_M = generate_m_from_indices(i_copy, k)
                        #
                        #if np.sum(test_M[col_indx, col_indx])<3:
                        #    if check_bad_lines(test_M,k): print('correct bad line check')
                        #    else: 
                        #        print(test_M)
                        #        print(i_copy)
                                
                                
                #if np.sum(M_p[col_indx, col_indx])<3: #removes only self-interacting triangles
                #    if check_bad_lines(M_p,k): #removes matrices that have intersecting line patterns
                #        new_int_matricies.append(get_element_info(M_p))
    
                #for new_m in new_int_matrices:
                #    print(new_m)
    
        #after looping through all existing Ms we replace the list with the new list
        interaction_matricies = new_int_matricies.copy()
        
        process = psutil.Process(os.getpid())
        print('memory used (MB)', process.memory_info().rss/1000000.)
        
        #interate col_indx
        col_set_indx += 1
        
        #for m in interaction_matricies:
        #    print(interaction_matricies)
    return interaction_matricies

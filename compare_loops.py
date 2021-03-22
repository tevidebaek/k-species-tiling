import numpy as np
import int_matrix_r_func as imf
import matplotlib.pyplot as plt

k = 9

valid_matrices = np.load('valid_mat_list_k'+str(k)+'.npy')
save_name = 'final_mat_list_k'+str(k)+'.npy'

save_src = './k'+str(k)+'-valid-tiling/'
save_csv = 'k'+str(k)+'-intmatrices.csv'
save_to_csv = True
save_to_img = True

print(len(valid_matrices))

#loop through the list and generate the loop properties for each matrix. The things we care about are:
#       number of loops
#       length of each loop
#       how many loops each color is in
'''
loop_props = []

for mat in valid_matrices:
    unique_loops, element_info = imf.find_all_loops(mat)

    loop_colors = np.zeros(k)
    loop_lengths = []

    for loop in unique_loops:
        loop_len, colors = imf.get_loop_props(loop, element_info)

        for c in colors:
            loop_colors[int(c)]+=1

        loop_lengths.append(loop_len)

    #properties to add to list: [number of loops, [loop lengths], [number of loops for each color]]

    loop_colors = list(loop_colors)

    loop_props.append([len(unique_loops), loop_lengths, loop_colors])

def check_same_props(p1, p2):
    #each is structured: [number of loops, [loop lengths], [number of loops for each color]]
    if p1[0]==p2[0]:
        p1[1].sort()
        p2[1].sort()
        if p1[1]==p2[1]:
            p1[2].sort()
            p2[2].sort()
            if p1[2]==p2[2]:
                return True
    return False

unique_mat_props = []
cnt = 0
for prop in loop_props:
    
    if len(unique_mat_props)==0: unique_mat_props.append(prop)

    new_prop = True
    cnt_check = 0
    for un_prop in unique_mat_props:
        
        if check_same_props(prop, un_prop):
            new_prop=False
            print(cnt_check)
            continue
        cnt_check+=1

    if new_prop: unique_mat_props.append(prop)
    print(cnt, new_prop)
    cnt+=1

print(unique_mat_props)
print(len(unique_mat_props))
'''

#another option is to look at the graph (these are technically called multigraphs) structure of the
#what colors can be next to each other to do this is pretty simple. Within each column_set just see 
# which row-sets are occupied. This gives the so if col_set=1 had elements in row_set=[0,1,2] then 
# there would be a node for color 1 with edges going to color [0,2] and a loop on the node for color 1. 
#We will be using `igraph` to check for isomorphisms of these graphs. `igraph` can only do simple 
# graphs so to account for loops we will color those nodes differently:
#       0 loops = black
#       1 loop = blue
#       2 loops = red
#       3 loops = purple #this case should not happen but included for robustness
#If there are multiple edges between two nodes (as in the case for multi colored lines) then
#the edge will be weighted by the sum of the number of edges.

import igraph
import time

start = time.time()

def construct_graph(element_info, k):
    #the element_info for each element contains the information of (c_set, r_set, side)
    #we can just use the c_set, r_set info

    loop_count = np.zeros(k).astype(int)

    G = igraph.Graph()

    G.add_vertices(k)

    for e in element_info:
        e = e.astype(int)
        if e[0]==e[1]: loop_count[int(e[0])]+=1
        else: G.add_edges([(e[0], e[1])])
    
    #add a weight attribute to each edge
    for i in range(len(G.es)):
        G.es[i]['weight']=1

    #print(loop_count)

    for i in range(k):
        G.vs[i]['color'] = int(loop_count[i])

    G.simplify(loops=False, combine_edges='sum')

    return G

matrix_graphs = []

for mat in valid_matrices:
    unique_loops, element_info = imf.find_all_loops(mat)

    matrix_graphs.append(construct_graph(element_info,k))


unique_graphs = []
unique_mats = []

for G, M in zip(matrix_graphs, valid_matrices):

    if len(unique_graphs)==0: 
        unique_graphs.append(G)
        unique_mats.append(M)

    else:
        new_mat = True
        for H in unique_graphs:
            if G.isomorphic_vf2(H, color1=G.vs["color"], color2=H.vs["color"], edge_color1=G.es["weight"], edge_color2=H.es["weight"]):
                new_mat = False
                continue
                
        if new_mat:
            unique_graphs.append(G)
            unique_mats.append(M)

print(len(unique_graphs))

def plot_graph(g):
    #get the total number of nodes
    colors = ['black', 'blue', 'red', 'purple']

    node_num = len(g.vs)

    #generate x,y points along a circle
    theta = np.linspace(0,2*np.pi, node_num+1)
    xs, ys = np.cos(theta), np.sin(theta)

    #plot all of the nodes:
    for i in range(len(g.vs)):
        color = colors[g.vs[i]['color']]
        plt.scatter(xs[i], ys[i], c=color)

    edges = g.get_edgelist()
    for i in range(len(edges)):
        weight = g.es[i]['weight']
        init, final = edges[i][0], edges[i][1]
        plt.plot([xs[init],xs[final]], [ys[init],ys[final]], lw=weight, c='k')

'''
#try to draw the graphs that are generated by this process
for g in unique_graphs:
    plt.figure(figsize=(3,3))
    plot_graph(g)
    plt.xlim(-1.1,1.1)
    plt.ylim(-1.1,1.1)

plt.show()
'''

final_list = unique_mats

for i in range(len(final_list)):
    imf.draw_int_matrix(final_list[i], save_src, i,save=save_to_img, show=True)

    fname = save_src+'k'+str(k)+'-intmatrix-'+str(i)+'.csv'
    #print(fname)
    if save_to_csv:
        np.savetxt(fname, final_list[i].astype(int),delimiter=',')

np.save(save_name, unique_mats)

print(time.time() - start)
plt.show()
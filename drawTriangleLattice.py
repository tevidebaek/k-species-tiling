import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import int_matrix_r_func as imf
import os

#first we need to impor the interaction matrix file

# I will try to do this using python classes to learn about them!

class Triangle:
    #a triangle has several properties
    #    color : this determines the species of triangle
    #    orientation : when drawing it can be up or down
    #           -this is either +1 or -1
    #    sides : each side is unique and changes directions depending on the orientation
    #           -this is [0,1,2], used to call the interaction matrix
    #    location : where the center of the triangle is located
    #    neighbors : a list of which sides have adjacent triangles
    #    neighbor_colors : given a side this tells what color can go there
    #
    #triangles have a few functions:
    #    draw : given a canvas, they can be drawn onto them
    #    add_neighbor : if they have an available side, they can create a new triangle on that side

    global intMatrix_elements #this interaction Matrix will be generated in an other section of code but is the same for all instances of the Triangle class

    def __init__(self, orientation, color, location, added_from=-1):
        self.orientation = orientation
        self.color = color
        self.location = location
        self.sides = np.array([0,1,2])
        if added_from==-1:
            self.neighbors =[]
        else:
            self.neighbors = [added_from]
        self.find_allowed_neighbors()

    def add_neighbors(self):
        #this looks at which sides don't have neighbors a side and create a new triangle there
        #if len(self.neighbors)==0: sides_to_add = np.array([0,1,2])
        #else:
        #    sides_to_add = np.setdiff1d(self.neighbors,np.array([0,1,2]))

        sides_to_add = np.array([0,1,2])

        new_triangles = []

        if self.orientation==1:
            side_vec = np.array([(0,-np.sqrt(3)/2.), (3/4., np.sqrt(3)/4.), (-3/4., np.sqrt(3)/4.)])
        else:
            side_vec = np.array([(0,np.sqrt(3)/2.), (-3/4., -np.sqrt(3)/4.), (3/4., -np.sqrt(3)/4.)])

        #print(sides_to_add)

        for sides in sides_to_add:
            new_orientation = int(-1*self.orientation)
            new_color = int(self.neighbor_colors[sides])
            new_location = (self.location[0] + side_vec[sides][0], self.location[1] + side_vec[sides][1])
            new_added_from = sides

            new_triangles.append(Triangle(new_orientation, new_color, new_location, new_added_from))

        return new_triangles

    def find_allowed_neighbors(self):
        #print('looking for neightbors')
        #this looks at the global interaction matrix and finds which colors are allowed on which sides
        loc_self = np.where(intMatrix_elements.T[0]==self.color)
        sides = intMatrix_elements.T[2][loc_self]
        neighbor_colors = intMatrix_elements.T[1][loc_self]

        allowed_neighbors = np.array([sides,neighbor_colors]).T

        n_colors = [-1,-1,-1]
        for a_n in allowed_neighbors:
            #print(a_n)
            n_colors[int(a_n[0])] = int(a_n[1])

        self.neighbor_colors = n_colors

    def print_data(self):
        print('orientation:', self.orientation)
        print('color:', self.color)
        print('location:', self.location)
        print('allowed neighbors:', self.neighbor_colors)
        #print('new triangles to add:', self.new_triangles)

    def draw_triangle(self, axs, color):
        #get the cooreds of the vertices
        if self.orientation==1:
            vecs = np.array([(0,np.sqrt(3)/2.), (-3/4., -np.sqrt(3)/4.), (3/4., -np.sqrt(3)/4.)])
        else:
            vecs = np.array([(0,-np.sqrt(3)/2.), (3/4., np.sqrt(3)/4.), (-3/4., np.sqrt(3)/4.)])
        
        x, y = self.location[0], self.location[1]
        v1 = [x + vecs[0][0], y + vecs[0][1]]
        v2 = [x + vecs[1][0], y + vecs[1][1]]
        v3 = [x + vecs[2][0], y + vecs[2][1]]

        tri_patch = patches.Polygon([v1, v2, v3], True, facecolor=color, edgecolor='k', lw=1, alpha=.8)
        axs.add_patch(tri_patch)

        #l1_x, l1_y = [v1[0], v2[0]], [v1[1], v2[1]]
        #l2_x, l2_y = [v2[0], v3[0]], [v2[1], v3[1]]
        #l3_x, l3_y = [v3[0], v1[0]], [v3[1], v1[1]]

        #axs.plot(l1_x, l1_y, 'k', )
        #axs.plot(l2_x, l2_y, 'k')
        #axs.plot(l3_x, l3_y, 'k')


    #need a way to check if an added triangle has any other neighbors

#first thing to do is load in the interaction matrix and get a list of elements from it
src = './k9-valid-tiling/'

intMat_files = []

for root, dir, files in os.walk(src):
    for filename in files:
        if '.csv' in filename:
            intMat_files.append(os.path.join(src,filename))

for filename in intMat_files:

    intMatrix = np.loadtxt(filename, delimiter = ',')

    intMatrix_elements = imf.get_element_info(intMatrix)
    #print(intMatrix_elements)

    #print(np.where(intMatrix_elements.T[0]==0))

    #make an initial Triangle object
    first_triangle = Triangle(1,0,(0,0))
    #first_triangle.print_data()

    triangle_list = []
    triangle_list.append(first_triangle)

    triangle_locations = []
    triangle_locations.append(first_triangle.location)

    new_triangles = first_triangle.add_neighbors()

    #for nt in new_triangles:
    #    nt.print_data()

    addition_count = 0
    addition_iterations=12

    while addition_count<addition_iterations:

        triangles_to_add = []  

        #for each new triangle check to see if any existing triangles are already in that location:
        for nt in new_triangles:
            new_location = nt.location

            add_it = True
            for t_l in triangle_locations:
                if (abs(new_location[0]-t_l[0])<0.05)&(abs(new_location[1]-t_l[1])<0.05): add_it=False
            for t_add in triangles_to_add:
                t_add_loc = t_add.location
                if (abs(new_location[0]-t_add_loc[0])<0.05)&(abs(new_location[1]-t_add_loc[1])<0.05): add_it=False
            if add_it: triangles_to_add.append(nt)
            #if new_location in set(triangle_locations): continue
            #else: triangles_to_add.append(nt)

        #now that we have a list of things to add, get a list of their neighbors and then add those to the list
        new_triangles = []

        for t_to_add in triangles_to_add:
            new_triangles.extend(t_to_add.add_neighbors())

            triangle_list.append(t_to_add)
            triangle_locations.append(t_to_add.location)

        addition_count+=1

    #should have the final list of triangles now
    #print(len(triangle_list))
    #print([c.color for c in triangle_list])

    #colors = ['r', 'b', 'g', 'orange']
    colors = ['#ff6666', '#4d79ff', '#00cc99', '#d279d2', '#ffad33', '#70dbdb', '#ff99bb', '#009999', 'k']

    fig, axs = plt.subplots(figsize=(4,4))
    for t in triangle_list:
        x, y = t.location[0], t.location[1]
        color = colors[int(t.color)]
        
        #plt.scatter(x,y,c=color)

        t.draw_triangle(axs, color)

    axs.set_xlim(-10,10)
    axs.set_ylim(-10,10)
    axs.axis('off')
    fig.savefig(filename[:-4]+'.png', dpi=400)
    print(filename)
    plt.close()


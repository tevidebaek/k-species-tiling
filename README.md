# k-species-tiling
This finds allowed colorings of triangular tiling patterns for k different colors. It does this by looking for valid interaction matrices for the different triangles. An interaction matrix is a (3k, 3k) matrix whose non-zero values indicate which sides of triangles are interacting, i.e. are allowed to touch each other in the tiling. There are additional restrictions that the triangles have distinguishable sides and these sides can only lie on a specific lattice vector of the triangular lattice.

## How to run the code

You will need to specify the following things at the beginning of the `int_matrix.py` file:
-`k` : how many species to calculate for
-`save_src` : the directory where you want to save any outputs
-`save_csv` : how you want to format the filename of output CSV file that has matrix elements
-`save_to_csv` : boolean you set if you want to save CSV output
-`save_img` : boolean you set if you want to save a PNG image of matrix output

## Restrictions on tilings

We want to find all allowed interaction matrices that follow restricitions that come about from requiring that the tiling can be used to generate self-assembling tubes. The physical restrictions are that:

1. All interactions are deterministic (one side has a unique interaction with another side)
2. Triangles have a distinguishable orientation (i.e. their sides are different) and the sides of trianlges must lie along a specific lattice vector (only rotations of 180 degrees allowed)
3. When walking around a vertex, you must come back to the original triangle.

### Impact of restrictions

In this code we generate interaction matrices which are (3k, 3k) matrices that have 3k non-zero elements. The restrictions above impose rules about which elements may or may not be non-zero. 

1. This means that each row or column must have one and only one element that is non-zero. It also imposes that this matrix is symmetric.
2. This forces that only like-numbered sides can interact. For the interaction matrix this means that in any sub-(3,3) block that corresponds to interactions between species i and j, only the diagonal elements can be non-zero.
3. We can translate walking around a vertex into some form of loops between matrix elements. Given an elements position we can find which other elements it is connected to. This restriction enforces that these loops must be length 3 or 6.

One last point that arises that is not directly from the restrictions but that we disallow are interaction matrices that contain sub-matrices that are valid tilings for m-species where `m<k`. We only want cases where all colors participate in a unique tiling.

## Code Structure

The code goes through the following steps to construct all valid, non-degenerate interaction matrices for a given number of colors.

1. Generate all possible matrices that obey restrictions 1 and 2. At the moment we do not check restriction 3 because the matrices are bulit up in a tree structure and we can only check 3 once all elements are filled. This step also has some checks that remove obviously incorrect tilings when they appear. This is done with the `generate_all_matrices` function.
2. Check that all loops in the remaining matrices are only lengths 3 or 6. This step also checks that all colors participate in the tiling. This is done with the `valid_coloring` function.
3. At this point only valid tilings remain, but many are degenerate with each other. This means they are the same up to a permutation of the colors or a shift of the indicies (corresponding to a +60 or -60 degree rotation of the pattern). Here we find all degenerate matrices and keep only one of them in our final list. This is done with the `remove_degenerate_matrices` function
4. Lastly, we go through the final matrices and save them to a CSV file and their images if applicable.

## Issues that need resolving

The biggest concern right now is that when k gets large the list of matrices generated in the first step becomes excessively large and the code becomes RAM limited. Possible ways to fix this:

1. Find a more clever algorithim that does not need to generate all possible matricies.
2. Reduce the memory needed to store all these matrices
3. Implement a system to store matrices on disk (this would slow down the code but actually allow it to finish)

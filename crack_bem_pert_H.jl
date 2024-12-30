# The function makes the simulation of Boundary Element Method (BEM)
# for a fault in a material
# It will calculates the stress and vector properties of the fault
# which will be used to analyse the mechanical behaviour of the faults

# Flow of the function:
# Calculate the midpoints and segment lengths
# Define patch properties such as the dimensions, strike-dip orientation
# normal, strike-slip, and dip-slip vectors
# Group all of the properties to a matrix
# Use the matrix for Green's functions and stress calculation
using LinearAlgebra

function crack_bem_pert_H(L, mu, nu, x, y)

    # inputs:
    # L = length of fault
    # mu = shear modulus
    # nu = Poisson's ration
    # x, y = coordinates of the fault

    # outputs:
    # G = matrix representing Green's functions for displacements caused by the fault movement
    # GN = Green's functions related to normal stress 
    # PU = means stress values (normal+shear) over each patch
    # m = matrix containing geometric properties, such as length, width
    # depth, strike, dip, midpoint coordinates
    # vn = normal vectors to the fault at each patch
    # vs = strike-slip direction vectors 
    # vd = dip-slip direction vectorss
    # xloc = 3D coordinates of patch centers

    # mindpoints of fault segments
    xc = 0.5 .* (x[2:end] .+ x[1:end-1]) # x-coordinates of midpoints
    yc = 0.5 .* (y[2:end] .+ y[1:end-1]) # y-coordinates of midpoints

    # segment lengths
    ls = sqrt.((x[2:end] .- x[1:end-1]).^2 .+ (y[2:end] .- y[1:end-1]).^2)

    # width and depth
    W = ls .* 0 .+ L * 100 # set width for each segment
    D = ls .* 0 .+ L * 100 # set depth for each segment

    # fault orientation (strike and dip)
    patch_strike = 90 .- atan.(diff(y) ./ diff(x)) .* (180 / pi) # strike angle in degrees
    patch_dip = 90 .+ ls .* 0 # dip angle is constant (90 degrees)

    # Lame's first parameter
    λ = 2 * mu * nu / (1 - 2 * nu)

    # number of fault elements
    N = length(xc)

    # initialise normal, strike-slip, and dip-slip vectors
    #vn = [cosd.(patch_strike) .* sind.(patch_dip); 
   # -sind.(patch_strike) .* sind.(patch_dip); 
   # cosd.(patch_dip)]
    vn = [
        cosd(patch_strike[1]) * sind(patch_dip[1]) cosd(patch_strike[2]) * sind(patch_dip[2]) cosd(patch_strike[3]) * sind(patch_dip[3])
        -sind(patch_strike[1]) * sind(patch_dip[1]) -sind(patch_strike[2]) * sind(patch_dip[2]) -sind(patch_strike[3]) * sind(patch_dip[3]) 
        cosd(patch_dip[1]) cosd(patch_dip[2]) cosd(patch_dip[3])
    ]
    vs = zeros(3,N) # strike-slip vectors
    vd = zeros(3,N) # dip-slip vectors

    # calculate normal, strike-slip, and dip-slip vectors for each segment 
    # loop to go over each fault segment N 
    # the product from the loop will be used to define the coordinate system of the fault 
    for iN in 1:N 
        # put the normal vector values to the fault plane
        # strike angle determines the orientation of the fault in the horizontal plane
        # dip angle determines the inclination of the fault plane
        # use trigonometric functions to change the value to vector unit in length 
        # x-component = cos(strike) sin(dip)
        # y-component = -sin(strike) sin(dip)
        # z-component = cos(dip)
        vn[:, iN] .= vn[:, iN] / norm(vn[:, iN])
        # variable vn will points outward, perpendicular to the fault plane which represents the fault plane's 
        # orientation

        # put the strike-slip vector to the fault plane (perpendicular to the normal vector)
        # take x and y components of the normal vector 
        # then rotate the 2D vector by 90 degrees counterclockwise
        # rotate 90 degrees counterclockwise
        vs[1:2, iN] .= [0 1; -1 0] * vn[1:2, iN] 
        # normalise and invert direction
        vs[:, iN] .= -vs[:, iN] / norm(vs[:, iN]) 
        # variable vs lies along the fault strike direction which defines the horizontal 
        # slip direction for strike-slip motion
        vn_vec = vn[:, iN]
        vs_vec = vs[:, iN]

        # put the dip-slip vector to the fault plane (orthogonal to normal vector and strike-slip vector)
        vd[:, iN] = cross(vn_vec, vs_vec) # dip-slip vector
        # variable vd lies along the fault dip direction which defines the vertical slip 
        # direction for the dip-slip motion 

        vd[:, iN] .= vd[:, iN] / norm(vd[:, iN])
        # normalize the vectors 

    end

# initilise a matrix m of size of 10 x N with zeros
m = zeros(10, N)

# Assign values to the rows of m
m[1, :] .= ls
m[2, :] .= W
m[3, :] .= D
m[4, :] .= patch_dip
m[5, :] .= patch_strike
m[6, :] .= xc
m[7, :] .= yc

# define xloc as a 3 x N matrix
xloc = [m[6, :]; m[7, :]; -m[3,:]]

# Initialise matrix G, GN and PU 
G = zeros(N, N)
GN = zeros(N, N)
PU = zeros(N, N)

# Create a vertical unit vector vss 
# vss represents the vertical direction of the components
vss = [0.0, 0.0, 1.0]

    # Loop over each segment of the fault with N as the total number of segments
    for i in 1:N 

        # Set single iteration loop where;
        for i2 in 0:0
            # m1 will extract the i-th column of m 
            # and the function copy will create new array to avoid modifying
            # the original matrix m
            m1 = copy(m[:,i]) 

            # Set the 8th element of m1 as 1 
            m1[8 + i2] = 1 

            # List stress components
            # The inputs to RDstressFS are
            # xloc[1:3, :] = 3D coordinates of the fault
            # m1[6:1] = fault parameters
            # mu, lambda = material properties 
            # Pc = type of stress calculation
            # The output will be matrix Stress that contain stress tensor components
            # at specific coordinates
            Stress, _ = RDstressFS(
                xloc[1,:], xloc[2,:], xloc[3,:], 
                m1[6], m1[7], m1[3], m1[1], m1[2],
                0, m1[4], m1[5], 90, 1, 0, mu, lambda, "Pc"
            )

            # Rearrange stress matrix
            # Function permutedims to transpose the rearranged stress tensor
            S = permutedims(Stress[:, [1,4,5,2,6,3]])

            # Compute traction vectors
            # First project the stress components to vn then combines all the projections from all the coordinates
            Tn = [
                sum(S[1:3, :] .* vn)
                sum(S[[2, 4, 5], :] .* vn)
                sum(S[[3, 5, 6], :] .* vn)
            ]

            # Update the matrix 
            # G contains the traction projected onto vertical direction
            G[:, 1] .= sum(Tn .* vss)
            # GN contains the traction projected onto fault-normal direction
            GN[:, i] .= sum(Tn .* vn)
            # PU stores the mean normal stress
            PU[:, i] .= -(S[1,:] + S[4,:] + S[6,:]) / 3


        end
    end 

end


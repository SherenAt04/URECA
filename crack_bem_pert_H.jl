function crack_bem_pert_H(L, mu, nu, x, y)
    # inputs:
    # L = length of fault
    # mu = shear modulus
    # nu = Poisson's ration
    # x, y = coordinates of the fault

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
    vn = zeros(3,N) # normal vectors
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
        vn[:, iN] = normalize([cosd(patch_strike[iN]) * sind(patch_dip[iN])
        - sind(patch_strike[iN]) * sind(patch_dip[iN]) cosd(patch_dip[iN])])
        # variable vn will points outward, perpendicular to the fault plane which represents the fault plane's 
        # orientation

        # put the strike-slip vector to the fault plane (perpendicular to the normal vector)
        # take x and y components of the normal vector 
        # then rotate the 2D vector by 90 degrees counterclockwise
        vs[1:2, iN] = [0 1; -1 0] * vn[1:2, iN] # rotate 90 degrees counterclockwise
        vs[:, iN] = -normalize(vs[:, iN]) # normalise and invert direction
        # variable vs lies along the fault strike direction which defines the horizontal 
        # slip direction for strike-slip motion

        # put the dip-slip vector to the fault plane (orthogonal to normal vector and strike-slip vector)
        vd[:, iN] = normalize(cross(vn[:, iN], vs[:, iN])) # dip-slip vector
        # variable vd lies along the fault dip direction which defines the vertical slip 
        # direction for the dip-slip motion 

    end

end
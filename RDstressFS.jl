function RDstressFS(X, Y, Z, XO, YO, depth, L, W, pluge, dip, strike, rake, slip, opening, mu, lambda, RefPoint)
    # Poisson's ratio
    nu = lambda / (lambda + mu) / 2 

    # Tensile-slip
    bx = opening 
    # Strike-slip
    by = slip * cosd(rake)
    # Dip-slip
    bz = slip * sind(rake)

    # Ensure x is a column vector
    X = X[:]
    # Ensure Y is a column vecctor
    Y = Y[:]
    # Ensure Z is a column vector
    Z = Z[:]

    # Rotation matrix
    Rz1 = [
        cosd(plunge) sind(plunge) 0
        -sind(plunge) cosd(plunge) 0
        0 0 1
    ]

    Ry = [
        cosd(dip) 0 sind(dip)
        0 1 0
        -sind(dip) 0 cosd(dip)
    ]

    Rz2 = [
        cosd(strike) sind(strike) 0 
        -sind(strike) cosd(strike) 0
        0 0 1
    ]

    Rt = Rz2 * Ry * Rz1

    # Coordinates for the points defining the fault
    Pt1 = [-W/2, L/2, 0]'
    Pt2 = [-W/2. - L/2, 0]'
    Pt3 = [W/2, - L/2, 0]'
    Pt4 = [W/2, L/2, 0]'

    

    


end 
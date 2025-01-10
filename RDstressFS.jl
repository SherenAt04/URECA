function RDstressFS(X, Y, Z, X0, Y0, depth, L, W, plunge, dip, strike, rake, slip, opening, mu, lambda, RefPoint)
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
    @show Rt
    
    # Coordinates for the points defining the fault
    Pt1 = [-W/2, L/2, 0]'
    Pt2 = [-W/2. - L/2, 0]'
    Pt3 = [W/2, - L/2, 0]'
    Pt4 = [W/2, L/2, 0]'

    # Set reference point based on the input argument 
    if RefPoint == "Pc"
        Ptr = [0, 0, 0]'
    elseif RefPoint == "P1"
        Ptr = [-W/2, L/2, 0]'
    elseif RefPoint == "P2"
        Ptr = [-W/2, - L/2, 0]'
    elseif RefPoint == "P3"
        Ptr = [W/2, - L/2, 0]'
    elseif RefPoint == "P4"
        Ptr = [W/2, L/2, 0]'
    elseif RefPoint in ["mP12", "mP21"]
        Ptr = [-W/2, 0, 0]'
    elseif RefPoint in ["mP23", "mP32"]
        Ptr = [0, - L/2, 0]'
    elseif RefPoint in ["mP34", "mP43"]
        Ptr = [W/2, 0, 0]'
    elseif RefPoint in ["mP14", "mP41"]
        Ptr = [0, L/2, 0]'
    else
        throw(ArgumentError("Undefined RD reference point!"))
    end 

    @show size(X0)
    @show size(Y0)
    @show size(depth)
    @show size(Ptr)
    @show size(Rt)

    # Calculate the displacement 
    Pr = [X0, Y0, -depth]' - Rt * Ptr
    P1 = Rt * Pt1 + Pr 
    P2 = Rt * Pt2 + Pr 
    P3 = Rt * Pt3 + Pr 
    P4 = Rt * Pt4 + Pr 

    return P1, P2, P3, P4

    # calculate unit strike, dip, and normal 
    eZ = [0.0, 0.0, 1.0]
    Vnorm = Rt * eZ
    Vstrike = [sin(deg2rad(strike)), cos(deg2rad(strike)), 0.0]
    Vdip = cross(Vnorm, Vstrike)

    Pm = (P1 + P2 + P3 + P4) / 4 

    # Transform coordinates and slip vector components from EFCS to RDCS
    p1 = zeros(3)
    p2 = zeros(3)
    p3 = zeros(3)
    p4 = zeros(3)

    At = hcat(Vnorm, Vstrike, Vdip)'
    x, y, z = CoordTrans(X .- Pm[1], Y .- Pm[2], Z .- Pm[3], At)
    p1 .= CoordTrans(P1[1] - Pm[1], P1[2] - Pm[2], P1[2] - Pm[3], At)
    p2 .= CoordTrans(P2[1] - Pm[1], P2[2] - Pm[2], P2[2] - Pm[3], At)
    p3 .= CoordTrans(P3[1] - Pm[1], P3[2] - Pm[2], P3[2] - Pm[3], At)
    p4 .= CoordTrans(P4[1] - Pm[1], P4[2] - Pm[2], P4[2] - Pm[3], At)

    # Calculate the unit vectors along RD sides in RDCS
    e12 = (p2 - p1) / norm(p2 - p1)
    e23 = (p3 - p2) / norm(p3 - p2)
    e34 = (p4 - p3) / norm(p4 - p3)
    e14 = (p4 - p1) / norm(p4 - p1)

    # Determine the best artifact-free configuration for each calculation point
    Rectmode = rectmodefinder(y, z, x, p1[2:3], p3[2:3], p4[2:3])
    casepLog = Rectmode .== 1
    casenLog = Rectmode .== -1
    casezLog = Rectmode .== 0

    xp = x[casepLog]
    yp = y[casepLog]
    zp = z[casepLog]

    xn = x[casenLog]
    yn = y[casenLog]
    zn = z[casenLog]

    # Configuration I 
    


end 
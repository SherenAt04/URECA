function RDstressFS(X, Y, Z, X0, Y0, depth, L, W, plunge, dip, strike, rake, slip, opening, mu, lambda, RefPoint)
   # This function transforms the strain and stress tensors for a fault system (FS) in a 3D system based on the rectangular dislocation (RD) model
   # The inputs are:
   # - fault geomtery parameters: coordinates of observation points (X, Y, Z), coordinates of the fault center (X0, Y0, depth), length and width of the fault plane (L, W), fault orientation angles in degrees (plunge, dip, strike), fault slip characteristics (rake, slip, opening)
   # - elastic properties: shear modulus (mu), first Lame parameter (lambda)
   # - Reference point to determine which points of the fault that can be the reference for calculations
   #   This includes the fault corners (P1, P2, P3, P4), center (PC) and midpoints of the fault edges

    # Poisson's ratio
    nu = lambda / (lambda + mu) / 2 

    # slip components 
    # Tensile-slip = normal opening of the fault
    bx = opening 
    # Strike-slip = horizontal movement along the fault strike
    by = slip * cosd(rake)
    # Dip-slip = vertical movement along the fault dip
    bz = slip * sind(rake)

    # Ensure x is a column vector
    X = X[:]
    # Ensure Y is a column vecctor
    Y = Y[:]
    # Ensure Z is a column vector
    Z = Z[:]

    # Rotation matrix
    # constructs the 3 rotation matrix to transform between coordinate systems
    # Rz1 = rotates about the z axis by the plunge angle
    # Ry = rotates about the y axis by the dip angle
    # Rz2 = rotates about the z-axis by the strike angle 
    # Rt (combined rotation matrix) transforms the fault coordinate system (FCS) to Earth-fixed coordinate system (EFCS)
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
    
    # fault geometry 
    # Coordinates for the points defining the fault
    # the fault is defined as a rectangle with Pt as the corners of the fault plane 
    Pt1 = [-W/2, L/2, 0]'
    Pt2 = [-W/2. - L/2, 0]'
    Pt3 = [W/2, - L/2, 0]'
    Pt4 = [W/2, L/2, 0]'

    # Set reference point based on the input argument 
    if RefPoint == "Pc"
        Ptr = [0, 0, 0]'
        # if the Refpoint is the center of he fault plane, then the reference point is set to the origin of the local fault coordinate system 
    elseif RefPoint == "P1"
        Ptr = [-W/2, L/2, 0]' # -W/2 = half the width to the left, L/2 = half the length upward
        # if the Refpoint is the top-left corner of the fault plane in the fault coordinate system, then the reference point will be relative to the fault center 
    elseif RefPoint == "P2" 
        Ptr = [-W/2, - L/2, 0]' # -W/2 = half the width to the left, -L/2 = half the length downward
        # if the Refpoint is the bottom-left corner of the fault plane, then the reference point will be relative to the fault center
    elseif RefPoint == "P3"
        Ptr = [W/2, - L/2, 0]' # W/2 = half the width to the right, -L/2 = half the length downward 
        # if the Refpoint is the bottom-right corner, the reference point will be relative to the fault center 
    elseif RefPoint == "P4"
        Ptr = [W/2, L/2, 0]' #W/2 = half the width to the right, L/2 = half the length upward 
        # if the Refpoint is the top-right corner, the reference point will be relative to the fault center 
    elseif RefPoint in ["mP12", "mP21"]
        Ptr = [-W/2, 0, 0]' # -W/2 = half the width to the left, 0 = center of the edge in the length direction 
        # if the refpoint is the midpoint of the left edge of the fault plane (either P1 or P2), then the reference point will be relative to the left side of the fault 
    elseif RefPoint in ["mP23", "mP32"]
        Ptr = [0, - L/2, 0]' # 0 = center of the edge in the width direction, -L/2 = half the length downward 
        # if the refpoint is the midpoint of the bottom edge of the fault plane (either P2 and P3), then the reference point will be relative to the downward of the fault 
    elseif RefPoint in ["mP34", "mP43"]
        Ptr = [W/2, 0, 0]' # W/2 = half the width to the right, 0 = center of the edge in the length direction 
        # if the refpoint is the midpoint of the right vertical edge of the fault plane (either P3 and P4), then the reference point will be relative to the right side of the fault
    elseif RefPoint in ["mP14", "mP41"]
        Ptr = [0, L/2, 0]' # 0 = center of the edge in the width direction, L/2 = half the length upward 
        # if the refpoint is the midpoint of the top horizontal edge of the fault plane (either P1 and P4), then the reference point will be relative to the upward of the fault 
    else
        throw(ArgumentError("Undefined RD reference point!"))
        # if the refpoint is not one of the options, the output will be error and stop the loop
    end 

    @show size(X0)
    @show size(Y0)
    @show size(depth)
    @show size(Ptr)
    @show size(Rt)

    # Transform the fault corners coordinate to global EFCS 
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

    # Transform coordinates and slip vector components from EFCS to RDCS by using CoordTrans function 
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

    # strain calculation
    # for each configuration, function RDSetupS will calculate the strain tensor components at observation points such as the contributions from angular dislocations along fault edges and summing the contribution from all the edges
    # for undefined points (casezLog), the strain components will be set to NaN 

    # Configuration I 
    if count(casepLog) != 0
        # calculate first angular dislocation contribution
        Exx1Tp, Eyy1Tp, Ezz1Tp, Exy1Tp, Exz1Tp, Eyz1Tp = RDSetupS(xp, yp, zp, bx, by, bz, nu, p1, e14)
        # calculate second angular dislocation contribution
        Exx2Tp, Eyy2Tp, Ezz2Tp, Exy2Tp, Exz2Tp, Eyz2Tp = RDSetupS(xp, yp, zp, bx, by, bz, nu, p1, -e12)
        # calculate third angular dislocation contribution
        Exx3Tp, Eyy3Tp, Ezz3Tp, Exy3Tp, Exz3Tp, Eyz3Tp = RDSetupS(xp, yp, zp, bx, by, bz, nu, p1, -e23)
        # calculate fourth angular dislocation contribution
        Exx4Tp, Eyy4Tp, Ezz4Tp, Exy4Tp, Exz4Tp, Eyz4Tp = RDSetupS(xp, yp, zp, bx, by, bz, nu, p1, -e34)
    end

    # Configuration II 
    if count(casenLog) != 0
        # calculate first angular dislocation contribution
        Exx1Tn, Eyy1Tn, Ezz1Tn, Exy1Tn, Exz1Tn, Eyz1Tn = RDSetupS(xn, yn, zn, bx, by, bz, nu, p1, -e14)
        # calculate second angular dislocation contribution
        Exx2Tn, Eyy2Tn, Ezz2Tn, Exy2Tn, Exz2Tn, Eyz2Tn = RDSetupS(xn, yn, zn, bx, by, bz, nu, p1, e12)
        # calculate third angular dislocation contribution
        Exx3Tn, Eyy3Tn, Ezz3Tn, Exy3Tn, Exz3Tn, Eyz3Tn = RDSetupS(xn, yn, zn, bx, by, bz, nu, p1, e23)
        # calculate fourth angular dislocation contribution
        Exx4Tn, Eyy4Tn, Ezz4Tn, Exy4Tn, Exz4Tn, Eyz4Tn = RDSetupS(xn, yn, zn, bx, by, bz, nu, p1, e34)
    end

    # calculate the strain tensor components in RDCS
    if count(casepLog) != 0 
        exx[casepLog] .= Exx1Tp .+ Exx2Tp .+ Exx3Tp .+ Exx4Tp
        eyy[casepLog] .= Eyy1Tp .+ Eyy2Tp .+ Eyy3Tp .+ Eyy4Tp
        ezz[casepLog] .= Ezz1Tp .+ Ezz2Tp .+ Ezz3Tp .+ Ezz4Tp
        exy[casepLog] .= Exy1Tp .+ Exy2Tp .+ Exy3Tp .+ Exy4Tp
        exz[casepLog] .= Exz1Tp .+ Exz2Tp .+ Exz3Tp .+ Exz4Tp
        eyz[casepLog] .= Eyz1Tp .+ Eyz2Tp .+ Eyz3Tp .+ Eyz4Tp
    end

    if count(casenLog) != 0
        exx[casenLog] .= Exx1Tn .+ Exx2Tn .+ Exx3Tn .+ Exx4Tn
        eyy[casenLog] .= Eyy1Tn .+ Eyy2Tn .+ Eyy3Tn .+ Eyy4Tn
        ezz[casenLog] .= Ezz1Tn .+ Ezz2Tn .+ Ezz3Tn .+ Ezz4Tn
        exy[casenLog] .= Exy1Tn .+ Exy2Tn .+ Exy3Tn .+ Exy4Tn
        exz[casenLog] .= Exz1Tn .+ Exz2Tn .+ Exz3Tn .+ Exz4Tn
        eyz[casenLog] .= Eyz1Tn .+ Eyz2Tn .+ Eyz3Tn .+ Eyz4Tn
    end

    if count(casezLog) != 0
        exx[casezLog] .= NaN
        eyy[casezLog] .= NaN
        ezz[casezLog] .= NaN
        exy[casezLog] .= NaN
        exz[casezLog] .= NaN
        eyz[casezLog] .= NaN
    end

    # transform the strain tensor components from RDCS into EFCS
    Exx, Eyy, Ezz, Exy, Eyz = TensTrans(exx, eyy, exy, exz, eyz, hcat(Vnorm, Vstrike, Vdip))

    # calculate the stress tensor components in EFCS 
    Sxx = 2 * mu * Exx .+ lambda * (Exx .+ Eyy .+ Ezz)
    Syy = 2 * mu * Eyy .+ lambda * (Exx .+ Eyy .+ Ezz)
    Szz = 2 * mu * Ezz .+ lambda * (Exx .+ Eyy .+ Ezz)
    Sxy = 2 * mu * Exy 
    Sxz = 2 * mu * Exz 
    Syz = 2 * mu * Eyz 

    # output 
    Strain = [Exx, Eyy, Ezz, Exy, Exz, Eyz]
    Stress = [Sxx, Syy, Szz, Sxy, Sxz, Syz]
    


end 
# CoordTrans function
function CoordTrans(x1, x2, x3, A)
    # This function transforms the coordinates of the vectors from x1x2x3 coordinates to X1X2X3 coordinates
    # The input will include x1, x2, x3 and A
    # A is the transformation matrix 

    x1 = vec(x1)
    x2 = vec(x2)
    x3 = vec(x3)

    mat = hcat(vec(x1), vec(x2), vec(x3))

    r = mat * A' 
    X1 = r[1,:]'
    X2 = r[2, :]'
    X3 = r[3, :]'

    return X1, X2, X3


end
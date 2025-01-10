# CoordTrans function
function CoordTrans(x1, x2, x3, A)
    x1 = vec(x1)
    x2 = vec(x2)
    x3 = vec(x3)

    mat = [
        x1'
        x2'
        x3'
    ]

    r = A * mat 
    X1 = r[1,:]'
    X2 = r[2, :]'
    X3 = r[3, :]'

    return X1, X2, X3


end
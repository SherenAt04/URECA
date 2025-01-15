# rectmodefinder function 

function rectmodefinder(x,y,z,p1,p2,p3,p4)
# This function shows the value of angular dislocations for the calculation
# The inputs are x, y, z which corresponds to y, z, x coordinates in RDCS
# p1, p2, p3, p4 are matrixes that represent the y and z coordinates of the RD vertices in RDCS

    x = vec[x]
    y = vec[y]
    z = vec[z]

    p1 = vec[p1]
    p2 = vec[p2]
    p3 = vec[p3]
    p4 = vec[p4]

    pm = [p1 + p2 + p3 + p4] / 4
    r21 = p1 - p2
    e21 = r21 / sqrt[r21' * r21]
    r41 = p1 - p4
    e41 = r41 / sqrt[r41' * r41]

    x = x - pm[1]
    y = y - pm[2]

    A = [e21 e41]
    S = [
        x'
        y'        
        ]

    r = A' * S 
    x = r[1, :]'
    y = r[2, :]'

    r = A' * [ p1 - pm, p2 - pm, p3 - pm, p4 - pm ]
    P1 = r[:, 1]
    P2 = r[:, 2]
    P3 = r[:, 3]
    P4 = r[:, 4]





end
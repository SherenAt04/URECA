function RDSetupS(x,y,z,bx,by,bz,nu,RDVertex,SideVec)
    # RDSetupS transforms coordinates of the calculation points and slip vector
    # Components from ADCS into RDCS, then calculates the strains in ADCS
    # and transforms them into RDCS

    #transformation matrix
    A = [SideVec[3] -SideVec[2]; SideVec[2:3]']
    
    #transform coordinates of the calculation points from RDCS into ADCS
    r1 = A * [y' .- RDVertex[2]; z' .- RDVertex[3]]
    y1 = r1[1, :]'
    z1 = r1[2, :]'

    #transform the in-plane slip vector components from RDCS into ADCS 
    r2 = A * [by; bz]
    by1 = r2[1, :]'
    bz1 = r2[2, :]'

    #calculate strains associated with an angular dislocation in ADCS
    exx, eyy, ezz, exy, exz, eyz = AngDisStrain(x, y1, z1, -pi/2, bx, by1, bz1, nu)

    #transform strains from ADCS into RDCS
    # 3x3 transformation matrix
    B = [1 0 0; zeros(2,1) A'] 
    exx, eyy, ezz, exy, exz, eyz = TensTrans(exx, eyy, ezz, exy, exz, eyz, B)

    return exx, eyy, ezz, exy, exz, eyz 
end
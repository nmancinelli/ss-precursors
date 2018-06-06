def sph_loc(flat1,flon1,delta,azi):
    from numpy import sin, cos, arccos, arctan2
    if delta==0.:
        flat2=flat1
        flon2=flon1
        return
    pi=3.141592654
    raddeg=pi/180.
    delr=delta*raddeg
    azr=azi*raddeg
    theta1=(90.-flat1)*raddeg
    phi1=flon1*raddeg      
    ctheta2=sin(delr)*sin(theta1)*cos(azr)+cos(theta1)*cos(delr)
    theta2=arccos(ctheta2)
    if theta1 == 0.:
         phi2=azr
    elif theta2 == 0.:
         phi2=0.0
    else:
        sphi2=sin(delr)*sin(azr)/sin(theta2)
        cphi2=(cos(delr)-cos(theta1)*ctheta2)/(sin(theta1)*sin(theta2))
        phi2=phi1+arctan2(sphi2,cphi2)
    flat2=90.-theta2/raddeg
    flon2=phi2/raddeg
    if (flon2 > 360.):
        flon2=flon2-360.
    if (flon2 < 0.):
        flon2=flon2+360.
        
    return flat2, flon2
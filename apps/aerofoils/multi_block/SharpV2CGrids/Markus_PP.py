bc4=bound4()
[dxdxi,dydxi,dxdeta,dydeta,corrf]=getmet(nxt,nyt,nzt,dz,TE,bc4)

dzdzeta=dz
det=dz*(dxdxi*dydeta-dxdeta*dydxi)
aidet=1/det
dydx=dydxi/dxdxi

ss=np.zeros(nxt)
dss=np.zeros(nxt)
th_=np.zeros(nxt)
S_=np.zeros(nxt)
S2_=np.zeros(nxt)
for i in range(0,nxt):
    th_[i]=np.arctan(dydx[i])
    if(dxdxi[i]>0):
        S_[i]=1
    else:
        S_[i]=-1
    if(dydxi[i]>0):
        S2_[i]=1
    else:
        S2_[i]=-1
    if(i>0):
        ss[i]=ss[i-1]+((x[i,0]-x[i-1,0])**2+(y[i,0]-y[i-1,0])**2)**0.5
        dss[i]=m.fabs(ss[i]-ss[i-1])

#    #-----------------------------------------------------------------  calc Coeff

    qin=np.zeros((nxt,nyt))
    qin[:,0]=q[:,0,5]/corrf
    dudxi=d1xi_2(qin,nxt,nyt,bc4)

    dvdeta=q[:,:,2]/corrf
    dudeta=q[:,:,1]/corrf
    dudy=np.zeros((nxt,nzt))
    dvdx=np.zeros((nxt,nzt))
    for k in range(0,nzt):
        dudy[:,k]=aidet[:]*(dudeta[:,k]*dxdxi[:]*dzdzeta)
        dvdx[:,k]=aidet[:]*(-dvdeta[:,k]*dydxi[:]*dzdzeta)

    sspsq=gam*(gam-1)*(q[:,:,4]/q[:,:,0])-0.5*0 # The velocities are neglected as the
    T=sspsq*(mach**2)
    mu=(T**(3/2)*(1+110.4/268.67))/(T+110.4/268.67)/re
    P=1/(gam*mach**2)*q[:,:,0]*T
    pinf=1/(gam*mach**2)

    tw=np.zeros((nxt,nzt))
    i=0
    for k in range(0,nzt):
        tw[i,k]=S_[i]*mu[i,k]*(dudy[i,k]*np.abs(np.cos(th_[i]))-dvdx[i,k]*np.abs(np.sin(th_[i])))

    for i in range(1,nxt):
        for k in range(0,nzt):
            dlts=ss[i]-ss[i-1]
            fa = -S_[i]*(P[i,k]-pinf)*np.abs(np.cos(th_[i]))
            fb = -S_[i-1]*(P[i-1,k]-pinf)*np.abs(np.cos(th_[i-1]))
            Cl[ind,k]=Cl[ind,k]+0.5*(fa+fb)*dlts

            fa = S2_[i]*P[i,k]*abs(np.sin(th_[i]))
            fb = S2_[i]*P[i-1,k]*abs(np.sin(th_[i-1]))
            Cdp[ind,k] = Cdp[ind,k]+0.5*(fa+fb)*dlts

            fa = S_[i]*mu[i,k]*(dudy[i,k]*np.abs(np.cos(th_[i]))-dvdx[i,k]*np.abs(np.sin(th_[i])))
            fb = S_[i]*mu[i-1,k]*( dudy[i-1,k]*np.abs(np.cos(th_[i-1]))-dvdx[i-1,k]*np.abs(np.sin(th_[i-1])) )
            Cdf[ind,k] = Cdf[ind,k]+0.5*(fa+fb)*dlts
            tw[i,k]=S_[i]*mu[i,k]*(dudy[i,k]*np.abs(np.cos(th_[i]))-dvdx[i,k]*np.abs(np.sin(th_[i])))

    #cp_=(P-pinf)/(gam*0.5*pinf*mach**2)
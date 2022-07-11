def getmet(nx,ny,nz,dz,TE,bc4):
    import struct
    import numpy as np
    import matplotlib.pyplot as plt
    from deriv import d1xi_2,d1eta_2 

#-----------------------------------------------read full grid
    f = open('Bl1.dat', 'r')
    header=f.readline()
    g = open('Bl3.dat', 'r')
    headerg=g.readline()
    if(headerg!=header):
        print('ERROR!!!!!!!!!!!!!!!!!!!!!!!!')
    else:
        print('read full grid')
    header=header.strip()
    columnsh = header.split()
    nxt=int(columnsh[0])
    nyt=int(columnsh[1])
    x1=np.zeros((1,nyt))
    y1=np.zeros((1,nyt))
    x3=np.zeros((1,nyt))
    y3=np.zeros((1,nyt))
    for indx in range(0,1):
        for indy in range(0,nyt):
            line = f.readline()
            line = line.strip()
            columns = line.split()
            x1[indx,indy] = float(columns[0])
            y1[indx,indy] = float(columns[2])
            line = g.readline()
            line = line.strip()
            columns = line.split()
            x3[indx,indy] = float(columns[0])
            y3[indx,indy] = float(columns[2])
    f.close()
    g.close()

    f = open('Bl2.dat', 'r')
    header=f.readline()
    header=header.strip()
    columnsh = header.split()
    nxt=int(columnsh[0])
    nyt=int(columnsh[1])
    x2=np.zeros((nxt,nyt))
    y2=np.zeros((nxt,nyt))
    for indx in range(0,nxt):
        for indy in range(0,nyt):
            line = f.readline()
            line = line.strip()
            columns = line.split()
            x2[indx,indy] = float(columns[0])
            y2[indx,indy] = float(columns[2])
    f.close()
    x=np.zeros((nx,nyt))
    y=np.zeros((nx,nyt))
    x[0,:]=x1[:,TE-1:]
    y[0,:]=y1[:,TE-1:]
    x[1:nx-1,:]=x2
    y[1:nx-1,:]=y2
    x[nx-1,:]=x3[:,TE-1:]
    y[nx-1,:]=y3[:,TE-1:]
#-----------------------------------------------Calculate metrics
    dxdxi=d1xi_2(x[:,:],nx,nyt,bc4)
    dydxi=d1xi_2(y[:,:],nx,nyt,bc4)
    dxdeta=d1eta_2(x[:,:],nx,nyt,bc4)
    dydeta=d1eta_2(y[:,:],nx,nyt,bc4)
    return dxdxi[:,0],dydxi[:,0],dxdeta[:,0],dydeta[:,0],nyt-1

def bound4():
    import math as m
    import numpy as np

    print('calculate carpenter parameters')
    r1=-(2177*m.sqrt(295369)-1166427)/25488
    r2=(66195*m.sqrt(53)*m.sqrt(5573)-35909375)/101952

    al4=np.zeros((4,4))
    al4[0,0]=-(216* r2+2160* r1-2125) /12960
    al4[0,1]= (81*  r2+675*  r1+415)  /540
    al4[0,2]=-(72*  r2+720*  r1+445)  /1440
    al4[0,3]=-(108* r2+756*  r1+421)  /1296
    al4[1,0]= (81 * r2+675*  r1+415)  /540
    al4[1,1]=-(4104*r2+32400*r1+11225)/4320
    al4[1,2]= (1836*r2+14580*r1+7295) /2160
    al4[1,3]=-(216 *r2+2160* r1+655)  /4320
    al4[2,0]=-(72*  r2+720*  r1+445)  /1440
    al4[2,1]= (1836*r2+14580*r1+7295) /2160
    al4[2,2]=-(4104*r2+32400*r1+12785)/4320
    al4[2,3]= (  81*r2+675*  r1+335)  /540
    al4[3,0]=-(108* r2+756*  r1+421)  /1296
    al4[3,1]=-(216* r2+2160* r1+655)  /4320
    al4[3,2]= (81*  r2+ 675* r1+335)  /540
    al4[3,3]=-(216* r2+2160* r1-12085)/12960

    ar4=np.zeros((4,6))

    ar4[0,0]= (-1)/2
    ar4[0,1]=-(864*r2+6480*r1+305)/4320
    ar4[0,2]= (216*r2+1620*r1+725)/540
    ar4[0,3]=-(864*r2+6480*r1+3335)/4320

    ar4[1,0]= (864*r2+6480*r1+305)/4320
    ar4[1,1]= 0
    ar4[1,2]=-(864*r2+6480*r1+2315)/1440
    ar4[1,3]= (108*r2+810*r1+415)/270

    ar4[2,0]=-(216*r2+1620*r1+725)/540
    ar4[2,1]= (864*r2+6480*r1+2315)/1440
    ar4[2,2]= 0
    ar4[2,3]=-(864*r2+6480*r1+785)/4320

    ar4[3,0]= (864*r2+6480*r1+3335)/4320
    ar4[3,1]=-(108*r2+810*r1+415)/270
    ar4[3,2]= (864*r2+6480*r1+785)/4320
    ar4[3,3]= 0

    ar4[0,4]= 0
    ar4[1,4]= 0
    ar4[2,4]=-1/12
    ar4[3,4]=8/12

    ar4[0,5]=0
    ar4[1,5]=0
    ar4[2,5]=0
    ar4[3,5]=-1/12

    minv=np.linalg.inv(al4)
    bc4=np.zeros((4,6))
    for i in range(0,4):
        for j in range(0,6):
            asum = 0.0
            for d in range(0,4):
                asum = asum + minv[i,d]*ar4[d,j]
            bc4[i,j]=asum
    return bc4
#===========================================================
def d1xi_2(fn_,nxp,nyp,bc4):
    import math as m
    import numpy as np

    fn=np.zeros((nxp+2,nyp+2))
    fn[1:nxp+1,1:nyp+1]=fn_
    del fn_

    ds=1
    factx=ds/12
    dfn=np.zeros((nxp+2,nyp+2))
    for i in range(3,nxp-3+1):
        dfn[i,:]=(fn[i-2,:]-fn[i+2,:]+8*(fn[i+1,:]-fn[i-1,:]))*factx

    for i in range(0,4):
        dfn[i+1,:]=ds*(bc4[i,0]*fn[1,:]+bc4[i,1]*fn[2,:]+bc4[i,2]*fn[3,:]+bc4[i,3]*fn[4,:]+bc4[i,4]*fn[5,:]+bc4[i,5]*fn[6,:])

    for i in range(0,4):  
        dfn[nxp-i,:]= -ds*( bc4[0+i,0]*fn[nxp-0,:]+bc4[0+i,1]*fn[nxp-1,:]+bc4[0+i,2]*fn[nxp-2,:]+bc4[0+i,3]*fn[nxp-3,:]+bc4[0+i,4]*fn[nxp-4,:]+bc4[0+i,5]*fn[nxp-5,:] )

    dfn_=np.zeros((nxp,nyp))
    dfn_= dfn[1:nxp+1,1:nyp+1]
    return dfn_
#===========================================================
def d1eta_2(fn_,nxp,nyp,bc4):
    import math as m
    import numpy as np

    fn=np.zeros((nxp+2,nyp+2))
    fn[1:nxp+1,1:nyp+1]=fn_
    del fn_

    ds=1
    factx=ds/12
    dfn=np.zeros((nxp+2,nyp+2))
    for j in range(3,nyp-3+1):
        dfn[:,j]=(fn[:,j-2]-fn[:,j+2]+8*(fn[:,j+1]-fn[:,j-1]))*factx

    for j in range(0,4):
        dfn[:,j+1]=ds*(bc4[j,0]*fn[:,1]+bc4[j,1]*fn[:,2]+bc4[j,2]*fn[:,3]+bc4[j,3]*fn[:,4]+bc4[j,4]*fn[:,5]+bc4[j,5]*fn[:,6])

#    for j in range(0,4):
#        dfn[:,nyp-j]= -ds*( bc4[0+j,0]*fn[:,nyp-0]+bc4[0+j,1]*fn[:,nyp-1]+bc4[0+j,2]*fn[:,nyp-2]+bc4[0+j,3]*fn[:,nyp-3]+bc4[0+j,4]*fn[:,nyp-4]+bc4[0+j,5]*fn[:,nyp-5] )

    dfn_=np.zeros((nxp,nyp))
    dfn_= dfn[1:nxp+1,1:nyp+1]
    return dfn_
 

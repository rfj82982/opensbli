from opensbli.utilities.numerical_functions import spline, splint
from sympy import Piecewise
from scipy.integrate import odeint
import numpy as np
import numpy.polynomial.polynomial as poly
import matplotlib.pyplot as plt
# from opensbli.initialisation import GridBasedInitialisation
from opensbli.core.opensbliobjects import DataObject, ConstantObject
from opensbli.core.grid import GridVariable
from opensbli.core.kernel import ConstantsToDeclare as CTD
# from opensbli.core.kernel import Kernel
import warnings
# from scipy.optimize import curve_fit
# from opensbli.equation_types.opensbliequations import OpenSBLIEq

plt.style.use('classic')


class Boundary_layer_profile(object):
    """ Performs a similarity solution (Viscous fluid flow, F.White 1974),
    to obtain u and T profiles for a laminar compressible boundary-layer.

    :arg float xmach: Free-stream Mach numnber.
    :arg float Pr: Prandtl number.
    :arg float gama: Ratio of specific heats.
    :arg float Tw: Wall temperature, use -1 for adibatic wall conditions.
    :arg float Re: Free-stream Reynolds number.
    :arg float Tinf: Dimensional free-stream temperature for Sutherland's law."""

    # def __init__(self, xmach, Pr, gama, Tw, Re, Tinf):
    # def __init__(self, Re, xmach, Tinf, Tw, Sc, adiabatic_condition, catalytic_condition, cN2e, cNe, cO2e, cOe, cNOe, pref, rhoref, uref, blthickness):
    def __init__(self, Re, xmach, Tinf, Tw, Sc, adiabatic_condition, catalytic_condition, cN2e, cNe, cO2e, cOe, cNOe, pref):

        self.y, self.u, self.T, self.scale = self.generate_boundary_layer_profile(Re, xmach, Tinf, Tw, Sc, adiabatic_condition, catalytic_condition, cN2e, cNe, cO2e, cOe, cNOe, pref)
        # self.uref, self.tref, self.rhoref, self.blthickness 
        self.uref, self.tinf, self.rhoref, self.blthickness, self.rho = self.extract_parameters()
        self.rhoref, self.cN2, self.cN, self.cO2, self.cO, self.cNO = self.extract_densities()
        
        self.Re = Re
        self.n = np.size(self.y)
        return

    def compbl(self, v, p=None):
        """ Sets up the system of equations to be integrated by odeint.

        :arg ndarray: v: Solution vector.
        :arg None: p: Empty dummy argument required by the odeint function.
        :returns: list: dv: System of equations."""
        suth = self.suth
        c = np.sqrt(v[3])*(1.0+suth)/(v[3]+suth)
        dcdg = 1.0/(2.0*np.sqrt(v[3])) - np.sqrt(v[3])/(v[3]+suth)
        dcdg *= (1.0+suth) / (v[3]+suth)
        cp = dcdg*v[4]
        dv = [v[1], v[2], -v[2]*(cp+v[0])/c, v[4],
              -v[4]*(cp+self.pr*v[0])/c - self.pr*(self.gama-1)*self.xmach**2 * v[2]**2]
        return dv
    
    def mucalc(self, T,cN2,cN,cO2,cO,cNO):
        mu=((cN2+cO2+cNO)*0.1*np.exp(-11.2202)*T**(0.021823*np.log(T)+0.343567)+(cN+cO)*0.1*np.exp(-11.7344)*T**(0.022652*np.log(T)+0.42509))*(1.0-np.exp(-0.010568*T))
        kappa=(1410.0*(cN2+cO2+cNO)*0.1*np.exp(-11.2202)*T**(0.021823*np.log(T)+0.343567)+2210.0*(cN+cO)*0.1*np.exp(-11.7344)*T**(0.022652*np.log(T)+0.42509))*(1.0-np.exp(-0.010568*T))   
        return mu, kappa 

    # frozen specific heat (thermal equilibrium)
    def cpcalc(self, T,cN2,cN,cO2,cO,cNO,MN2,MN,MO2,MO,MNO): 
        Rhat=8314.3
        thetavO2, thetavN2, thetavNO = 2270.0, 3390.0, 2740.0

        cpN=5.0/2.0*Rhat/MN
        cpO=5.0/2.0*Rhat/MO

        cpN2=(7.0/2.0+((thetavN2/T)**2*np.exp(thetavN2/T)/(np.exp(thetavN2/T)-1.0)**2))*Rhat/MN2
        cpO2=(7.0/2.0+((thetavO2/T)**2*np.exp(thetavO2/T)/(np.exp(thetavO2/T)-1.0)**2))*Rhat/MO2
        cpNO=(7.0/2.0+((thetavNO/T)**2*np.exp(thetavNO/T)/(np.exp(thetavNO/T)-1.0)**2))*Rhat/MNO

        cp=cN2*cpN2+cN*cpN+cO2*cpO2+cO*cpO+cNO*cpNO
        Mhat=1.0/(cN2/MN2+cN/MN+cO2/MO2+cO/MO+cNO/MNO); # average molar mass
        af=np.sqrt(cp*Rhat*T/(Mhat*cp-Rhat)); # frozen flow specific heat

        return cp,cpN2,cpN,cpO2,cpO,cpNO,af

    def hcalc(self, T,cN2,cN,cO2,cO,cNO,MN2,MN,MO2,MO,MNO):
        Rhat=8314.3
        thetavO2, thetavN2, thetavNO = 2270.0, 3390.0, 2740.0
        dhO, dhN, dhNO = 59.544, 112.951, 21.6009

        hN=5.0/2.0*Rhat/MN*T
        hO=5.0/2.0*Rhat/MO*T
        hN2=7.0/2.0*Rhat/MN2*T+thetavN2*Rhat/(MN2*(np.exp(thetavN2/T)-1.0))
        hO2=7.0/2.0*Rhat/MO2*T+thetavO2*Rhat/(MO2*(np.exp(thetavO2/T)-1.0))
        hNO=7.0/2.0*Rhat/MNO*T+thetavNO*Rhat/(MNO*(np.exp(thetavNO/T)-1.0))
        h=cN2*hN2+cN*hN+cO2*hO2+cO*hO+cNO*hNO
        hform=4.1868e6*(dhO*cO/MO+dhN*cN/MN+dhNO*cNO/MNO) # enthalpy for formation (including unit conversion)

        return h, hform

    def derivmatrices(self, d1,d2,d3,h,n):

        d1[0,0]=-1.5/h; d1[0,1]=2.0/h; d1[0,2]=-0.5/h
        d1[n-1,n-3]=0.5/h; d1[n-1,n-2]=-2.0/h; d1[n-1,n-1]=1.5/h

        d2[0,0]=1.0/h**2; d2[0,1]=-2.0/h**2; d2[0,2]=1.0/h**2
        d2[n-1,n-3]=1.0/h**2; d2[n-1,n-2]=-2.0/h**2; d2[n-1,n-1]=1.0/h**2

        for i in range(1,n-1):
            d1[i,i-1]=-1.0/(2.0*h)
            d1[i,i+1]=1.0/(2.0*h)

            d2[i,i-1]=1.0/h**2
            d2[i,i]=-2.0/h**2
            d2[i,i+1]=1.0/h**2

        d3[0,0]=1.0; d3[1,1]=1.0; d3[n-2,n-2]=1.0; d3[n-1,n-1]=1.0; # the first and last 2 rows will anyway be overwritten later

        for i in range(2,n-2):
            d3[i,i-2]=-1.0/(2.0*h**3)
            d3[i,i-1]=2.0/(2.0*h**3)
            d3[i,i+1]=-2.0/(2.0*h**3)
            d3[i,i+2]=1.0/(2.0*h**3)

        return d1, d2, d3

    def generate_boundary_layer_profile(self, Re, xmach, Tinf, Tw, Sc, adiabatic_condition, catalytic_condition, cN2e, cNe, cO2e, cOe, cNOe, pref):
        """ Generates a boundary layer initial profile. Solves the mean flow
        in a compressible boundary layer. (Equations 7.32) in White (1974).

        :arg float xmach: Mach number.
        :arg float pr: Prandtl number.
        :arg float gama: Ratio of specific heats.
        :arg float Tw: Wall temperature Tw/Tinf (< 0 for adiabatic)
        :arg float Tinf: Freestream reference temperature (Kelvin)
        :arg float Re: Reynolds number."""

        # molar masses
        MN2, MN, MO2, MO, MNO = 28.0, 14.0, 32.0, 16.0, 30.0 # molar mass
        Rhat = 8314.3 

        # parameters
        etamax, n, niter = 10.0, 1001, 100
        n_step = n - 1

        # flow concentrations
        # cN2e, cNe, cO2e, cOe, cNOe = 0.64, 0.05, 0.06, 0.15, 0.1 # free stream mass fractions (note that c_s are mass fractions here)
        # cN2e, cNe, cO2e, cOe, cNOe = 1.0, 0.0, 0.0, 0.0, 0.0 # free stream mass fractions (note that c_s are mass fractions here)

        self.xmach, self.Re, self.Tw = xmach, Re, Tw

        # freestream mach number, pressure and temperature along with wall temperature and SChmidt number (currently assumed constant for all species)
        Me, pref, Te, Tw, Sc = xmach, pref, Tinf, Tw, Sc

        # wall properties 
        adiabatic, catalytic = adiabatic_condition, catalytic_condition 

        # derived quantities in free stream
        # sigOtoNe = (2.0*cO2e/MO2+cOe/MO+cNOe/MNO)/(2.0*cN2e/MN2+cNe/MN+cNOe/MNO) # ratio of moles of O to N in freestream
        sigOtoNe = 0.0
        mue, kappae = self.mucalc(Te,cN2e,cNe,cO2e,cOe,cNOe)

        # kappave=kappavcalc(Te,cN2e,cNe,cO2e,cOe,cNOe); # not used in this
        # formulation
        cpe, dum1, dum2, dum3, dum4, dum5, afe = self.cpcalc(Te,cN2e,cNe,cO2e,cOe,cNOe,MN2,MN,MO2,MO,MNO)
        Pre=cpe*mue/kappae # not used

        he, hforme = self.hcalc(Te,cN2e,cNe,cO2e,cOe,cNOe,MN2,MN,MO2,MO,MNO)
        Ue=afe*Me # work out a Ue based on the frozen sound speed

        # initialise arrays
        A=np.zeros([n,n]); B=np.zeros([n,n]); C=np.zeros([n,n]); eta=np.zeros([n,1]); y=np.zeros([n,1]); d1=np.zeros([n,n]); d2=np.zeros([n,n]); d3=np.zeros([n,n])
        lvar=np.zeros([n,1]); mvar=np.zeros([n,1]); nvar=np.zeros([n,1]); cvar=np.zeros([n,1]); dl=np.zeros([n,1]); dm=np.zeros([n,1]); rhs=np.zeros([n,1])
        f=np.zeros([n,1]); df=np.zeros([n,1]); d2f=np.zeros([n,1]); theta=np.zeros([n,1]); dtheta=np.zeros([n,1]); rho=np.zeros([n,1])
        cN2=np.zeros([n,1]); cN=np.zeros([n,1]); cO2=np.zeros([n,1]); cO=np.zeros([n,1]); cNO=np.zeros([n,1])
        dcN2=np.zeros([n,1]); dcN=np.zeros([n,1]); dcO2=np.zeros([n,1]); dcO=np.zeros([n,1]); dcNO=np.zeros([n,1])

        # grid size (uniform grid) and set derivative matrices
        deta=etamax/(n-1)
        d1,d2,d3=self.derivmatrices(d1,d2,d3,deta,n)

        

        # set initial variation of f,theta and cs
        for i in range(0,n):
            eta[i]=(i)*deta
            f[i]=eta[i] # simple estimate works fine
            theta[i]=1.0
            cN2[i]=cN2e
            cN[i]=cNe
            cO2[i]=cO2e
            cO[i]=cOe
            cNO[i]=cNOe

        df=d1@f

        # set boundary conditions for momentum equation
        A[0,0]=1.0 # first row = stream function boundary condition
        A[1,0]=-1.5/deta; A[1,1]=2.0/deta; A[1,2]=-0.5/deta; # second row = velocity no slip condition
        A[n-1,n-1]=1.5/deta; A[n-1,n-2]=-2.0/deta; A[n-1,n-3]=0.5/deta;# nth row = external velocity boundary condition
        A[n-2,n-1]=1.0/deta**2; A[n-2,n-2]=-2.0/deta**2; A[n-2,n-3]=1.0/deta**2;# n-1 row = zero derivative boundary condition (is this needed?)

        #set Dirichlet bc in matrix B
        B[0,0]=1.0; # Dirichlet (T=Tw)
        B[n-1,n-1]=1.0; # Dirichlet (T=Te)

        # set Neumann wall bc in matrix C
        C[0,0]=-1.5/deta; C[0,1]=2.0/deta; C[0,2]=-0.5/deta
        C[n-1,n-1]=1.0

        
        tol=1.e-6; Hklast=0.0 # we will iterate on f until Hk converges to tol

        # main iteration loop
        for iter in range(0, niter):
            # print(iter)
            # set local variables and required derivatives
            # dcN2=d1*cN2; dcN=d1*cN; dcO2=d1*cO2; dcO=d1*cO; dcNO=d1*cNO

            dcN2, dcN, dcO2, dcO, dcNO = np.matmul(d1,cN2), d1 @ cN, d1 @ cO2, d1 @ cO, d1 @ cNO

            for i in range(0,n):
                T = theta[i]*Te
                mu,kappa = self.mucalc(T,cN2[i],cN[i],cO2[i],cO[i],cNO[i])
                cp,cpN2,cpN,cpO2,cpO,cpNO,dum = self.cpcalc(T,cN2[i],cN[i],cO2[i],cO[i],cNO[i],MN2,MN,MO2,MO,MNO)
                lvar[i]=mu/(mue*theta[i]) # Chapman-Rubesin l=(rho.mu)/(rhoe.mue)
                cvar[i]=cp/cpe
                mvar[i]=kappa/(theta[i]*cpe*mue)
                nvar[i]=lvar[i]/(Sc*cpe)*(cpN2*dcN2[i] + cpN*dcN[i] + cpO2*dcO2[i] + cpO*dcO[i] + cpNO*dcNO[i])

            dl=d1 @lvar
            dm=d1 @mvar

            # first solve the momentum equation
            for i in range(2,n-2):
                for j in range(0,n):
                    A[i,j]=lvar[i]*d3[i,j]+(dl[i]+f[i])*d2[i,j]

            rhs[0:n-2]=0.0; # wall bc
            rhs[n-1]=1.0 # free stream bc (i.e. f'=1)
            f= np.linalg.lstsq(A, rhs, rcond=None)[0] # solve system for f

            # next solve the energy equation
            d2f = d2 @f 


            for i in range(1,n-1): # exclude first and last rows where bc are applied
                for j in range(0,n):
                    B[i,j]=mvar[i]*d2[i,j]+(dm[i]+cvar[i]*f[i]+nvar[i])*d1[i,j]
                    C[i,j]=B[i,j] 

            rhs[1:n-2]=-lvar[1:n-2]*Ue**2/(cpe*Te)*d2f[1:n-2]**2
            rhs[n-1]=1.0

            if adiabatic:
                rhs[0] = 0.0
                theta = np.linalg.lstsq(C,rhs, rcond=None)[0] # solve system for theta
            else:
                rhs[0] = Tw/Te
                theta = np.linalg.lstsq(B,rhs, rcond=None)[0] # solve system for theta

            # solve atom mass fraction equations
            for i in range(1,n-1):
                for j in range(0,n):
                    B[i,j]=lvar[i]/Sc*d2[i,j]+(dl[i]/Sc+f[i])*d1[i,j]
                    C[i,j]=B[i,j]

            if catalytic:
                sigOtoNe = (2.0*cO2e/MO2+cOe/MO+cNOe/MNO)/(2.0*cN2e/MN2+cNe/MN+cNOe/MNO) # ratio of moles of O to N in freestream
                
                rhs[0:n]=0.0
                rhs[n-1]=cNe
                cN=np.linalg.lstsq(B,rhs, rcond=None)[0] # cNw=0.0

                rhs[n-1]=cOe
                cO=np.linalg.lstsq(B,rhs, rcond=None)[0] #cOw=0.0

                rhs[0:n-1]=0.0
                rhs[n-1]=cNOe
                cNO=np.linalg.lstsq(C,rhs, rcond=None)[0] # Neumann condition on cNO 

                cN2w=(1.0-cNO[0]+cNO[0]/2.0*MO2/MNO*(1.0-sigOtoNe))/(1.0+sigOtoNe*MO2/MN2) # using mole ratio O:N
                #cN2w=cN2e+cNe; # simpler (only for frozen flow cases?)
                rhs[0]=cN2w
                rhs[n-1]=cN2e
                cN2=np.linalg.lstsq(B,rhs, rcond=None)[0]

                rhs[0]=1.0-cN2w-cNO[0]
                rhs[n-1]=cO2e
                cO2=np.linalg.lstsq(B,rhs, rcond=None)[0]

            else:
                rhs[1:n-1]=0.0
                rhs[n-1]=cN2e; cN2=np.linalg.lstsq(C,rhs, rcond=None)[0] 
                rhs[n-1]=cNe;  cN=np.linalg.lstsq(C,rhs, rcond=None)[0] 
                rhs[n-1]=cO2e; cO2=np.linalg.lstsq(C,rhs, rcond=None)[0] 
                rhs[n-1]=cOe;  cO=np.linalg.lstsq(C,rhs, rcond=None)[0] 
                rhs[n-1]=cNOe; cNO=np.linalg.lstsq(C,rhs, rcond=None)[0] 
    
            # monitor shape factor (in eta space) for convergence
            df=d1@f
            # delta1=abs(np.trapz(eta,(1.0-df))) # kinematic displacement thickness
            delta1=np.trapz( 1.0 -  df.T , eta.T)
            delta2=np.trapz( df.T*(1.0 -  df.T) , eta.T)
            # delta2=abs(np.trapz(eta,df*(1.0-df))) # kinematic momemtum thickness
            Hk=delta1/delta2; # kinematic shape factor
            err=abs(Hk-Hklast)
            # print('shape factor', Hk)

            # print('simulation error', err)
            if err < tol:
                break 

            Hklast=Hk

        # velocity value
        df = d1@f
        rho = ((cN2e/MN2+cNe/MN+cO2e/MO2+cOe/MO+cNOe/MNO)/(cN2/MN2+cN/MN+cO2/MO2+cO/MO+cNO/MNO))/theta #rho/rhoe
        
        df, theta = np.reshape(df, np.size(df)), np.reshape(theta, np.size(theta))
        cN2, cN, cO2, cO, cNO = np.reshape(cN2, np.size(cN2)), np.reshape(cN, np.size(cN)), np.reshape(cO2, np.size(cO2)), np.reshape(cO, np.size(cO)), np.reshape(cNO, np.size(cNO))
        rho = np.reshape(rho, np.size(rho))

        self.rho = rho 
    
        delta1 = np.trapz( (1.0-df.T*rho.T) , y.T)
        delta2 = np.trapz( df.T*rho.T*(1.0-df.T) , y.T)

        rhoref = (pref/(Rhat*Te)) * (cN2e/MN2 +cNe/MN +cO2e/MO2 +cOe/MO +cNOe/MNO)**(-1)
        rhof = (cN2/MN2+cN/MN+cO2/MO2+cO/MO+cNO/MNO)/theta #rho/rhoe

        print('freestream velocity, uref :', Ue)
        print('freestream density, rhoref :', rhoref)

        self.Twall = theta[0]
        self.Twall = self.Twall
        print("The wall temperature is :", self.Twall)

        self.deta = deta 
        self.df, self.theta = df, theta
        Rex = 950 
        deltastar = Rex*mue/(Ue*rhoref)
        print('inlet boundary layer thickness :', deltastar)

        

        self.ue, self.tinf, self.rhoref, self.blthickness = Ue, Te, rhoref, deltastar 
        # yO=rhoO/(MO*sumy)
        # sumy=rhoO/MO+rhoO2/MO2+rhoN/MN+rhoN2/MN2+rhoNO/MNO
        self.cN2, self.cN, self.cO2, self.cO, self.cNO = cN2, cN, cO2, cO, cNO

       
        

        # blthickness, uref, tinf, rhoref = 
        print('input variables line one: ', Re, Me, self.tinf, self.Twall, Sc)
        print('input variables: ', pref, rhoref, self.ue, deltastar)
        
        # define constants
        Twall_cpp = ConstantObject('Twall')
        Twall_cpp.value = self.Twall
        CTD.add_constant(Twall_cpp)

        # define constants
        pref_cpp, rhoref_cpp, uref_cpp, blthickness_cpp = ConstantObject('pref'), ConstantObject('rhoref'), ConstantObject('uref'), ConstantObject('blthicknesss')
        pref_cpp.value, rhoref_cpp.value, uref_cpp.value, blthickness_cpp.value = pref, rhoref, self.ue, deltastar

        CTD.add_constant(pref_cpp); CTD.add_constant(rhoref_cpp); CTD.add_constant(uref_cpp); CTD.add_constant(blthickness_cpp)
        

        # self.Twall = v[1]
        # print("The wall temperature is :", self.Twall)
        # y, u, T, scale_factor = self.integrate_boundary_layer(nstep)
        y, u, T, scale_factor = self.integrate_boundary_layer(n_step)

        self.scale_factor = scale_factor
        print("The scale factor is :", self.scale_factor)
        # print("Wall normal derivative of velocity at the wall is :", self.dudy)
        return y, u, T, scale_factor
    
    def extract_parameters(self):

        # reference velocity
        uref = self.ue
        rhoref = self.rhoref
        tinf = self.tinf

        blthickness = self.blthickness
        rho = self.rho 

        return uref, tinf, rhoref, blthickness, rho
    
    def extract_densities(self):

        # reference velocity
        # rhoN2, rhoN, rhoO2, rhoO, rhoNO = self.cN2*self.rho, self.cN*self.rho, self.cO2*self.rho, self.cO*self.rho, self.cNO*self.rho
        cN2, cN, cO2, cO, cNO = self.cN2, self.cN, self.cO2, self.cO, self.cNO
        
        return self.rhoref, cN2, cN, cO2, cO, cNO
        
    def integrate_boundary_layer(self, n):
        """ Integrates the boundary-layer and calculates the scale factor from displacement thickness.

        :arg int n: Iteration number from the iterative solver.
        :returns: ndarray: y: Wall normal coordinates.
        :returns: ndarray: u: Streamwise velocity component profile.
        :returns: ndarray: T: Temperature profile.
        :returns: float: scale: Scale factor of the boundary-layer."""
        sumd, record_z = 0, 0
        z = np.zeros(n+1)
        # d_eta = self.eta[1]*0.5
        d_eta = self.deta*0.5

        # self.soln[1,:] is the u velocity, should be 1 in free stream
        for i in range(1, n+1):
            z[i] = z[i-1] + d_eta*(self.theta[i] + self.theta[i-1])
            dm = self.theta[i-1]- self.df[i-1]
            dd = self.theta[i] - self.df[i]
            sumd += d_eta*(dd+dm)
            if(self.df[i] > 0.999 and record_z < 1.0):
                # print "recording at iteration: ", i
                # dlta = z[i]
                record_z = 2.0
            scale = sumd
        # print("delta is :", dlta)
        print("conversion factor is: ", scale)
        # print("scaled delta is: ", dlta/scale)
        # Rescale with displacement thickness and convert to FLOWER variable normalisation
        y, u, T = z/scale, self.df[:], self.theta[:]
        # Calculate du/dy at the wall
        dy = y[1]
        # self.dudy = (-3*u[0]+4*u[1]-u[2])/(2.0*dy)
        self.dudy = (-1.83333333333334*u[0]+3.00000000000002*u[1]-1.50000000000003*u[2]+0.333333333333356*u[3]-8.34657956545823e-15*u[4]+1.06910315192207e-15*u[5])/dy
        self.dTdy = (-1.83333333333334*T[0]+3.00000000000002*T[1]-1.50000000000003*T[2]+0.333333333333356*T[3]-8.34657956545823e-15*T[4]+1.06910315192207e-15*T[5])/dy
        return y, u, T, scale

from opensbli.initialisation import GridBasedInitialisation
class Initialise_Flatmix(GridBasedInitialisation):
    """ Generates the initialiastion equations for the boundary-layer profile.

    :arg list npoints: Numerical values of the number of points in each direction.
    :arg list lengths: Numerical values of the problem dimensions.
    :arg list directions: Integer values of the problem directions.
    :arg list betas: Stretching factors for stretched grids.
    :arg int n_coeffs: Desired number of coefficients for the polynomial fit.
    :arg float Re: Reynolds number.
    :arg float xMach: Free-stream Mach number"""
    def __new__(cls, bl_directions, n_coeffs, Re, xMach, Tinf, Twall, Sc, adiabatic_condition, catalytic_condition, cN2e, cNe, cO2e, cOe, cNOe, pref, rhoref, uref, blthickness, coordinate_evaluations=None):
        ret = super(Initialise_Flatmix, cls).__new__(cls)
        print("Polynomial boundary-layer initialiastion called with Re = %f, Mach = %f, T_inf = %f." % (Re, xMach, Tinf))
        print("                                                     cN2e = %f, cNe = %f, cO2e = %f, cOe = %f, cNOe = %f." % (cN2e, cNe, cO2e, cOe, cNOe))
        # print("                    with given dimensional constants pref = %f, rhoref = %f, uref = %f, blthickness = %f." % (pref, rhoref, uref, blthickness,))

        ret.coordinates = [x[1] for x in bl_directions]
        ret.bl_directions = bl_directions
        ret.n_coeffs = n_coeffs
        ret.coordinate_evaluations = coordinate_evaluations
        ret.Re = ret.find_constant_values([Re])[0]
        ret.Tinf = ret.find_constant_values([Tinf])[0]
        ret.Tw = ret.find_constant_values([Twall])[0]
        ret.equations = []
        ret.xMach = ret.find_constant_values([xMach])[0]

        # get specified adabatic and catalytic conditions
        ret.adiabatic_condition = ret.find_constant_values([adiabatic_condition])[0]
        ret.catalytic_condition = ret.find_constant_values([catalytic_condition])[0]

        # get concentrations of species
        ret.cN2e = ret.find_constant_values([cN2e])[0]
        ret.cNe = ret.find_constant_values([cNe])[0]
        ret.cO2e = ret.find_constant_values([cO2e])[0]
        ret.cOe = ret.find_constant_values([cOe])[0]
        ret.cNOe = ret.find_constant_values([cNOe])[0]

        # get dimensional constants
        ret.Sc = ret.find_constant_values([Sc])[0]
        ret.pref = ret.find_constant_values([pref])[0]
        # ret.rhoref = ret.find_constant_values([rhoref])[0]
        # ret.uref = ret.find_constant_values([uref])[0]
        # ret.blthickness = ret.find_constant_values([blthickness])[0]

        return ret

    def find_constant_values(self, input):
        outlist = []
        for l in input:
            if isinstance(l, ConstantObject):
                if isinstance(l.value, str):
                    raise ValueError("")
                else:
                    outlist += [l.value]
            else:
                outlist += [l]
        return outlist

    def check_inputs(self, block):
        bl_directions = [x[0] for x in self.bl_directions]
        if sum(bl_directions) < 1:
            raise ValueError("Provide the directions to apply a boundary layer profile in.")
        if len(bl_directions) != block.ndim:
            raise ValueError("The list of polynomial directions must match the dimensions of the problem.")
        if self.n_coeffs < 10:
            raise ValueError("Higher number of polynomial coefficients are required for a good polynomial fit.")
        return

    def spatial_discretisation(self, block):
        self.equations = []
        self.block = block
        self.idxs = block.grid_indexes
        self.check_inputs(block)

        # Check if user has passed equations to evaluate coordinates, and add them to the kernel
        if self.coordinate_evaluations:
            self.equations += self.coordinate_evaluations
        self.initial = self.generate_initial_condition()
        # Add polynomial equations to initialise the solution
        self.equations += self.eqns

        self.equations = block.dataobjects_to_datasets_on_block(self.equations)
        # print(self.order)
        kernel1 = Kernel(block, computation_name="Grid_based_initialisation%d" % self.order)
        kernel1.set_grid_range(block)
        schemes = block.discretisation_schemes
        for d in range(block.ndim):
            for sc in schemes:
                if schemes[sc].schemetype == "Spatial":
                    kernel1.set_halo_range(d, 0, schemes[sc].halotype)
                    kernel1.set_halo_range(d, 1, schemes[sc].halotype)
        kernel1.add_equation(self.equations)
        kernel1.update_block_datasets(block)
        self.Kernels = [kernel1]
        return

    def generate_initial_condition(self):
        n_coeffs = self.n_coeffs
        # Load from similarity solution class
        # y, u, T, rho, n, uref, Tinf, rhoref, blthickness, rhoN2, rhoN, rhoO2, rhoO, rhoNO = self.load_similarity()
        y, u, T, rho, n, uref, Tinf, rhoref, blthickness, cN2, cN, cO2, cO, cNO = self.load_similarity()
        self.blthickness = blthickness
        
        y0 = y 
        # rho = rho[:]
        # y0, u0, T0, rho0, n0 = self.load_similarity()
        MN2, MN, MO2, MO, MNO = 28.0, 14.0, 32.0, 16.0, 30.0 # molar mass
        Rhat = 8314.3
        # rho = (1.0/(( cN2/MN2 + cN/MN + cO2/MO2 + cO/MO + cNO/MNO)))*(1.0/T)
        rho = (self.pref/(Rhat*( self.cN2e/MN2 + self.cNe/MN + self.cO2e/MO2 + self.cOe/MO + self.cNOe/MNO)))*(1.0/T)
        v0, dvdy = self.solved_continuity(y, u, rho)
        v0 = v0/rho

        # uref, Tinf, rhoref, blthickness = self.load_similarity_reference_variables()
        # rhoN2, rhoN, rhoO2, rhoO, rhoNO = self.load_similarity_densities()

        y, u, v, T, rho, n = y*blthickness , u*uref, v0*uref, T*Tinf, rho*rhoref, n # dimensionalise all variables - gnsa1e21, 2022
        
        # rho = (self.pref/(Rhat*( cN2/MN2 + cN/MN + cO2/MO2 + cO/MO + cNO/MNO)))*(1.0/T)

        print('edge concentrations of species')
        print('cN2', cN2[-1])
        print('cN', cN[-1])
        print('cO2', cO2[-1])
        print('cO', cO[-1])
        print('cNO', cNO[-1])

        # rhoN2, rhoN, rhoO2, rhoO, rhoNO = rhoN2*rhoref, rhoN*rhoref, rhoO2*rhoref, rhoO*rhoref, rhoNO*rhoref
        # rhoN2, rhoN, rhoO2, rhoO, rhoNO = cN2*rho, cN*rho, cO2*rho, cO*rho, cNO*rho
        self.dudy, dvdy = self.dudy*uref/blthickness, dvdy*uref/blthickness # dimensionalise all variables - gnsa1e21, 2022

        # fig3,ax3=plt.subplots(2, 3)

        # ax3[0,0].plot(cN2, y)
        # ax3[0,1].plot(cN, y)
        # ax3[0,2].plot(cO2, y)
        # ax3[1,0].plot(cO, y)
        # ax3[1,1].plot(cNO, y)



        # plt.show()

        dy0 = y0[1]
        dcN2dy = (-1.83333333333334*cN2[0]+3.00000000000002*cN2[1]-1.50000000000003*cN2[2]+0.333333333333356*cN2[3]-8.34657956545823e-15*cN2[4]+1.06910315192207e-15*cN2[5])/dy0
        dcNdy = (-1.83333333333334*cN[0]+3.00000000000002*cN[1]-1.50000000000003*cN[2]+0.333333333333356*cN[3]-8.34657956545823e-15*cN[4]+1.06910315192207e-15*cN[5])/dy0
        dcO2dy = (-1.83333333333334*cO2[0]+3.00000000000002*cO2[1]-1.50000000000003*cO2[2]+0.333333333333356*cO2[3]-8.34657956545823e-15*cO2[4]+1.06910315192207e-15*cO2[5])/dy0
        dcOdy = (-1.83333333333334*cO[0]+3.00000000000002*cO[1]-1.50000000000003*cO[2]+0.333333333333356*cO[3]-8.34657956545823e-15*cO[4]+1.06910315192207e-15*cO[5])/dy0
        dcNOdy = (-1.83333333333334*cNO[0]+3.00000000000002*cNO[1]-1.50000000000003*cNO[2]+0.333333333333356*cNO[3]-8.34657956545823e-15*cNO[4]+1.06910315192207e-15*cNO[5])/dy0

        # Tolerance for finding the edge of the boundary layer.
        tolerance = 1e-10
        bl_directions = [x[0] for x in self.bl_directions]
        

        if sum(bl_directions) == 1:  # 2D Flatmix and 3D spanwise periodic Flatmix, boundary layer in one direction
            # Create a large array of coordinates for this direction to interpolate the profile onto
            dire = [i for i, x in enumerate(self.bl_directions) if x[0]][0]
            # poly_coordinates = self.uniform_1d_coordinate()
            poly_coordinates = self.uniform_1d_coordinate_dimensional()
            # Interpolate u, T, rho onto the grid
            u_new = self.interpolate_onto_grid(y, poly_coordinates, u, self.dudy, 0)
            T_new = self.interpolate_onto_grid(y, poly_coordinates, T, 0, 0)
            v_new = self.interpolate_onto_grid(y, poly_coordinates, v, dvdy, 0)

            # rho_new = self.interpolate_onto_grid(y, poly_coordinates, rho, 0, 0)

            cN2_new = self.interpolate_onto_grid(y, poly_coordinates, cN2, dcN2dy, 0)
            cN_new = self.interpolate_onto_grid(y, poly_coordinates, cN, dcNdy, 0)
            cO2_new = self.interpolate_onto_grid(y, poly_coordinates, cO2, dcO2dy, 0)
            cO_new = self.interpolate_onto_grid(y, poly_coordinates, cO, dcOdy, 0)
            cNO_new = self.interpolate_onto_grid(y, poly_coordinates, cNO, 0, 0)

            

            rho_new = (self.pref/(Rhat*(cN2_new/MN2 + cN_new/MN + cO2_new/MO2 + cO_new/MO + cNO_new/MNO)))*(1.0/T_new) # create array of varying rhoref values - gnsa1e21, 2023
            
            # fig3,ax3=plt.subplots(2, 3)

            # ax3[0,0].plot(cN2_new, poly_coordinates)
            # ax3[0,1].plot(cN_new, poly_coordinates)
            # ax3[0,2].plot(cO2_new, poly_coordinates)
            # ax3[1,0].plot(cO_new, poly_coordinates)
            # ax3[1,1].plot(cNO_new, poly_coordinates)
            # ax3[1,2].plot(rho_new, poly_coordinates)

            

            # plt.show()

            rhou_new = rho_new*u_new
            rhov_new = rho_new*v_new

            # Solve continuity equation to obtain rhov
            # rhov_new = self.solve_continuity(poly_coordinates, u_new, rho_new)
            print('velocity profile, u0', u_new)

            edge = self.find_edge_of_bl(u_new, tolerance)
            # Obtain polynomial fit coefficients
            rhou_coeffs = self.fit_polynomial(poly_coordinates, rhou_new, edge, n_coeffs)
            rhov_coeffs = self.fit_polynomial(poly_coordinates, rhov_new, edge, n_coeffs)
            T_coeffs = self.fit_polynomial(poly_coordinates, T_new, edge, n_coeffs)

            cN2_coeffs = self.fit_polynomial(poly_coordinates, cN2_new, edge, n_coeffs)
            cN_coeffs = self.fit_polynomial(poly_coordinates, cN_new, edge, n_coeffs)
            cO2_coeffs = self.fit_polynomial(poly_coordinates, cO2_new, edge, n_coeffs)
            cO_coeffs = self.fit_polynomial(poly_coordinates, cO_new, edge, n_coeffs)
            cNO_coeffs = self.fit_polynomial(poly_coordinates, cNO_new, edge, n_coeffs)

            rho_coeffs = self.fit_polynomial(poly_coordinates, rho_new, edge, n_coeffs)

            self.generate_one_wall_equations([rhou_new, rhov_new, T_new, cN2_new, cN_new, cO2_new, cO_new, cNO_new, rho_new], [rhou_coeffs, rhov_coeffs, T_coeffs, cN2_coeffs, cN_coeffs, cO2_coeffs, cO_coeffs, cNO_coeffs, rho_coeffs], dire, edge, poly_coordinates)

        elif sum(bl_directions) == 2:  # 3D with one side wall in x2 # WARNING: This currently does the same thing twice,
            edges, coeffs, profiles, normal_coeffs, normal_profiles = [], [], [], [], []
            directions = [i for i, x in enumerate(self.bl_directions) if x[0]]
            for dire in directions:  # For different bl thicknesses on different walls this still needs to be a loop
                poly_coordinates = self.uniform_1d_coordinate()
                # Interpolate u, and T onto the grid
                u_new = self.interpolate_onto_grid(y, poly_coordinates, u, self.dudy, 0)
                T_new = self.interpolate_onto_grid(y, poly_coordinates, T, 0, 0)

                
                # Temperature scaling function in region [0, 1]
                Tw, Tinf = self.Twall, 1.0  # Free-stream normalised temperature of 1.0, Wall temp taken from similarity solution
                g = self.temperature_scaling(T_new, Tw, Tinf)
                rho_new = 1.0/T_new
                rhou_new = rho_new*u_new
                # Solve continuity equation to obtain rho*wall_normal_velocity
                rho_vel_normal = self.solve_continuity(poly_coordinates, u_new, rho_new)
                profiles.append([rhou_new, T_new])
                normal_profiles.append(rho_vel_normal)
                edge = self.find_edge_of_bl(u_new, tolerance)
                edges.append(edge)
                # Obtain polynomial fit coefficients
                rhou_coeffs = self.fit_polynomial(poly_coordinates, rhou_new, edge, n_coeffs)
                g_coeffs = self.fit_polynomial(poly_coordinates, g, edge, n_coeffs)
                rho_vel_normal_coeffs = self.fit_polynomial(poly_coordinates, rho_vel_normal, edge, n_coeffs)
                coeffs.append([rhou_coeffs, g_coeffs])
                normal_coeffs.append(rho_vel_normal_coeffs)
            self.generate_two_wall_equations(profiles, coeffs, directions, edges, normal_profiles, normal_coeffs, poly_coordinates)
        else:
            raise NotImplementedError("Boundary layer initialisation is not implemented for walls in 3 dimensions.")
        return

    def uniform_1d_coordinate(self):
        n_elem = 10000
        return np.linspace(0, 20.0, n_elem)
    
    def uniform_1d_coordinate_dimensional(self):
        n_elem = 10000
        return np.linspace(0, 20.0*10**self.find_exp(self.blthickness), n_elem)

    def temperature_scaling(self, temp_profile, Tw, Tinf):
        """ Computes the temperature profile between [0,1]

        :arg ndarray temp_profile: Temperature profile values between wall temperature and freestream.
        :arg float Tw: Wall temperature.
        :arg float Tinf: Free-stream temperature.
        :returns: ndarray: g: Temperature profile ranging from 0 to 1. """
        g = (temp_profile - Tw)/(Tinf - Tw)
        return g

    def form_equation(self, variable, name, coefficients, direction, edge, poly_coordinates):
        """ Creates the piecewise equations for the cases of 2D and 3D span-periodic boundary-layer profiles.

        :arg ndarray variable: Array of values for a given flow variable, used to obtain the free-stream value.
        :arg string name: Name of the variable.
        :arg ndarray coefficients: Coefficients for the polynomial fit.
        :arg int direction: Spatial direction to apply the equation to.
        :arg int edge: Grid index for the edge of the boundary-layer.
        returns: Eq: eqn: OpenSBLI equation to add to the initialisation kernel."""
        bl_edge_coordinate = poly_coordinates[edge]
        powers = [i for i in range(np.size(coefficients))][::-1]
        eqn = sum([coeff*self.coordinates[direction]**power for (coeff, power) in zip(coefficients, powers)])  # TODO set to exactl 1.0 if required
        eqn = OpenSBLIEq(GridVariable('%s' % name), Piecewise((eqn, self.coordinates[direction] < bl_edge_coordinate), (variable[edge], True)))
        return eqn

    def form_mixed_equation(self, profiles, names, coefficients, directions, edges, normal_profiles, normal_coeffs, poly_coordinates):
        """ Generates the equations for the 3D SBLI sidewall case.

        :arg list profiles: Arrays of values for the rhou and [0,1] temperature profiles.
        :arg list names: Variable names as strings.
        :arg list coefficients: Coefficients for the polynomial fit for rhou and temperature profiles.
        :arg list directions: Directions that contain a wall.
        :arg list edges: Indices for the boundary layer edges in each direction.
        :arg list normal_profiles: Arrays of values for the rhov and rhow profiles.
        :arg list normal_coeffs: Coefficients for the polynomial fit for rhov and rhow.
        :returns: list: piecewise_eqns: Piecewise initialisation equations to be added to the initialisation class."""

        # Assuming we have the same number of poly coefficients in each direction, change later if required
        powers = [i for i in range(np.size(coefficients[0][0]))][::-1]
        # Loop over rhou, and T profiles
        piecewise_eqns = []
        direction1 = directions[0]
        direction2 = directions[1]
        # x variables and coordinates at the edge of the boundary layer in each
        coord1, coord2 = self.coordinates[direction1], self.coordinates[direction2]
        bl_coord1, bl_coord2 = poly_coordinates[edges[0]], poly_coordinates[edges[1]]
        # Create rhou profiles
        eqn1 = sum([coeff*coord1**power for (coeff, power) in zip(coefficients[0][0], powers)])
        eqn2 = sum([coeff*self.local_coordinate**power for (coeff, power) in zip(coefficients[1][0], powers)])  # profiles[0][max(edges)]
        u_var1, u_var2 = GridVariable('%s_1' % names[0]), GridVariable('%s_2' % names[0])
        # freestream_value = np.max([profiles[0][0][edges[0]], profiles[1][0][edges[1]]])
        freestream_value = 1.0
        piecewise_eqns.append(OpenSBLIEq(u_var1, Piecewise((eqn1, coord1 < bl_coord1), (freestream_value, True))))
        piecewise_eqns.append(OpenSBLIEq(u_var2, Piecewise((eqn2, self.local_coordinate < bl_coord2), (freestream_value, True))))
        piecewise_eqns.append(OpenSBLIEq(GridVariable('%s' % names[0]), u_var1*u_var2))
        # Create T profiles, g = (T-Tw)/(Tinf - Tw) ---> T = g*(Tinf-Tw) + Tw
        Tw = ConstantObject('Twall')
        eqn1 = sum([coeff*coord1**power for (coeff, power) in zip(coefficients[0][1], powers)])
        eqn2 = sum([coeff*self.local_coordinate**power for (coeff, power) in zip(coefficients[1][1], powers)])  # profiles[0][max(edges)]
        T_var1, T_var2 = GridVariable('%s_1' % names[1]), GridVariable('%s_2' % names[1])
        freestream_value = np.max([profiles[0][1][edges[0]], profiles[1][1][edges[1]]])
        freestream_value = 1.0  # CHECK THIS VALUE
        piecewise_eqns.append(OpenSBLIEq(T_var1, Piecewise((eqn1, coord1 < bl_coord1), (freestream_value, True))))
        piecewise_eqns.append(OpenSBLIEq(T_var2, Piecewise((eqn2, self.local_coordinate < bl_coord2), (freestream_value, True))))
        piecewise_eqns.append(OpenSBLIEq(GridVariable('%s' % names[1]), T_var1*T_var2*(freestream_value - Tw) + Tw))
        # Create normal velocity component profiles, rhov:
        eqn1 = sum([coeff*coord1**power for (coeff, power) in zip(normal_coeffs[0], powers)])
        temp1 = GridVariable('%s' % names[2])
        rhov_inf = normal_profiles[0][edges[0]]*u_var2  # Multiplying by rhou from the other direction
        piecewise_eqns.append(OpenSBLIEq(temp1, Piecewise((eqn1*u_var2, coord1 < bl_coord1), (rhov_inf, True))))
        # rhow:
        eqn2 = sum([coeff*self.local_coordinate**power for (coeff, power) in zip(normal_coeffs[1], powers)])
        temp2 = GridVariable('%s' % names[3])
        # rhow should reduce to zero at the symmetry plane
        rhow_inf = normal_profiles[1][edges[1]]*(1 - (self.local_coordinate-bl_coord2)/(0.5*self.domain_length-bl_coord2))*u_var1
        piecewise_eqns.append(OpenSBLIEq(temp2, self.side_fact*Piecewise((eqn2*u_var1, self.local_coordinate < bl_coord2), (rhow_inf, True))))
        return piecewise_eqns

    def generate_one_wall_equations(self, data, coeffs, direction, edge, poly_coordinates):
        """ Generates the equations for 2D SBLI and 3D span-periodic cases.

        :arg list data: Profile arrays for rhou0, rhou1 and temperature.
        :arg list coeffs: Coefficients for the polynomial fits.
        :arg list direction: Direction normal to the wall.
        :arg int edge: Grid index for the edge of the boundary-layer."""

        self.eqns = []
        rhou0_eqn = self.form_equation(data[0], 'rhou0', coeffs[0], direction, edge, poly_coordinates)
        rhou1_eqn = self.form_equation(data[1], 'rhou1', coeffs[1], direction, edge, poly_coordinates)
        T_eqn = self.form_equation(data[2], 'T', coeffs[2], direction, edge, poly_coordinates)

        # densities of species
        cN2_eqn = self.form_equation(data[3], 'cN2', coeffs[3], direction, edge, poly_coordinates)
        cN_eqn = self.form_equation(data[4], 'cN', coeffs[4], direction, edge, poly_coordinates)
        cO2_eqn = self.form_equation(data[5], 'cO2', coeffs[5], direction, edge, poly_coordinates)
        cO_eqn = self.form_equation(data[6], 'cO', coeffs[6], direction, edge, poly_coordinates)
        cNO_eqn = self.form_equation(data[7], 'cNO', coeffs[7], direction, edge, poly_coordinates)

        rho_eqn = self.form_equation(data[8], 'rho', coeffs[8], direction, edge, poly_coordinates)

        # Set conservative values
        rho, rhou0, rhou1, T = GridVariable('rho'), GridVariable('rhou0'), GridVariable('rhou1'), GridVariable('T')
        cN2, cN, cO2, cO, cNO = GridVariable('cN2'), GridVariable('cN'), GridVariable('cO2'), GridVariable('cO'), GridVariable('cNO')
        rhoN2, rhoN, rhoO2, rhoO, rhoNO = GridVariable('rhoN2'), GridVariable('rhoN'), GridVariable('rhoO2'), GridVariable('rhoO'), GridVariable('rhoNO')
        
        rhoev = GridVariable('rhoev')
        evequilO2, evequilN2, evequilNO = GridVariable('evequilO2'), GridVariable('evequilN2'), GridVariable('evequilNO')

        # set dimensional grid constants -------------------------------------------------------------------------------------------------------------------------------
        uref, rhoref, pref, Rhat, Twall = ConstantObject('uref'), ConstantObject('rhoref'), ConstantObject('pref'), ConstantObject('Rhat'), ConstantObject('Twall')
        thetavO, thetavN, thetavO2, thetavN2, thetavNO = ConstantObject('thetavO'),ConstantObject('thetavN'),ConstantObject('thetavO2'),ConstantObject('thetavN2'),ConstantObject('thetavNO')
        MO, MN, MO2, MN2, MNO = ConstantObject('MO'), ConstantObject('MN'), ConstantObject('MO2'), ConstantObject('MN2'), ConstantObject('MNO')
        dhO, dhN, dhO2, dhN2, dhNO = ConstantObject('dhO'), ConstantObject('dhN'), ConstantObject('dhO2'), ConstantObject('dhN2'), ConstantObject('dhNO')

        cOe, cNe, cO2e, cN2e, cNOe = ConstantObject('cOe'), ConstantObject('cNe'), ConstantObject('cO2e'), ConstantObject('cN2e'), ConstantObject('cNOe')
        sigOtoNe = ConstantObject('sigOtoNe')
        CTD.add_constant(sigOtoNe)

        # rho_eqn = OpenSBLIEq(rho, 1.0/T)
        # rho_store = OpenSBLIEq(DataObject('rho'), rhoN2 + rhoN + rhoO2 + rhoO + rhoNO )
        rho_store = OpenSBLIEq(rho,  (pref/(Rhat*(cN2/MN2 + cN/MN + cO2/MO2 + cO/MO + cNO/MNO)))*(1.0/T))

        rhoN2_pre = OpenSBLIEq(GridVariable('rhoN2'), rho*cN2)
        rhoN_pre = OpenSBLIEq(GridVariable('rhoN'), rho*cN)
        rhoO2_pre = OpenSBLIEq(GridVariable('rhoO2'), rho*cO2)
        rhoO_pre = OpenSBLIEq(GridVariable('rhoO'), rho*cO)
        rhoNO_pre = OpenSBLIEq(GridVariable('rhoNO'), rho*cNO)

        # rho_store = OpenSBLIEq(DataObject('rho'), rho)
        rhou0_store = OpenSBLIEq(DataObject('rhou0'), rhou0)
        rhou1_store = OpenSBLIEq(DataObject('rhou1'), rhou1)

        rhoN2_store = OpenSBLIEq(DataObject('rhoN2'), rhoN2)
        rhoN_store = OpenSBLIEq(DataObject('rhoN'), rhoN)
        rhoO2_store = OpenSBLIEq(DataObject('rhoO2'), rhoO2)
        rhoO_store = OpenSBLIEq(DataObject('rhoO'), rhoO)
        rhoNO_store = OpenSBLIEq(DataObject('rhoNO'), rhoNO)

        # Tinf = GridVariable('T')
        # tempeq = OpenSBLIEq(Tinf, pref/(Rhat*(rhoO/MO+rhoO2/MO2+rhoN/MN+rhoN2/MN2+rhoNO/MNO)))
        Tempwall = GridVariable('Twalld')
        tempwall = OpenSBLIEq(Tempwall,  pref/(Rhat*(rhoO/MO+rhoO2/MO2+rhoN/MN+rhoN2/MN2+rhoNO/MNO)) + Twall + cOe+ cNe+ cO2e+ cN2e+ cNOe + sigOtoNe)
        # equilibrium quantities
        evequilO2eq = OpenSBLIEq(evequilO2, thetavO2*Rhat/(MO2*(exp(thetavO2/T)-1.0)))
        evequilN2eq = OpenSBLIEq(evequilN2, thetavN2*Rhat/(MN2*(exp(thetavN2/T)-1.0)))
        evequilNOeq = OpenSBLIEq(evequilNO, thetavNO*Rhat/(MNO*(exp(thetavNO/T)-1.0)))
        rhoev_store = OpenSBLIEq(DataObject('rhoev'), evequilN2*rhoN2 + evequilNO*rhoNO+ evequilO2*rhoO2)
        

        rhof = OpenSBLIEq(DataObject('rhof'), rhoN)
        # gama, Minf = ConstantObject("gama"), ConstantObject("Minf")
        rhoE_store = OpenSBLIEq(DataObject('rhoE'), Rhat*T*(3.0/2.0*(rhoO/MO+rhoN/MN)+5.0/2.0*(rhoO2/MO2+rhoN2/MN2+rhoNO/MNO))+4.1868e6*(dhO*rhoO/MO+dhN*rhoN/MN+dhNO*rhoNO/MNO)+ (rhoO2*evequilO2+rhoN2*evequilN2+rhoNO*evequilNO )+0.5*((rhou0/rho)**2+(rhou1/rho)**2)*rho)
        # rhoE_store = OpenSBLIEq(DataObject('rhoE'), Rhat*T*(3.0/2.0*(rhoOi/MO+rhoNi/MN)+5.0/2.0*(rhoO2i/MO2+rhoN2i/MN2+rhoNOi/MNO))+4.1868e6*(dhO*rhoOi/MO+dhN*rhoNi/MN+dhNO*rhoNOi/MNO)+rhoO2i*evequilO2+rhoN2i*evequilN2+rhoNOi*evequilNO+0.5*((rhou0/rho)**2+(rhou1/rho)**2)*rho)
        


        # gama, Minf = ConstantObject("gama"), ConstantObject("Minf")
        # rhoE_store = OpenSBLIEq(DataObject('rhoE'), rho*T/(gama*(gama-1)*Minf**2) + 0.5*(rhou0**2 + rhou1**2)/rho)


        self.eqns += [rhou0_eqn, rhou1_eqn, T_eqn, cN2_eqn, cN_eqn, cO2_eqn, cO_eqn, cNO_eqn, rho_store, rhou0_store, rhou1_store, rhoN2_pre,  rhoN_pre, rhoO2_pre, rhoO_pre, rhoNO_pre, rhoN2_store, rhoN_store, rhoO2_store, rhoO_store, rhoNO_store, evequilO2eq, evequilN2eq, evequilNOeq, rhoev_store, rhoE_store, rhof, tempwall]
        if self.block.ndim == 3:  # Periodic case, rhow = 0
            self.eqns += [OpenSBLIEq(DataObject('rhou2'), 0.0)]
        return

    def generate_two_wall_equations(self, profiles, coeffs, directions, edges, normal_profile, normal_coeffs, poly_coordinates):
        """ Generates the equations for 3D case with sidewalls.

        :arg list profiles: Profile arrays for rhou0 and temperature.
        :arg list coeffs: Coefficients for the polynomial fits.
        :arg list directions: Directions normal to the wall.
        :arg list edges: Grid indexes for the edge of the boundary-layer.
        :arg list normal_profile: Profile for the wall normal velocity components.
        :arg list normal_coeffs: Coefficients for the wall normal polynomial fit."""
        self.eqns = []
        names = ['rhou0', 'T', 'rhou1', 'rhou2']
        # Hard code to spanwise direction for 2 walls
        direction = 2
        array_coord = self.coordinates[direction]
        self.domain_length = self.block.deltas[direction]*(self.block.shape[direction]-1)
        self.side_fact = GridVariable('side_fact')
        side_fact_eqn = OpenSBLIEq(self.side_fact, Piecewise((1, array_coord <= self.domain_length/2.0), (-1, True)))

        self.local_coordinate = GridVariable('local_coord')
        local_coord_eqn = OpenSBLIEq(self.local_coordinate, Piecewise((array_coord, array_coord <= self.domain_length/2.0), (self.domain_length-array_coord, True)))
        self.eqns += [local_coord_eqn]
        self.eqns += [side_fact_eqn]
        # Create the piecewise equations formed from the boundary layers
        bl_equations = self.form_mixed_equation(profiles, names, coeffs, directions, edges, normal_profile, normal_coeffs, poly_coordinates)
        # Set conservative values
        rho, rhou0, rhou1, rhou2, T = GridVariable('rho'), GridVariable('rhou0'), GridVariable('rhou1'), GridVariable('rhou2'), GridVariable('T')
        rho_eqn = OpenSBLIEq(rho, 1.0/T)
        rho_store = OpenSBLIEq(DataObject('rho'), rho)
        rhou0_store = OpenSBLIEq(DataObject('rhou0'), rhou0)
        rhou1_store = OpenSBLIEq(DataObject('rhou1'), rhou1)
        rhou2_store = OpenSBLIEq(DataObject('rhou2'), rhou2)
        gama, Minf = ConstantObject("gama"), ConstantObject("Minf")
        rhoE_store = OpenSBLIEq(DataObject('rhoE'), rho*T/(gama*(gama-1)*Minf**2) + 0.5*(rhou0**2 + rhou1**2 + rhou2**2)/rho)
        self.eqns += bl_equations + [rho_eqn, rho_store, rhou0_store, rhou1_store, rhou2_store, rhoE_store]
        return

    def fit_polynomial(self, coords, variable, bl_edge, n_coeffs):
        """ Fits a polynomial to the input data, coefficients are returned.

        :arg ndarray coords: Independent variable of the input data.
        :arg ndarray variable: Dependent variable of the input data.
        :arg int bl_edge: Array index at the edge of the boundary-layer.
        :arg int n_coeffs: Desired number of coefficients for the polynomial.
        :returns: ndarray: coeffs: Coefficients of the polynomial fit."""
        coords = coords[0:bl_edge]
        variable = variable[0:bl_edge]
        # with warnings.catch_warnings():
        #     warnings.filterwarnings('error')
        #     try:
        warnings.filterwarnings("ignore")
        coeffs = poly.polyfit(coords, variable, n_coeffs)
        # except np.RankWarning:
        # print "Poorly conditioned fit"
        # ffit = poly.polyval(coords, coeffs)
        # plt.plot(coords, ffit, label='fit')
        # plt.plot(coords, variable, label='original_data')
        # plt.legend(loc="best")
        # plt.show()
        # Reverse coefficients so they are in descending order
        return coeffs[::-1]
    
    def find_exp(self,number) -> int:
        base10 = np.log10(abs(number))
        return np.floor(base10)

    def solved_continuity(self, y, u, rho):
        """ Solves the continuity equation to obtain the wall normal velocity profile.

        :arg ndarray y: Dependent coordinate values.
        :arg ndarray u: Streamwise velocity component values.
        :arg ndarray rho: Density values.
        :returns: ndarray: rhov: Array of values for the wall normal velocity components. """
        # Grid offset delta to form derivative approximation
        n = np.size(y)
        ya2 = y[:]
        # delta, scale, re = 0.001, self.scale_factor, self.Re
        delta, scale, re = 0.01*10**self.find_exp(self.blthickness), self.scale_factor, self.Re

        print('-------------------------------------------- printing out delta value :',delta)
        rex0 = 0.5*(re/scale)**2
        x0 = 0.5*re/scale**2
       
        drudx, rhov = np.zeros_like(y), np.zeros_like(y)
        # Local Reynolds number scaling to obtain a v profile
        sqrex = np.sqrt(rex0)
        delsx = np.sqrt(2.0)*scale*(x0)/sqrex
        ya2 = delsx*ya2
        
        d2y_u = spline(ya2, u, n, self.dudy, 0)
        d2y_rho = spline(ya2, rho, n, 0, 0)
        for j in range(0, n):
            ya2[j] = delsx*ya2[j]
            dstarp = delsx*np.sqrt((x0+delta)/(x0))
            dstarm = delsx*np.sqrt((x0-delta)/(x0))
            yp, ym = ya2[j]/dstarp, ya2[j]/dstarm
            uxp = splint(ya2, u, d2y_u, n, yp)
            uxm = splint(ya2, u, d2y_u, n, ym)
            rhoxp = splint(ya2, rho, d2y_rho, n, yp)
            rhoxm = splint(ya2, rho, d2y_rho, n, ym)
            drudx[j] = (rhoxp*uxp-rhoxm*uxm)/(2.0*delta)
            rhov[j] = rhov[j-1]-0.5*(ya2[j]-ya2[j-1])*(drudx[j]+drudx[j-1])

        
        dy = y[1]
        dvdy = (-1.83333333333334*rhov[0]+3.00000000000002*rhov[1]-1.50000000000003*rhov[2]+0.333333333333356*rhov[3]-8.34657956545823e-15*rhov[4]+1.06910315192207e-15*rhov[5])/dy
        
        v_n = rhov 
        dvdy = (-1.83333333333334*v_n[0]+3.00000000000002*v_n[1]-1.50000000000003*v_n[2]+0.333333333333356*v_n[3]-8.34657956545823e-15*v_n[4]+1.06910315192207e-15*v_n[5])/dy

        return v_n, dvdy
    
    def solve_continuity(self, y, u, rho):
        """ Solves the continuity equation to obtain the wall normal velocity profile.

        :arg ndarray y: Dependent coordinate values.
        :arg ndarray u: Streamwise velocity component values.
        :arg ndarray rho: Density values.
        :returns: ndarray: rhov: Array of values for the wall normal velocity components. """
        # Grid offset delta to form derivative approximation
        n = np.size(y)
        ya2 = y[:]
        delta, scale, re = 0.001, self.scale_factor, self.Re
        rex0 = 0.5*(re/scale)**2
        x0 = 0.5*re/scale**2
        # print "Domain inlet is when the boundary-layer has developed for a length of x0 = %.10f" % x0
        drudx, rhov = np.zeros_like(y), np.zeros_like(y)
        # Local Reynolds number scaling to obtain a v profile
        sqrex = np.sqrt(rex0)
        delsx = np.sqrt(2.0)*scale*(x0)/sqrex
        ya2 = delsx*ya2
        d2y_u = spline(ya2, u, n, self.dudy, 0)
        d2y_rho = spline(ya2, rho, n, 0, 0)
        for j in range(0, n):
            ya2[j] = delsx*ya2[j]
            dstarp = delsx*np.sqrt((x0+delta)/(x0))
            dstarm = delsx*np.sqrt((x0-delta)/(x0))
            yp, ym = ya2[j]/dstarp, ya2[j]/dstarm
            uxp = splint(ya2, u, d2y_u, n, yp)
            uxm = splint(ya2, u, d2y_u, n, ym)
            rhoxp = splint(ya2, rho, d2y_rho, n, yp)
            rhoxm = splint(ya2, rho, d2y_rho, n, ym)
            drudx[j] = (rhoxp*uxp-rhoxm*uxm)/(2.0*delta)
            rhov[j] = rhov[j-1]-0.5*(ya2[j]-ya2[j-1])*(drudx[j]+drudx[j-1])
        return rhov

    def load_similarity(self):
        """ Solves the compressible boundary-layer equations via similarity solution."""
        Re, xMach, Tinf, Tw, Sc = self.Re, self.xMach, self.Tinf, self.Tw, self.Sc
        cN2e, cNe, cO2e, cOe, cNOe = self.cN2e, self.cNe, self.cO2e, self.cOe, self.cNOe
        # pref, rhoref, uref, blthickness = self.pref, self.rhoref, self.uref, self.blthickness
        pref = self.pref 

        adiabatic_condition, catalytic_condition = self.adiabatic_condition, self.catalytic_condition



        Pr, gama = 0.72, 1.4  # Prandtl number, ratio of specific heats
        # bl = Boundary_layer_profile(xMach, Pr, gama, -1, Re, Tinf)  # -1 for Tw sets an adiabatic wall
        bl = Boundary_layer_profile(Re, xMach, Tinf, Tw, Sc, adiabatic_condition, catalytic_condition, cN2e, cNe, cO2e, cOe, cNOe, pref)
        y, u, T, rho, n = bl.y, bl.u, bl.T, bl.rho, np.size(bl.y)
        uref, Tinf, rhoref, blthickness, rho = bl.uref, bl.tinf, bl.rhoref, bl.blthickness, bl.rho
        cN2, cN, cO2, cO, cNO = bl.cN2, bl.cN, bl.cO2, bl.cO, bl.cNO

        # self.uref, self.tinf, self.rhoref, self.blthickness 

        self.Twall, self.scale_factor = bl.Twall, bl.scale_factor  # Wall temperature and scale factor from the similarity solution
        self.dudy = bl.dudy  # du/dy at the wall
        return y, u, T, rho, n, uref, Tinf, rhoref, blthickness, cN2, cN, cO2, cO, cNO
    
    def load_similarity_reference_variables(self):
        """ Solves the compressible boundary-layer equations via similarity solution."""
        Re, xMach, Tinf = self.Re, self.xMach, self.Tinf
        Pr, gama = 0.72, 1.4  # Prandtl number, ratio of specific heats
        bl = Boundary_layer_profile(xMach, Pr, gama, -1, Re, Tinf)  # -1 for Tw sets an adiabatic wall
        # y, u, T, rho, n = bl.y, bl.u, bl.T, 1.0/bl.T, np.size(bl.y)
        uref, tinf, rhoref, blthickness, rho = bl.uref, bl.tinf, bl.rhoref, bl.blthickness, bl.rho
        cN2, cN, cO2, cO, cNO = bl.cN2, bl.cN, bl.cO2, bl.cO, bl.cNO

        # self.uref, self.tinf, self.rhoref, self.blthickness 

        # self.Twall, self.scale_factor = bl.Twall, bl.scale_factor  # Wall temperature and scale factor from the similarity solution
        # self.dudy = bl.dudy  # du/dy at the wall
        return uref, tinf, rhoref, blthickness

    def interpolate_onto_grid(self, y_in, y_out, var_in, y0, yn):
        # Create interpolating second derivative spline
        # n = size of the original data
        n = np.size(y_in)
        n_out = np.size(y_out)
        d2y = spline(y_in, var_in, n, y0, yn)
        # Array for variable interpolated onto the grid
        var_out = np.zeros_like(y_out)
        for i in range(n_out):
            var_out[i] = splint(y_in, var_in, d2y, n, y_out[i])
        return var_out

    def find_edge_of_bl(self, variable, tolerance):
        """ Finds the edge of the boundary layer and returns the index of that grid point.

        :arg ndarray variable: Array of values for a given flow variable.
        :arg float tolerance: Stopping tolerance for the difference between two successive grid points.
        :returns: int: index: Index of the boundary-layer edge."""
        index = 1
        while np.abs(variable[index]-variable[index-1]) > tolerance:
            index += 1
        return index

    def freestream_value(self, variable, index):
        """ Returns the value of a flow variable for a given grid index.

        :arg ndarray variable: Array of values for a given flow variable.
        :arg index Array index to use.
        :returns: float: variable[index]: Value of the flow variable at that index."""
        return variable[index]

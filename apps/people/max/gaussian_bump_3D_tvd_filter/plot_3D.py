import numpy
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import h5py
import os.path
import matplotlib.cm as cm
import os


class plotFunctions(object):
    def __init__(self):
        return

    def contour_local(self, fig, levels0, label, variable):
        ax1 = fig.add_subplot(1, 1, 1, aspect='equal')
        ax1.set_xlabel(r"$x_0$", fontsize=20)
        ax1.set_ylabel(r"$x_1$", fontsize=20)
        # ax1.set_xlim([0, 0.1])
        # ax1.set_ylim([0, 0.1])
        CS = ax1.contourf(self.x[0,:,:], self.y[0,:,:], variable[0,:,:], levels=levels0, cmap=cm.jet)
        divider = make_axes_locatable(ax1)
        cax1 = divider.append_axes("right", size="5%", pad=0.05)
        ticks_at = numpy.linspace(levels0[0], levels0[-1], 10)
        cbar = plt.colorbar(CS, cax=cax1, ticks=ticks_at, format='%.3f')
        cbar.ax.set_ylabel(r"$%s$" % label, fontsize=20)
        return

    def read_file(self, fname):
        f = h5py.File(fname, 'r')
        group = f["opensbliblock00"]
        return f, group

    def read_dataset(self, group, dataset):
        d_m = group["%s" % (dataset)].attrs['d_m']
        size = group["%s" % (dataset)].shape
        read_start = [abs(d) for d in d_m]
        read_end = [s-abs(d) for d, s in zip(d_m, size)]
        if len(read_end) == 2:
            read_data = group["%s" % (dataset)][read_start[0]:read_end[0], read_start[1]:read_end[1]]
        elif len(read_end) == 3:
            read_data = group["%s" % (dataset)][read_start[0]:read_end[0], read_start[1]:read_end[1], read_start[2]:read_end[2]]
        else:
            raise NotImplementedError("")
        return read_data


class plots_3D(plotFunctions):
    def __init__(self):
        self.Minf = 4.0
        self.Re = 4000
        self.RefT = 439.0
        self.SuthT = 110.4
        self.Ly = 115.0
        self.Lx = 400.0
        self.scale = 4.959043
        self.D11 = self.extract_metrics()
        self.Pr = 0.72
        return

    def extract_coordinates(self):
        f, group1 = self.read_file(fname)
        x = self.read_dataset(group1, "x0_B0")
        y = self.read_dataset(group1, "x1_B0")
        z = self.read_dataset(group1, "x2_B0")
        dx, dy, dz = x[0, 0, 1], y[0, 1, 0], z[1,0,0]
        print("Grid size (x,y)  is: (%f, %f)" % (x.shape[1], x.shape[0]))
        print("First grid point dx: %f, dy: %f" % (dx, dy))
        dx, dy, dz = x[0,0, -1]-x[0,0, -2], y[0,-1, 0]-y[0,-2, 0], z[-1,0, 0]-z[-2,0, 0]
        print("Last grid point dx: %f, dy: %f" % (dx, dy))
        return x, y, z

    def viscous_length_scale(self, fname):
        f, group1 = self.read_file(fname)
        rho, u, v, rhoE, p, T, M, mu = self.extract_flow_variables(group1)
        Cf, tau_wall = self.compute_skin_friction(u, mu)
        # Wall viscosity all x points
        mu_wall = mu[0, 0, :]
        # print(rho[0,0,:])
        u_tau = numpy.sqrt(abs(tau_wall / rho[0,0,:])) # friction velocity
        # print(u_tau)
        dv = mu_wall / (rho[0,0,:]*u_tau)  # viscous length scale
        return dv, u_tau

    def extract_metrics(self):
        fname = 'opensbli_output.h5'
        f, group1 = self.read_file(fname)
        D11 = self.read_dataset(group1, "D11_B0")
        return D11

    def extract_flow_variables(self, group):
        rho = self.read_dataset(group, "rho_B0")
        rhou = self.read_dataset(group, "rhou0_B0")
        rhov = self.read_dataset(group, "rhou1_B0")
        rhow = self.read_dataset(group, "rhou2_B0")
        rhoE = self.read_dataset(group, "rhoE_B0")
        u = rhou/rho
        v = rhov/rho
        w = rhow/rho
        p = (0.4)*(rhoE - 0.5*(u**2+v**2+w**2)*rho)
        a = numpy.sqrt(1.4*p/rho)
        M = numpy.sqrt(u**2 + v**2 + w**2)/a
        T = 1.4*(self.Minf**2)*p/rho
        mu = self.compute_viscosity(T)
        return rho, u, v, rhoE, p, T, M, mu

    def compute_wall_derivative(self, variable):
        ny = numpy.size(self.y[:, 0])
        Ly = self.Ly
        delta = Ly/(ny-1.0)
        D11 = self.D11[0, 0:6, :]
        var = variable[0, 0:6, :]
        coeffs = numpy.array([-1.83333333333334, 3.00000000000002, -1.50000000000003, 0.333333333333356, -8.34617916606957e-15, 1.06910884386911e-15])
        coeffs = coeffs.reshape([6, 1])
        dudy = sum(D11*var*coeffs)/delta
        return dudy

    def compute_viscosity(self, T):
        mu = (T**(1.5)*(1.0+self.SuthT/self.RefT)/(T+self.SuthT/self.RefT))
        return mu

    def compute_skin_friction(self, u, mu):
        # Wall viscosity all x points
        mu_wall = mu[0, 0, :]
        dudy = self.compute_wall_derivative(u)
        tau_wall = dudy*mu_wall
        Cf = tau_wall/(0.5*self.Re)
        return Cf, tau_wall
    
    def compute_heat_flux(self,T,mu):
        dTdy = self.compute_wall_derivative(T)
        mu_wall = mu[0, 0, :]
        q = (mu_wall/(0.4*self.Minf**2*self.Pr*self.Re)*dTdy)
        return q
    
    def St_Rex(self,Cf):
        Rex = 0.5*(self.Re / self.scale)**2 + self.x[0, 1, :]*self.Re
        St = Cf / (2*self.Pr**(2/3))
        return Rex,St

    def SBLI_comparison(self, Cf, P):

        fig = plt.figure()
        ax = fig.add_subplot(111)

        dv, u_tau = self.viscous_length_scale(fname)
        dx, dy, dz = self.x[0, 0, 1], self.y[0, 1, 0], self.z[1,0,0]
        dx_plus, dy_plus, dz_plus = dx/ dv, dy/ dv, dz/ dv

        loc=346
        ax.plot(self.x[0, 1, :], Cf, color='k', label='OpenSBLI')
        ax.axhline(y=0.0, linestyle='--', color='k')
        ax.annotate("At x=%d \n$dx^+$: %.2f \n$dy^+$: %.2f \n$dz^+$: %.2f" % (loc, dx_plus[loc], dy_plus[loc], dz_plus[loc]), xy =(348,0.31), xycoords='data', xytext=(-90,-5),
            textcoords='offset points',bbox=dict(boxstyle="round", fc="0.8"), arrowprops=dict(arrowstyle="->"))
        ax.set_xlabel(r'$x_0$', fontsize=20)
        ax.set_ylabel(r'$C_f$', fontsize=20)
        # ax.set_ylim([-0.004, 0.004])
        ax.set_xlim([min(self.x[0, 1, :]), max(self.x[0, 1, :])])
        ax.set_title('Skin friction')
        plt.legend(loc="best")
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        fig.savefig(directory + 'skin_friction.pdf', bbox_inches='tight')
        fig.clf()

        plt.plot(self.x[0, 1, :], P[0, 0, :]/P[0, 0, 0], color='k', label="OpenSBLI")
        # linestyle='', marker='o',markevery=10)
        plt.xlabel(r'$x_0$', fontsize=20)
        plt.xlim(-10)
        plt.grid()
        plt.ylabel(r'$\frac{P_w}{P_1}$', fontsize=22)
        plt.title('Normalized wall pressure')
        plt.legend(loc="best")
        plt.savefig(directory + "wall_pressure.pdf", bbox_inches='tight')

        Rex,St = self.St_Rex(Cf)
        fig1, ax1 = plt.subplots()
        ax1.plot(Rex, St)
        ax1.set_xlabel(r'$Re_x$', fontsize=20)
        ax1.set_ylabel(r'$St$', fontsize=20)
        ax1.set_title('Reynolds number vs Stanton number', fontsize=15)
        fig1.savefig(directory+'Rex_St.pdf', bbox_inches='tight')
        plt.clf()
        return

    def main_plot(self, fname, n_levels):
        f, group1 = self.read_file(fname)

        self.x, self.y, self.z = self.extract_coordinates()
        rho, u, v, rhoE, p, T, M, mu = self.extract_flow_variables(group1)
        
        variables = [rho, u, v, rhoE, p, M, T, mu]
        names = ["rho", "u", "v", "rhoE", "P", "M", "T", "mu"]

        # Contour plots
        for var, name in zip(variables, names):
            min_val = numpy.min(var)
            max_val = numpy.max(var)
            levels = numpy.linspace(min_val, max_val, n_levels)
            print("%s" % name)
            print(levels)
            fig = plt.figure()
            self.contour_local(fig, levels, "%s" % name, var)
            plt.savefig(directory + "gaussian_bump_3D_%s.pdf" % name, bbox_inches='tight')
            plt.clf()
        # Compare to SBLI
        Cf, tau_wall = self.compute_skin_friction(u, mu)
        self.SBLI_comparison(Cf, p)

        q = self.compute_heat_flux(T,mu)
        fig, ax = plt.subplots()
        ax.plot(self.x[0, 1, :], q, color='k') #, label='OpenSBLI'
        ax.set_xlabel(r'$x_0$', fontsize=15)
        ax.set_ylabel(r'$\dot{q}$', fontsize=15)
        ax.set_xlim([min(self.x[0, 1, :]), max(self.x[0, 1, :])])
        ax.set_title('Wall heat flux')
        # plt.legend(loc="best")
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        fig.savefig(directory + 'heat_flux.pdf', bbox_inches='tight')


        f.close()

    def turbulence_plots(self, fname):
        f, group1 = self.read_file(fname)
        # self.x, self.y, self.z = self.extract_coordinates()
        dv, u_tau = self.viscous_length_scale(fname)
        rho, u, v, rhoE, p, T, M, mu = self.extract_flow_variables(group1)
        dx, dy, dz = self.x[0, 0, 1], self.y[0, 1, 0], self.z[1,0,0]
        dx_plus, dy_plus, dz_plus = dx/ dv, dy/ dv, dz/ dv
        loc = 200
        print(r'dx+: %f, dy+: %f, dz+: %f' % (dx_plus[loc], dy_plus[loc], dz_plus[loc]))
        
        u_plus = u[0,1] / u_tau
        y_plus = dy*u_tau*rho[0,0,:] / mu[0,0,:]
        fig, ax = plt.subplots()
        ax.plot(self.x[0,0,:],y_plus)
        ax.set_ylabel(r'$Y^+$', fontsize=20)
        ax.set_xlabel(r'$x_0$', fontsize=20)
        ax.set_xlim([0,400])
        ax.grid()
        # ax.set_ylim([-0.004, 0.004])
        ax.set_title(r'$Y^+$ value along the wall',fontsize=15)
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        fig.savefig(directory + "y_plus.pdf", bbox_inches='tight')

fname = "opensbli_output.h5"
n_contour_levels = 25
directory = './simulation_plots/'

if not os.path.exists(directory):
    os.makedirs(directory)

plots = plots_3D()
plots.main_plot(fname, n_contour_levels)
plots.turbulence_plots(fname)
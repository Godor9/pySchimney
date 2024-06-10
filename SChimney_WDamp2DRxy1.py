# -*- coding: utf-8 -*-
"""
Created on Wed May  1 16:14:58 2024
@author: Lawrence H.W
This is the python version of the 2D full coupled porosity wave code

"""
from scipy.io import savemat
import numpy as np
import matplotlib.pyplot as plt
import sys
#from time import time

# Constants and Parameters
nt = 3000
phi0 = 0.1
rhos = 2
rhof = 1
mut = 0.1  # total rock (solid+fluid) shear viscosity
muf = 1     # shear viscosity of the fluid
kperm0 = 1  # permeability coefficient
#eta2mut = 10
etas0   = 1  #eta2mut * mut / phi0  # background bulk viscosity
eta2mut = etas0*phi0/mut
npower = 3.0  # power exponent of the permeability-porosity relationship
Lx = 20  #10*np.sqrt(etas0 * kperm0 / muf)
Rxy = 2
Ly = Lx * Rxy
gravity = 1  # simplified gravity value
clength = 1
ctime_unit = clength * phi0 / kperm0 / (rhos - rhof) / gravity * muf
tsc        = ctime_unit
ttotal = 8 * ctime_unit
R = 100
Vplastic_on = 0
pkc = 3.0
p0 = 4.05
tau0 = 18
ntau = 2
pk0 = pkc
alpha_tau = 0
mp = 2
Lambda = 10
relax_eta = 0.1
relax_kp = 0.1

# Numerics
nx = 25 * 4
ny = 25 * 8
dx = Lx / nx
dy = Ly / ny
x = np.linspace(0.5 * dx, Lx - 0.5 * dx, nx)
y = np.linspace(0.5 * dy, Ly - 0.5 * dy, ny)
xx, yy = np.meshgrid(x, y,indexing='ij')
nout= 10
ploton=1

CN         = 0.5;# Crank-Nicolson CN=0.5, Backward Euler CN=0.0
force_iter = 0.5*ny;##1*ny;
iterMax1   = 100*ny;
Vdmp       = 3.14*2;
eta_b      = 0.5;

cnt  = 50;    #check/calculate PT residule every cnt iteration
cnt1  = 500;    #output residule every cnt iteration
epsi = 1e-5;
Xd   = 5.1
Ptsc = 3.0

Vsc  = 1.5
lambda_p= 0.01 #2*rhos*gravity*dy;

# Initial Conditions
phi = phi0*np.ones((nx, ny))
#phi = phi0 * (1 + 3 * np.exp(-((yy - 5) / 0.5) ** 2))
#initialize phi with a inclusion that has high porosity.
xa  = 1
phiA= 3
radc=((xx - 0.5*Lx) /4/xa) ** 2 + ((yy - 0.2*Ly)/xa)**2;
phi[radc<1]=phi0+phi0*phiA
#phi = phi0 * (1 + 3 * np.exp(-(((xx - 5) / 1.5) ** 2 + ((yy - 10) / 0.8) ** 2))) 
# initialize phi with an elliptical area that has high value
 
#phi[yy < 5] = phi0 * (1 + 3)
Pt = np.zeros((nx, ny))
Pe = np.zeros((nx, ny))
Pd = np.zeros((nx, ny))
Vx = np.zeros((nx + 1, ny))
Vy = np.zeros((nx, ny + 1))
div = np.zeros((nx, ny))
div_old = np.zeros((nx, ny))
Tauxx = np.zeros((nx, ny))
Tauyy = np.zeros((nx, ny))
Tauxy = np.zeros((nx + 1, ny + 1))
RVx = np.zeros((nx-1, ny))
RVy = np.zeros((nx, ny-1))
qDx = np.zeros((nx+1, ny))
qDy = np.zeros((nx, ny+1))
divqD=np.zeros((nx,ny))
dphidt = np.zeros((nx, ny))
tau = np.zeros((nx, ny))
phi_old = phi.copy();
F   = np.zeros((nx, ny))
ap  = np.ones((nx, ny))
kperm = kperm0 * np.ones((nx, ny))
etas  = etas0 * np.ones((nx, ny))
eta_phi = etas * phi0 / phi
kp_muf  =kperm/muf;

rhoBG = rhos * (1 - phi0) + rhof * phi0

#etas(cond2)=9*etas0; % not for shear viscosity, but only for bulk viscosity
#kperm(cond3)=0.1*kperm0;
#eta_phi  = etas0*ones(nx,ny);


phi0f   = phi.copy();
dt_fact = 1e-1*ctime_unit;
dt_init = 1e-5*ctime_unit;  #initial guess of dt
dt	    = dt_init;

phi0bc   = phi0 #mean(phi(:,end));
phi0bot  = phi[:,0];
# Constants and parameters initialization (define these as needed)
# Example: nt, ny, nx, CN, mut, eta_b, Ptsc, etc.

ndim = 2
betaqD = muf * (1 - phi0) / kperm0

# Calculate dtV
dtV = min(dx, dy)**2 / 2.1 / (1.0 + eta_b) / mut / ndim / Vsc
# Calculate dmp1 and dmp2
dmp1 = -1 / (2 * betaqD) + np.sqrt(1 + 4 * betaqD) / (2 * betaqD)
dmp2 = 1 / dmp1
qDy[:, [0, -1]] = dmp1 * (rhos - rhof) * gravity * (1.0 - phi0bc) * kperm0 / muf * (phi0bc / phi0)**npower
resid = 1e5;
residRVy = 1e5;
residRPt = 1e5;
residRPe = 1e5;
time  = 0;
TimeS=np.zeros(nt);
# Time integration
#tstart = time()
#plt.show()
for it in range(1,nt+1):
    # Update equations and other operations
    # Placeholder for actual simulation logic

# Time-stepping loop (assuming this is inside a loop over time steps)
    if it < 2:
        iterMax = 200 * ny
    else:
        iterMax = iterMax1
    
    np.copyto(phi_old, phi)
    np.copyto(div_old, div)
    #div_old = div.copy() % it will create new array every time!
    #div_old = div.copy()

    for iter in range(1, iterMax + 1):
        Tauxy1 = 0.5 * (Tauxy[:-1, :] + Tauxy[1:, :])
        Tauxy2 = 0.5 * (Tauxy1[:, :-1] + Tauxy1[:, 1:])
        tau = np.sqrt(0.5 * (Tauxx**2 + Tauyy**2 + 2.0 * Tauxy2**2))
        div = (np.diff(Vx, axis=0) / dx + np.diff(Vy, axis=1) / dy)
        phi = phi_old + dt * (1 - phi) * (CN * div_old + (1 - CN) * div)
    
        if Vplastic_on == 1:
            F = np.exp((np.abs(-alpha_tau * tau + Pe - p0) / pk0) ** mp - 1.0) * (1.0 + (tau / tau0) ** ntau) * pk0 - pk0
            f = pk0 * np.log(F / pk0 + 1.0)
            ind1 = F >= 0
            f[~ind1] = 0.0
            A = 2.0 * Lambda * f * np.exp(f / pk0) / pk0
            ap = np.ones_like(phi)
            Pd = np.zeros_like(phi)
            ap[ind1] = A[ind1] + 1.0
            Pd[ind1] = A[ind1] * (alpha_tau * tau[ind1] + p0)
            Pd = Pd / ap
            eta_phi = phi0 * etas / phi / ap * relax_eta + eta_phi * (1.0 - relax_eta)
        else:
            eta_phi = etas * phi0 / phi * (1.0 + 0.5 * (1.0 / R - 1.0) * (1.0 - np.tanh((Pe) / lambda_p))) * relax_eta + eta_phi * (1.0 - relax_eta)
    
        rhot = rhof * phi + rhos * (1 - phi) - rhoBG
        rhoPe = 1 / (1 - phi) / eta_phi
        kp_muf = np.exp(np.log(kperm / muf * (phi / phi0) ** npower) * relax_kp + np.log(kp_muf) * (1 - relax_kp))
        Vpref = np.sqrt(kp_muf * eta_phi)
        Xd1=np.minimum(Xd * (1 + 1.5*(iter / ny)), 2 * Xd)
        dtau = dx / Vpref / Xd1
        dtauPe = dtau / rhoPe
        dtauPe[1:-1, 1:-1] = np.min([dtauPe[1:-1, 1:-1], dtauPe[:-2, 1:-1], dtauPe[2:, 1:-1], dtauPe[1:-1, :-2], dtauPe[1:-1, 2:]], axis=0)
        dtauqDx = np.min([dtau[:-1, :], dtau[1:, :]], axis=0)
        dtauqDy = np.min([dtau[:, :-1], dtau[:, 1:]], axis=0)
    
        RPt = -div - (Pe - Pd) / eta_phi / (1 - phi)
        Ptsc1 = np.minimum(Ptsc * (1 + (iter / ny)), 5 * Ptsc)
        dtPt = 4.1 * mut * (1 + eta_b) / max(nx, ny) / Ptsc1
        Pt += dtPt * RPt
    
        Tauxx = 2.0 * mut * (np.diff(Vx, axis=0) / dx - div / 3.0 - eta_b * RPt)
        Tauyy = 2.0 * mut * (np.diff(Vy, axis=1) / dy - div / 3.0 - eta_b * RPt)
        Tauxy[1:-1, 1:-1] = mut * (np.diff(Vx[1:-1, :], axis=1) / dy + np.diff(Vy[:, 1:-1], axis=0) / dx)
    
        RVx = RVx * (1 - Vdmp / nx) + (np.diff(Tauxx - Pt, axis=0) / dx + np.diff(Tauxy[1:-1, :], axis=1) / dy)
        rho_avy = 0.5 * (rhot[:, :-1] + rhot[:, 1:])
        RVy = RVy * (1 - Vdmp / ny) + (np.diff(Tauyy - Pt, axis=1) / dy + np.diff(Tauxy[:, 1:-1], axis=0) / dx - rho_avy * gravity)
    
        Vx[1:-1, :] += dtV * RVx
        Vy[:, 1:-1] += dtV * RVy
    
        kpx_muf = 0.5 * (kp_muf[:-1, :] + kp_muf[1:, :])
        kpy_muf = 0.5 * (kp_muf[:, :-1] + kp_muf[:, 1:])
        RqDx = (kpx_muf * (np.diff(Pe - Pt, axis=0) / dx) - dmp2 * qDx[1:-1, :])
        RqDy = (kpy_muf * (np.diff(Pe - Pt, axis=1) / dy - (rhof - rhoBG) * gravity) - dmp2 * qDy[:, 1:-1])
        qDx[1:-1, :] += dtauqDx * RqDx / betaqD
        qDy[:, 1:-1] += dtauqDy * RqDy / betaqD
    
        divqD = np.diff(qDx, axis=0) / dx + np.diff(qDy, axis=1) / dy
        RPe = (divqD - dmp1 * (Pe - Pd) / (1 - phi) / eta_phi)
        Pe += dtauPe * RPe
        
        if iter % cnt == 0:
            residRVy = np.linalg.norm(RVy) / len(RVy.ravel())
            residRPe = np.linalg.norm(RPe) / len(RPe.ravel())
            residRPt = np.linalg.norm(RPt) / len(RPt.ravel())
            resid = max(residRVy, residRPe, residRPt)
        if iter % cnt1 == 0:
            #print(f'At Step {it},iter={iter} residRVy={residRVy:5.3e},residRPe={residRPe:5.3e},residRPt={residRPt:5.3e}')
            print('at step {:4d}, iter={}, residRVy={:5.3e}, residRPe={:5.3e}, residRPt={:5.3e}'.format(it, iter, residRVy, residRPe, residRPt))
    
        if residRVy==np.nan:
            print('At step {it},iter={:5d}. resid=nan. Exit now.'.format(it,iter))
            sys.exit()
        if resid > 1e5:
            #print(f'At Step {it}, Error is larger than 1e5 at iteration {iter}, resid={resid}')
            print(f'At Step {it}, Error is larger than 1e5 at iteration {iter}, resid={resid}')
            sys.exit()#return   
        if resid < epsi and iter >= force_iter:
            break
    time = time+ dt        
    TimeS[it] = time
    max_div = np.max(np.abs(div))
    dt = min(0.01 * tsc, dt_fact / (1e-10 + max_div))
    #print(f'Step={it:4d} converge , resid={resid:e}, dt={dt:e}, max_div={max_div:e} at iterations={iter:4d}(={iter/ny:3.1f} *ny)')
    print('Step={:4d} converge,time={:5.2f}({:5.2f} δt), resid={:e}, dt={:e}, max_div={:e} at iterations={:4d}(={:4.1f} *ny)'.format(it, time,time/tsc,resid, dt, max_div, iter, iter/ny))
    print('****maxporosity={:e},minporosity={:e},maxbulkv={:e},minbulk={:e},,maxperm={:e},minperm={:e}'.format(phi.max(),phi.min(),eta_phi.max(),eta_phi.min(),kp_muf.max(),kp_muf.min()))
    phi[phi<0.2*phi0]=0.2*phi0


# Assuming 'it', 'nout', 'x', 'y', 'Pe', 'phi0f', 'nx', 'phi', 'Ly', 'eta_phi', 'iter', 'iterMax1', 'resid', and 'epsi' are defined

    if it % nout == 0 and 1==0:  # not plotting
        plt.figure(2)
        plt.suptitle(f'step = {it}')
        # Subplot 1: Peclet number
        plt.subplot(2, 2, 1)
        plt.gca().set_aspect(2)  # pbaspect([1 2 1])
        plt.pcolor(xx, yy, Pe)  # Transpose Pe to match MATLAB's behavior
        plt.axis('image')
        plt.colorbar()
        plt.title(f'Pe at step = {it}')

        # Subplot 2: Porosity profiles
        plt.subplot(2, 2, 2)
        plt.plot(phi0f[int(nx/2), :], y, 'g', phi[int(nx/2), :], y, 'r')
        plt.axis([0, 0.5, 0, Ly])
        #'plt.gca().set_aspect(0.5)  # pbaspect([0.5 1 0.5])
        plt.title('Porosity profiles')
        plt.xlabel('X axis')

        # Subplot 3: Bulk Viscosity
        plt.subplot(2, 2, 3)
        plt.contourf(xx, yy, eta_phi)  # Transpose eta_phi to match MATLAB's behavior
        plt.axis('image')
        plt.colorbar()
        plt.title('Bulk Viscosity')

        # Subplot 4: Porosity
        plt.subplot(2, 2, 4)
        plt.pcolor(xx, yy, phi)  # Transpose phi to match MATLAB's behavior
        plt.axis('image')
        plt.colorbar()
        plt.title('Porosity')
        plt.xlabel('X axis')
        #plt.draw()
    # Save the status
    if it % nout == 0 or it==1:
       # np.save(f'step{it}', {'x': x, 'y': y, 'Pe': Pe, 'phi': phi, 'eta_phi': eta_phi}) # save npy file
       #filename = 'step' + str(it) + '.mat'
        filename = 'step{}.mat'.format(it)
        savemat(filename, {'x': x, 'y': y, 'Pe': Pe, 'phi': phi, 'eta_phi': eta_phi})
        if ploton==1:
            #plt.figure();plt.pcolor(x,y,phi.T);plt.axis("scaled");plt.colorbar()#plt.ion();plt.show();
            plt.figure()#%figsize=(10, 5)
            plt.subplot(1, 2, 1)
            #plt.contour(x, y, phi.T, cmap='viridis')
            plt.pcolor(x, y, phi.T,shading='auto', cmap='viridis')
            plt.axis("scaled")
            plt.colorbar()
            plt.title('step {}:({:5.2f}δt)'.format(it,time/tsc))
            plt.xlabel('X')
            plt.ylabel('Y')
            plt.subplot(1, 2, 2)
            plt.plot(phi[int(nx / 2), :], y, 'r')
            plt.title('Porosity Profile')
            plt.savefig(f'step{it}.jpg')
            plt.close()
    if it > 10 and iter == iterMax1 and resid > epsi:
        #print(f'Not converged at {iter:5d} iteration. Save and return now.')
        print('Not converged at {:5d} iteration. Save and return now.'.format(iter))
        filename = 'step{}.mat'.format(it)
        savemat(filename, {'x': x, 'y': y, 'Pe': Pe, 'phi': phi, 'eta_phi': eta_phi})
        #return
        sys.exit()
        # Check for convergence and exit if necessary
    if time>ttotal:
        print('The time reachs ttotal={:5.2f}({:5.2f} δt).Now exit'.format(time,time/tsc))
        sys.exit()


#telapsed = time() - tstart
#£print('Elapsed time is {:.4f} seconds'.format(telapsed))

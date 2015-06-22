import time
import sys
from pylab import *
sys.path.append('../..')
from numpad import *
from pdb import set_trace
from scipy.optimize import minimize
import re, string, os


def degrade(_adarray_):
    if isinstance(_adarray_, adarray):
        return _adarray_._value
    return _adarray_

def upgrade(_ndarray_):
    if isinstance(_ndarray_, np.ndarray):
        return array(_ndarray_)
    return _ndarray_

def genMesh(xin, yin, xout, yout, layer, grid):
    """
    generate 2D Ubend blockMesh for openFoam using inner and
    outer wall vertices coordinates.
	inputs ----
    	xin, yin:	inner wall coordinates	array
        xout, yout:	outer wall coordinates, array
        layer: 		radial partition, 	array
        grid:		each partition's grid, 	int array
    	outputs ----
        Xg, Yg: mesh coordinates
        write file blockMeshDict
    """
    # ======== INPUT CHECK ========
    assert(layer.size == grid.size+1)
    assert(layer[0]==1. and layer[-1]==0.)
    assert(all( np.diff(layer)<0 ))
    assert(all(grid>=1))
    assert(os.path.isfile('blockMeshDict_template'))

    f = open('blockMeshDict_template','r')
    fstr = ''.join( f.readlines() )
    f.close()
    
    # ======== COORDINATE ========
    X, Y = [], []			# block mesh
    Xg, Yg = [], []			# grid mesh
    
    for ii in layer:
        X.append( (xout-xin)*ii + xin )
        Y.append( (yout-yin)*ii + yin )

    X = np.vstack(X).transpose()
    Y = np.vstack(Y).transpose()
    Z = np.zeros(X.shape)    

    spacing = - np.diff(layer) / grid	# grid's radial spacing (normalized)
    spacing = np.hstack( [ np.repeat(spacing[ii], grid[ii]) for ii in range(spacing.size) ] )
    spacing = layer.copy()
    for ii in spacing:
        Xg.append( (xout-xin)*ii + xin )
        Yg.append( (yout-yin)*ii + yin )

    Xg = vstack(Xg).transpose()
    Yg = vstack(Yg).transpose()
   
    # Z-direction extrusion
    # Each row of X,Y,Z corresponds to a pizza slice
    # Each slice has 2*layer.size vertices
    X = np.hstack([X, X])
    Y = np.hstack([Y, Y])
    Z = np.hstack([Z, 0.1*ones(Z.shape)])
    
    XYZ = np.vstack([X, Y, Z])
    
    # ======== VERTEX ========
    vertex_str = ''
    for ii in range(XYZ.shape[1]):
        vertex_str += '('+array_str( XYZ[:,ii] )[1:-1]+')\n\t'
    pattern = '_VERTEX_'
    fstr = re.sub(pattern, vertex_str, fstr)
    
    # ======== BLOCKS ========
    endv = (xin.size-1)*2*layer.size    # first vertex of last slice
    # blocks by slice
    # simpleGrading depends on jj
    Hex = []
    Grade = []
    for ii in r_[0 : endv : 2*layer.size]:
        for jj in range(layer.size-1):
            hexij = np.array( [ii, ii+1, ii+1+layer.size*2, ii+layer.size*2, \
                            ii+layer.size, ii+1+layer.size, \
                            ii+layer.size*3+1, ii+layer.size*3] , dtype=int)
            hexij += jj
            Hex.append(hexij)
            Grade.append(jj)    # mesh refinement
    
    Hex_str = ''
    for ii in range(len(Hex)):
        Hex_str += 'hex ' + '(' + array_str(Hex[ii])[1:-1] + ')\n\t'
        Hex_str += '(' + str(grid[Grade[ii]]) + ' 1 1)\n\t'
        Hex_str += 'simpleGrading (1 1 1)\n\n\t'
    pattern = '_HEX_'
    fstr = re.sub(pattern, Hex_str, fstr)
    
    # ======== EDGE ========
    pass
    
    # ======== BOUNDARY ========
    # inlet
    inlet_str = ''
    for ii in range(layer.size-1):
        faces = np.array([ii, ii+1, ii+1+layer.size, ii+layer.size],dtype=int)
        inlet_str += '('+array_str(faces)[1:-1]+')\n\t\t'
    pattern = '_INLET_FACES_'
    fstr = re.sub(pattern, inlet_str, fstr)
    
    # outlet
    outlet_str = ''
    for ii in range(layer.size-1):
        faces = np.array([ii+1, ii, ii+layer.size, ii+layer.size+1],dtype=int)+endv
        outlet_str += '('+array_str(faces)[1:-1]+')\n\t\t'
    pattern = '_OUTLET_FACES_'
    fstr = re.sub(pattern, outlet_str, fstr)
    
    # wall_inner
    wallin_str = ''
    for ii in range(xin.size-1):
        faces = np.array( [layer.size*(2*ii+1) -1, layer.size*(2*ii+3) -1,
                        layer.size*(2*ii+4) -1, layer.size*(2*ii+2) -1], dtype=int)
        wallin_str += '('+array_str(faces)[1:-1]+')\n\t\t'
    pattern = '_WALL_IN_FACES_'
    fstr = re.sub(pattern, wallin_str, fstr)
    
    # wall_outer
    wallout_str = ''
    for ii in range(xin.size-1):
        faces = np.array( [layer.size*(2*ii), layer.size*(2*ii+1),
                        layer.size*(2*ii+3), layer.size*(2*ii+2)], dtype=int)
        wallout_str += '('+array_str(faces)[1:-1]+')\n\t\t'
    pattern = '_WALL_OUT_FACES_'
    fstr = re.sub(pattern, wallout_str, fstr)
    
    # ======== WRITE TO FILE ========
    f = open('blockMeshDict','w')
    f.write(fstr)
    f.close()
    return Xg, Yg

def showgrid(Xg, Yg):
    [plot(x,y,'black') for x,y in zip(Xg,Yg)]
    [plot(x,y,'black') for x,y in zip(Xg.transpose(),Yg.transpose())]
    axis('equal')
    show(block=True)

def MikeMesh(pert_xin, pert_yin, pert_xout, pert_yout):
    # grid count in streamwise direction
    n_1 = 15   # 1st straight section
    n_2 = 20   # 2nd straight section
    n_b = 20   # bend section
    n_s = 4    # sponge region
    
    # straight section lengths
    L_1 = 2.5
    L_2 = 2.5
    L_s = 1.

    # straight section streamwise grid size
    dx1 = L_1/n_1
    dx2 = L_2/n_2

    # duct width, bend inner and outer radius, grid angle
    Width = 1.
    rin  = .25
    rout = Width + rin
    theta = np.linspace(np.pi/2, -np.pi/2, n_b+1)

    # turn off or on sponge at exit, 1 = on, 0 = off
    sponge = 1
    
    # bend region, interior
    xin = np.cos(theta) * rin
    yin = np.sin(theta) * rin
    # boundary bend, exterior
    xout = np.cos(theta) * rout
    yout = np.sin(theta) * rout

    # sponge region
    Rsponge = 4.
    cratio_sponge = pow(Rsponge, 1./(n_s-1.))
    dstart_sponge = L_s*(1-cratio_sponge)/(1-pow(cratio_sponge,n_s))
    xin_test = [dstart_sponge]
    xin_test2 = [dstart_sponge]
    for i in range(0,n_s-1):
       xin_new = xin_test[i]*cratio_sponge
       xin_test = np.hstack( [xin_test, xin_new] )
       xin_test2 = np.hstack( [xin_test2, sum(xin_test)] )
    xin_test2 = -L_2 - xin_test2
    
    xin_sponge = xin_test2
    xout_sponge = xin_sponge
    yin_sponge = -rin*np.ones(n_s)
    yout_sponge = -rout*np.ones(n_s)

    # complete boundary
    if sponge == 1:
       # interior boundary
       xin = np.hstack( [np.linspace(-L_1, -dx1, n_1), xin, 
                         np.linspace(-dx2, -L_2, n_2), xin_sponge] )
       yin = np.hstack( [rin*np.ones(n_1), yin, -rin*np.ones(n_2), yin_sponge] )
       # exterior boundary
       xout = np.hstack( [np.linspace(-L_1, -dx1, n_1), xout, 
                          np.linspace(-dx2, -L_2, n_2), xout_sponge] )
       yout = np.hstack( [rout*np.ones(n_1), yout, -rout*np.ones(n_2), yout_sponge] )
    else:
       n_s = 0
       # interior boundary
       xin = np.hstack( [np.linspace(-L_1, -dx1, n_1), xin, 
                         np.linspace(-dx2, -L_2, n_2)] )
       yin = np.hstack( [rin*np.ones(n_1), yin, -rin*np.ones(n_2)] )
       # exterior boundary
       xout = np.hstack( [np.linspace(-L_1, -dx1, n_1), xout, np.linspace(-dx2, -L_2, n_2)] )
       yout = np.hstack( [rout*np.ones(n_1), yout, -rout*np.ones(n_2)] )

    xin  += degrade(pert_xin)
    yin  += degrade(pert_yin)
    xout += degrade(pert_xout)
    yout += degrade(pert_yout)
    xin, yin, xout, yout = upgrade(xin), upgrade(yin), upgrade(xout), upgrade(yout)

    # spanwise layers
    # n = full width number of cells (must be even for symmetric grading)
    # Rspan = overall grading ratio (max cell size to min cell size)
    # cratio = cell expansion ratio (adjacent cell size ratio)
    nspan = 30
    n_halfspan = nspan/2		
    Rspan = 10.
    cratio = pow(Rspan, 1./(n_halfspan-1.))
    dstart = Width/2 * (1-cratio)/(1-pow(cratio,n_halfspan))
    span_test = [dstart]
    span_test2 = [dstart]

    for i in range(0,int(n_halfspan)-1):
       span_new = span_test[i]*cratio
       span_test = np.hstack( [span_test, span_new] )
       span_test2 = np.hstack( [span_test2, sum(span_test)] )
    span_test2 = np.hstack( [0, span_test2] )
    layer_in = span_test2

    n_cells = nspan*(n_1 + n_2 + n_b + n_s)
    
    # layer parameterization, # of layers == layer.size-1
    ratio = np.logspace(.1,.5,10)
    layer_out = -layer_in[:-1][::-1] + 1.
    layer = np.hstack([layer_in, layer_out])[::-1]
    grid = np.array( np.ones(layer.size-1) ,dtype=int)

    Xg, Yg = genMesh(xin, yin, xout, yout, layer, grid)
    Xg = Xg[:,::-1]
    Yg = Yg[:,::-1]
    # showgrid(Xg, Yg)
    return Xg, Yg, xin, yin, xout, yout

def visualize(w, geo, state):
    # Visualize Mach number, non-dimensionalized stagnation and static pressure
    def avg(a):
        return 0.25 * (a[1:,1:] + a[1:,:-1] + a[:-1,1:] + a[:-1,:-1])

    import numpy as np
    rho, u, v, E, p = primative(value(extend(w, geo, state)), state)
    rho, u, v, E, p = degrade(rho), degrade(u), degrade(v), degrade(E), degrade(p)

    x, y = value(geo.xy)
    xc, yc = value(geo.xyc)
    
    c2 = 1.4 * p / rho
    M = sqrt((u**2 + v**2) / c2)
    pt = p * (1 + 0.2 * M**2)**3.5

    subplot(2,2,1)
    contourf(x, y, degrade(avg(M)), 100)
    colorbar()
    quiver(xc, yc, u[1:-1,1:-1], v[1:-1,1:-1])
    axis('scaled')
    xlabel('x')
    ylabel('y')
    title('Mach')
    draw()

    subplot(2,2,2)
    pt_frac = (pt - p_out) / (pt_in - p_out)
    contourf(x, y, degrade(avg(pt_frac)), 100)
    colorbar()
    axis('scaled')
    xlabel('x')
    ylabel('y')
    title('pt')
    draw()
    
    subplot(2,2,3)
    p_frac = (p - p_out) / (pt_in - p_out)
    contourf(x, y, degrade(avg(p_frac)), 100)
    colorbar()
    axis('scaled')
    xlabel('x')
    ylabel('y')
    title('p')
    draw()

    subplot(2,2,4)
    contourf(x, y, degrade(avg(rho)), 100)
    colorbar()
    axis('scaled')
    xlabel('x')
    ylabel('y')
    title('density')
    draw()
    show(block=True)


def visual_compare(w, geo, state, standard_sol):
    # Visualize Mach number, non-dimensionalized stagnation and static pressure
    # compare with true state equation's solution
    def avg(a):
        return 0.25 * (a[1:,1:] + a[1:,:-1] + a[:-1,1:] + a[:-1,:-1])

    import numpy as np
    rho, u, v, E, p = primative(value(extend(w, geo, state)), state)
    rho, u, v, E, p = degrade(rho), degrade(u), degrade(v), degrade(E), degrade(p)
    x, y = value(geo.xy)
    xc, yc = value(geo.xyc)
    rho_s, u_s, v_s, E_s = standard_sol[0], standard_sol[1], standard_sol[2], standard_sol[3]
   
    Vel = sqrt(u**2 + v**2)
    Vel_s = sqrt(u_s**2 + v_s**2)

    subplot(2,2,1)
    contourf(x, y, degrade(abs(avg(Vel-Vel_s))), 100)
    colorbar()
    quiver(xc, yc, (u-u_s)[1:-1,1:-1], (v-v_s)[1:-1,1:-1])
    axis('scaled')
    xlabel('x')
    ylabel('y')
    title('Velocity')
    draw()

    subplot(2,2,2)
    contourf(x, y, degrade(abs(avg(E-E_s))), 100)
    colorbar()
    axis('scaled')
    xlabel('x')
    ylabel('y')
    title('Energy')
    draw()
    
    subplot(2,2,3)
    contourf(x, y, degrade(abs(avg(rho-rho_s))), 100)
    colorbar()
    axis('scaled')
    xlabel('x')
    ylabel('y')
    title('density')
    draw()
    show(block=True)

def extend(w_interior, geo, state):
    '''
    Extend the conservative variables into ghost cells using boundary condition
    '''
    w = zeros([4, Ni+2, Nj+2])
    w[:,1:-1,1:-1] = w_interior.reshape([4, Ni, Nj])

    # inlet
    rho, u, v, E, p = primative(w[:,1,1:-1], state)
    c2 = 1.4 * p / rho
    c = sqrt(c2)
    mach2 = (u**2+v**2) / c2
    rhot = rho * (1 + 0.2 * mach2)**2.5
    pt = p * (1 + 0.2 * mach2)**3.5

    d_rho = 1 - rho
    d_pt = pt_in - pt
    d_u = d_pt / (rho * (u + c))
    d_p = rho * c * d_u

    rho = rho + d_rho
    u = u + d_u
    p = p + d_p
    w[0,0,1:-1] = rho
    w[1,0,1:-1] = rho * u
    w[2,0,1:-1] = 0
    inlet_U = state.U(p,rho,'inlet')
    if isinstance(inlet_U, bool):
        return False
    else:
        w[3,0,1:-1] = inlet_U + .5*rho*(u**2+v**2)

    # outlet
    w[:,-1,1:-1] = w[:,-2,1:-1]
    rho, u, v, E, p = primative(w[:,-1,1:-1], state)
    p = p_out*ones(rho.shape)
    outlet_U = state.U(p,rho,'outlet')
    if isinstance(outlet_U, bool):
        return False
    else:
        w[3,-1,1:-1] = outlet_U + .5*rho*(u**2+v**2)

    # walls
    w[:,:,0] = w[:,:,1]
    rhoU_n = sum(w[1:3,1:-1,0] * geo.normal_j[:,:,0], 0)
    w[1:3,1:-1,0] -= 2 * rhoU_n * geo.normal_j[:,:,0]

    w[:,:,-1] = w[:,:,-2]
    rhoU_n = sum(w[1:3,1:-1,-1] * geo.normal_j[:,:,-1], 0)
    w[1:3,1:-1,-1] -= 2 * rhoU_n * geo.normal_j[:,:,-1]

    return w
    
def primative(w, state):
    '''
    Transform conservative variables into primative ones
    '''
    rho = w[0]
    u = w[1] / rho
    v = w[2] / rho
    E = w[3]
    p = state.P( E - .5*(u*w[1]+v*w[2]), rho)
    return rho, u, v, E, p

def grad_dual(phi, geo):
    '''
    Gradient on the nodes assuming zero boundary conditions
    '''
    dxy_i = 0.5 * (geo.dxy_i[:,1:,:] + geo.dxy_i[:,:-1,:])
    dxy_j = 0.5 * (geo.dxy_j[:,:,1:] + geo.dxy_j[:,:,:-1])
    phi_i = array([phi * dxy_i[1], -phi * dxy_i[0]])
    phi_j = array([-phi * dxy_j[1], phi * dxy_j[0]])

    grad_phi = zeros(geo.xy.shape)
    grad_phi[:,:-1,:-1] += phi_i + phi_j
    grad_phi[:,1:,:-1] += -phi_i + phi_j
    grad_phi[:,:-1,1:] += phi_i - phi_j
    grad_phi[:,1:,1:] += -phi_i - phi_j

    area = zeros(geo.xy.shape[1:])
    area[:-1,:-1] += geo.area
    area[1:,:-1] += geo.area
    area[:-1,1:] += geo.area
    area[1:,1:] += geo.area

    return 2 * grad_phi / area

def ns_flux(rho, u, v, E, p, grad_u):
    # viscous stress
    dudx, dudy, dvdx, dvdy = grad_u
    mu_t = 0 # blabla
    sigma_xx = (mu + mu_t) * (2 * dudx - 2./3 * (dudx + dvdy))
    sigma_yy = (mu + mu_t) * (2 * dvdy - 2./3 * (dudx + dvdy))
    sigma_xy = (mu + mu_t) * (dudy + dvdx)

    F = array([rho * u, rho * u**2 + p - sigma_xx,
                        rho * u * v    - sigma_xy,
                        u * (E + p)    - sigma_xx * u - sigma_xy * v])
    G = array([rho * v, rho * u * v    - sigma_xy,
                        rho * v**2 + p - sigma_yy,
                        v * (E + p)    - sigma_xy * u - sigma_yy * v])
    return F, G

def sponge_flux(c_ext, w_ext, geo):
    ci = 0.5 * (c_ext[1:,1:-1] + c_ext[:-1,1:-1])
    cj = 0.5 * (c_ext[1:-1,1:] + c_ext[1:-1,:-1])

    a = geo.area
    ai = vstack([a[:1,:], (a[1:,:] + a[:-1,:]) / 2, a[-1:,:]])
    aj = hstack([a[:,:1], (a[:,1:] + a[:,:-1]) / 2, a[:,-1:]])

    wxx = (w_ext[:,2:,1:-1] + w_ext[:,:-2,1:-1] - 2 * w_ext[:,1:-1,1:-1]) / 3.
    wyy = (w_ext[:,1:-1,2:] + w_ext[:,1:-1,:-2] - 2 * w_ext[:,1:-1,1:-1]) / 3.
    # second order dissipation at boundary, fourth order in the interior
    Fi = -0.5 * ci * ai * (w_ext[:,1:,1:-1] - w_ext[:,:-1,1:-1])
    Fi[:,1:-1,:] = 0.5 * (ci * ai)[1:-1,:] * (wxx[:,1:,:] - wxx[:,:-1,:])
    Fj = -0.5 * cj * aj * (w_ext[:,1:-1,1:] - w_ext[:,1:-1,:-1])
    Fj[:,:,1:-1] = 0.5 * (cj * aj)[:,1:-1] * (wyy[:,:,1:] - wyy[:,:,:-1])
    return Fi, Fj

def ns_kec(w, w0, geo, dt, state):
    '''
    Kinetic energy conserving scheme with no numerical viscosity
    '''
    w_ext = extend(w, geo, state)
    # check if inlet/outlet U solve blow up
    if isinstance(w_ext, bool):
        return ones(w0.shape)*1e5
    rho, u, v, E, p = primative(w_ext, state)
    c = sqrt(1.4 * p / rho)
    # velocity gradient on nodes
    dudx, dudy = grad_dual(u[1:-1,1:-1], geo)
    dvdx, dvdy = grad_dual(v[1:-1,1:-1], geo)
    duv_dxy = array([dudx, dudy, dvdx, dvdy])
    # interface average
    rho_i = 0.5 * (rho[1:,1:-1] + rho[:-1,1:-1])
    rho_j = 0.5 * (rho[1:-1,1:] + rho[1:-1,:-1])
    u_i = 0.5 * (u[1:,1:-1] + u[:-1,1:-1])
    u_j = 0.5 * (u[1:-1,1:] + u[1:-1,:-1])
    v_i = 0.5 * (v[1:,1:-1] + v[:-1,1:-1])
    v_j = 0.5 * (v[1:-1,1:] + v[1:-1,:-1])
    E_i = 0.5 * (E[1:,1:-1] + E[:-1,1:-1])
    E_j = 0.5 * (E[1:-1,1:] + E[1:-1,:-1])
    p_i = 0.5 * (p[1:,1:-1] + p[:-1,1:-1])
    p_j = 0.5 * (p[1:-1,1:] + p[1:-1,:-1])
    # interface strain rate averged from dual mesh (nodal) values
    duv_dxy_i = 0.5 * (duv_dxy[:,:,1:] + duv_dxy[:,:,:-1])
    duv_dxy_j = 0.5 * (duv_dxy[:,1:,:] + duv_dxy[:,:-1,:])
    # inlet and outlet have no viscous stress
    duv_dxy_i[:,[0,-1],:] = 0
    # interface flux
    f_i, g_i = ns_flux(rho_i, u_i, v_i, E_i, p_i, duv_dxy_i)
    f_j, g_j = ns_flux(rho_j, u_j, v_j, E_j, p_j, duv_dxy_j)
    Fi = + f_i * geo.dxy_i[1] - g_i * geo.dxy_i[0]
    Fj = - f_j * geo.dxy_j[1] + g_j * geo.dxy_j[0]
    # sponge
    Fi_s, Fj_s = sponge_flux(c, w_ext, geo)
    Fi += 0.5 * Fi_s
    Fj += 0.5 * Fj_s
    # residual
    divF = (Fi[:,1:,:] - Fi[:,:-1,:] + Fj[:,:,1:] - Fj[:,:,:-1]) / geo.area
    # E.diff(w) gives dense matrix, but E[i].diff(w) gives sparse matrix
    return (w - w0) / dt + ravel(divF)


# -------------------------- geometry ------------------------- #
class geo2d:
    def __init__(self, xy):
        xy = array(xy)
        self.xy = xy
        self.xyc = (xy[:,1:,1:]  + xy[:,:-1,1:] + \
                    xy[:,1:,:-1] + xy[:,:-1,:-1]) / 4

        self.dxy_i = xy[:,:,1:] - xy[:,:,:-1]
        self.dxy_j = xy[:,1:,:] - xy[:,:-1,:]

        self.L_j = sqrt(self.dxy_j[0]**2 + self.dxy_j[1]**2)
        self.normal_j = array([self.dxy_j[1] / self.L_j,
                              -self.dxy_j[0] / self.L_j])

        self.area = self.tri_area(self.dxy_i[:,:-1,:], self.dxy_j[:,:,1:]) \
                  + self.tri_area(self.dxy_i[:,1:,:], self.dxy_j[:,:,:-1]) \

    def tri_area(self, xy0, xy1):
        return 0.5 * (xy0[1] * xy1[0] - xy0[0] * xy1[1])


# ----------------------- flow solver ------------------------- #
def flowsolver(xy, state):
    geo = geo2d(xy)
    t, dt = 0, 1./Nj/100.

    w0 = zeros([4, Ni, Nj])
    w0[0] = 1
    w0[3] = 1E5 / (1.4 - 1)
    w0 = ravel(w0)

    for i in range(100):
        w = solve(ns_kec, w0, args=(w0, geo, dt, state), rel_tol=1E-5, abs_tol=1E-4, max_iter=15,
                  verbose=False)
        print('i = ', i, 't = ', format(t,'.4f'), 'res = ', format(w._res_norm,'.2f'))
        if isinf(w._n_Newton) or isnan(w._res_norm):
            set_trace()
            return None

        if w._n_Newton == 1 and t>500.:
            w0 = w
            w0.obliviate()
            dt *= 2
            break
        elif w._n_Newton < 5:
            w0 = w
            dt *= 2
        elif w._n_Newton < 7:
            w0 = w
        else:
            dt *= 0.5
            continue
        t += dt
        w0.obliviate()

    dt = np.inf
    w = solve(ns_kec, w0, args=(w0, geo, dt, state), rel_tol=1E-5, abs_tol=1E-5, max_iter=30,
              verbose=False)

    # visualize(w, geo, state)
    # visual_compare(w, geo, state, standard_sol)
    return extend(w, geo, state)

'p = p(U,rho)'
class StateEqn:
   
    def __init__(self, Un, Rn, Urange, Rrange):
        self.Un, self.Rn = Un, Rn
        self.Uis = np.linspace(Urange[0], Urange[1], Un)
        self.Ris = np.linspace(Rrange[0], Rrange[1], Rn)
        self.Ubeta = 1.2*Un/np.diff(Urange)
        self.Rbeta = 1.2*Rn/np.diff(Rrange)
        self.coef = None
        self._last_inlet_U = None
        self._last_outlet_U = None

    def setcoef(self, coef):
        assert(coef.shape == (self.Un, self.Rn,))
        self.coef = upgrade(coef)

    def P(self, U, R):
        assert(self.coef is not None)
        assert(U.shape == R.shape)
        result = zeros(U.shape)
        for i in range(self.Un):
            for j in range(self.Rn):
                result += sigmoid(self.Ubeta*(U-self.Uis[i])) * \
                          sigmoid(self.Rbeta*(R-self.Ris[j])) * \
                          self.coef[i,j]
        result += 0.04*U + 4e4*R
        return result

    def _P_res_(self, U, P, R):
        return self.P(U, R) - P

    def U(self, P, R, which_boundary):
        assert(self.coef is not None)
        assert(P.shape == R.shape)
        if which_boundary == 'inlet':
            if self._last_inlet_U is None:
                for relax_factor in [.75, .5, .25]:
                    result = solve(self._P_res_, np.zeros(degrade(R).shape), 
                             args=(P,R), rel_tol=1e-7, abs_tol=1e-6, 
                             max_iter=30, verbose=False, lam=relax_factor)
                    if not isinf(result._n_Newton):
                        break
            else:
                for relax_factor in [1., .75, .5]:
                    result = solve(self._P_res_, self._last_inlet_U,
                             args=(P,R), rel_tol=1e-7, abs_tol=1e-6, 
                             max_iter=20, verbose=False, lam=relax_factor)
                    if not isinf(result._n_Newton):
                        break
            if result._res_norm < 1.:
                self._last_inlet_U = result.copy()
                self._last_inlet_U.obliviate()

        if which_boundary == 'outlet':
            if self._last_outlet_U is None:
                for relax_factor in [.5, .35, .2]:
                    result = solve(self._P_res_, np.zeros(degrade(R).shape), 
                             args=(P,R), rel_tol=1e-7, abs_tol=1e-6, 
                             max_iter=100, verbose=False, lam=relax_factor)
                    if not isinf(result._n_Newton):
                        break
            else:
                for relax_factor in [1., .75, .5]:
                    result = solve(self._P_res_, self._last_outlet_U, 
                             args=(P,R), rel_tol=1e-7, abs_tol=1e-6, 
                             max_iter=30, verbose=False, lam=relax_factor)
                    if not isinf(result._n_Newton):
                       break
            if result._res_norm < 1.:
                self._last_outlet_U = result.copy()
                self._last_outlet_U.obliviate()

        # print(' '*7, 'U.resid: ', result._res_norm)
        if result._res_norm > 1.:
            return False
        return result

    def clean(self):
        self._last_inlet_U = None
        self._last_outlet_U = None


# ----------------- normal direction of geo at pert_section ---------------------- #
def normal_geo(geo):
    normal = geo.dxy_i[:,pert_section,0]._value
    L = np.sqrt(normal[0]**2 + normal[1]**2)
    normal[0] /= L
    normal[1] /= L
    return degrade(normal)


# ---------- convert normal perturbation to xy perturbation at pert_section ------ #
def normal_to_xy(pertin, pertout, geo):
    normal = normal_geo(geo)
    return pertin*normal, pertout*normal


# ----------------- objective and gradient--------------------- #
def objective(w, geo_base):
    # total outlet mass flux
    w = w.reshape([4,Ni,Nj])
    outflux = sum( w[1][-1,:] * np.diff(geo_base.xy[1][-1,:]) )
    # nodal gradient
    ginx = asarray(outflux.diff(xin))[0]
    giny = asarray(outflux.diff(yin))[0]
    goutx = asarray(outflux.diff(xout))[0]
    gouty = asarray(outflux.diff(yout))[0]
    # normal direction
    normal = normal_geo(geo_base)
    # projected gradient
    gradin = ginx[pert_section] * normal[0] + giny[pert_section] * normal[1]
    gradout = goutx[pert_section] * normal[0] + gouty[pert_section] * normal[1]
    return outflux, gradin, gradout


# ----------------- solution mismatch and gradient ---------------- #
def mismatch_var_grad(coef, xy, state, standard_sol):
    state._last_inlet_U = None
    state._last_outlet_U = None
    state.setcoef(coef.reshape([NE,NR])*1e3)
    w_ext = flowsolver(xy, state)
    var, grad = mismatch(w_ext, state, standard_sol)
    return var, grad
    

def mismatch(w_ext, state, standard_sol):
    global last_valid_state
    if w_ext is None:
        return 1e3, 1e2 * (state.coef/1e3 - last_valid_state)._value 

    w_ext = w_ext.reshape([4, Ni+2, Nj+2])
    rho, u, v, E = primative(w_ext, state)[:-1]

    rho_s, u_s, v_s, E_s = standard_sol[0], standard_sol[1], standard_sol[2], standard_sol[3]

    weights = np.array( [rho_s.max() - rho_s.min(),
                         u_s.max()-u_s.min(),
                         v_s.max()-v_s.min(),
                         E_s.max()-E_s.min()] )
               
    diffs = array( [linalg.norm(rho-rho_s,2),
                    linalg.norm(u-u_s,2),
                    linalg.norm(v-v_s,2),
                    linalg.norm(E-E_s,2)] )

    err = sum(diffs/weights) + linalg.norm(state.coef, L=1)*1e-5
    grad = asarray( err.diff(state.coef) ) * 1e3
 
    print('coef', state.coef)
    print('grad', grad)
    print('err', err)

    w_ext.obliviate()
    err.obliviate()
    state.coef.obliviate()
    last_valid_state = state.coef._value.copy() / 1e3
    return err._value, grad
        

if __name__ == '__main__':
    # global variables
    Ni, Nj = 59, 30
    pt_in = 1.2E5
    p_out = 1E5
    mu = 1
    last_valid_state = None
    iter_count = 0
    pert_section = np.r_[16:35]
    
    # base geometry
    x0, y0, xin0, yin0, xout0, yout0 = \
        MikeMesh(np.zeros(Ni+1), np.zeros(Ni+1), np.zeros(Ni+1), np.zeros(Ni+1))
    geo_base = geo2d(array([x0,y0]))

    # add perturbation to base geometry
    # gradin = np.loadtxt('gradin')
    # gradout = np.loadtxt('gradout')

    # pertin, pertout = \
    #     normal_to_xy(np.ones_like(gradin)*1e-5, np.ones_like(gradout)*1e-5, geo_base)

    pert_xin = np.zeros(Ni+1)
    pert_yin = np.zeros(Ni+1)
    # pert_xin[pert_section] += pertin[0]
    # pert_yin[pert_section] += pertin[1]

    pert_xout = np.zeros(Ni+1)
    pert_yout = np.zeros(Ni+1)
    # pert_xout[pert_section] += pertout[0]
    # pert_yout[pert_section] += pertout[1]

    # solve flow
    x, y, xin, yin, xout, yout = MikeMesh(pert_xin, pert_yin, pert_xout, pert_yout)
    NE, NR = 10, 10
    state = StateEqn(NE, NR,[0e5,4e5],[0.,1])

    # coef = np.ravel( np.ones([NE,NR])*1.5 )
    coef = np.loadtxt('sol/checkpoint/statecoef1')
    
    bnds = tuple( map(tuple, np.tile( np.array([0.,1e2]), [coef.size,1] )) )
    case = 0
    DIR = '/home/voila/Documents/2014GRAD/numpad/test/sol/standard_sol/'
    standard_sol = [np.loadtxt(DIR+'rho'+str(case)),
                    np.loadtxt(DIR+'u'+str(case)),
                    np.loadtxt(DIR+'v'+str(case)),
                    np.loadtxt(DIR+'E'+str(case))]
    
    coef = minimize(mismatch_var_grad, coef, args=(array([x,y]), state, standard_sol), 
                    method='SLSQP', bounds=bnds, jac=True)
    set_trace()
    


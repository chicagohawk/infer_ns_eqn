def vis(w, geo, coef):
    # Visualize Mach number, non-dimensionalized stagnation and static pressure
    def avg(a):
        return 0.25 * (a[1:,1:] + a[1:,:-1] + a[:-1,1:] + a[:-1,:-1])

    import numpy as np
    rho, u, v, E, p = primative(value(extend(w, geo, coef)), coef)
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
    contourf(x, y, avg(pt_frac), 100)
    colorbar()
    axis('scaled')
    xlabel('x')
    ylabel('y')
    title('pt')
    draw()
    
    subplot(2,2,3)
    p_frac = (p - p_out) / (pt_in - p_out)
    contourf(x, y, avg(p_frac), 100)
    colorbar()
    axis('scaled')
    xlabel('x')
    ylabel('y')
    title('p')
    draw()

'''
def vis(w, geo, coef):
    # Visualize Mach number, non-dimensionalized stagnation and static pressure
    def avg(a):
        return 0.25 * (a[1:,1:] + a[1:,:-1] + a[:-1,1:] + a[:-1,:-1])

    import numpy as np
    rho, u, v, E, p = primative(value(extend(w, geo, coef)), coef)
    rho, u, v, E, p = degrade(rho), degrade(u), degrade(v), degrade(E), degrade(p)

    x, y = value(geo.xy)
    xc, yc = value(geo.xyc)
    
    c2 = 1.4 * p / rho
    M = sqrt((u**2 + v**2) / c2)
    pt = p * (1 + 0.2 * M**2)**3.5

    avgM_s     = np.loadtxt('sol/vis_avgM_s')
    u_s        = np.loadtxt('sol/vis_u_s')
    v_s        = np.loadtxt('sol/vis_v_s')
    avgpt_frac = np.loadtxt('sol/vis_avg_pt_frac')
    avgp_frac  = np.loadtxt('sol/vis_avg_p_frac')

    subplot(2,2,1)
    contourf(x, y, degrade(avg(M)) - avgM_s, 100)
    colorbar()
    quiver(xc, yc, u[1:-1,1:-1] - u_s, v[1:-1,1:-1] - v_s)
    axis('scaled')
    xlabel('x')
    ylabel('y')
    title('Mach')
    draw()
    
    subplot(2,2,2)
    pt_frac = (pt - p_out) / (pt_in - p_out)
    contourf(x, y, avg(pt_frac) - avgpt_frac, 100)
    colorbar()
    axis('scaled')
    xlabel('x')
    ylabel('y')
    title('pt')
    draw()
    
    subplot(2,2,3)
    p_frac = (p - p_out) / (pt_in - p_out)
    contourf(x, y, avg(p_frac) - avgp_frac, 100)
    colorbar()
    axis('scaled')
    xlabel('x')
    ylabel('y')
    title('p')
    draw()

'''




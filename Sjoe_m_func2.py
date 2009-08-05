from math import *
from scipy import *
from pylab import *

## Physical constants
g = 9.82; # Gravity constant
c_el = 1.065e10;
c_sq1 = 0.08;
c_sq2 = 5.e-4;
c_ro1 = 2.1e-6;
c_ro2 = 5.e-4;
mu0 = 0.1;
gamma = 8.e3;

## Misc constants
a0 = 5.e-3;  # Radius of rolling element (RE)
b0 = 30.e-3; # Radius of outer ring
c0 = 20.e-3; # Radius of inner ring
#d0 = 4.e-5; # Offset between the origins of the ring coordinate systems
d0 = 0; # Offset between the origins of the ring coordinate systems

rho = 7.82e3; # Density of the RE (iron)
m = 4.*pi/3.*a0**3*rho; # Mass of the RE (it's a sphere)
I = 2.*m*a0**2/5.; # Inertia of the RE
w0 = 628.32; # Angular speed of the outer ring (~6000 rpm)
#w0 = 52.36; # Angular speed of the outer ring (~500 rpm)
m0 = 1 # Mass of the outer ring
F_a = 1000 # External force on the outer ring
#F_a = 0


# This function calculates the right hand side of the equations of motion for
# rolling element j and the forces acting on the outer ring.
def rhs_ball(ball_state, outer_ring_state):
	r, th, rdot, thdot, phidot = ball_state
	x, y, xdot, ydot = outer_ring_state
# Update the deltas
	ro = sqrt((r*sin(th)-y)**2 + (r*cos(th)-x)**2);
	tho = atan2(r*sin(th)-y, r*cos(th)-x);
	rdoto = (r*rdot - (r*ydot+rdot*y)*sin(th) - (r*xdot+rdot*x)*cos(th) + y*ydot + x*xdot + r*thdot*(x*sin(th)-y*cos(th)))/ro;
	thdoto = ((r*xdot-rdot*x)*sin(th) + (rdot*y-r*ydot)*cos(th) - r*thdot*(x*cos(th)+y*sin(th)) + x*ydot - xdot*y + r*r*thdot)/ro**2;

	di = (r-a0)-c0;
	do = b0-(ro+a0);

	if di > 0:
		ai = a0;
	else:
		ai = a0+di/2.;
	
	if do > 0:
		ao = a0;
	else:
		ao = a0+do/2.;
	
	if di > 0:
		c = c0;
	else:
		c = c0+di/2.;

	if do > 0:
		b = b0;
	else:
		b = b0-do/2.;

	drdoti = rdot;
	drdoto = -rdoto;

	dthdoti = r*thdot - ai*phidot;
	dthdoto = ro*thdoto + ao*phidot - b*w0;

	si = (a0-c0)*thdot - ai*phidot;
	so = -(a0+b0)*thdoto + ao*phidot + b*w0;

# Update forces and torques
	if di <= 0:
		F_eli = c_el*abs(di)**(3./2.);
	else:
		F_eli = 0;
	
	if do <= 0:
		F_elo = c_el*abs(do)**(3./2.);
	else:
		F_elo = 0;

	F_sqi = -F_eli*c_sq1*atan(drdoti/c_sq2)*2/pi;
	F_sqo = -F_elo*c_sq1*atan(drdoto/c_sq2)*2/pi;
	
	Ni = F_eli + F_sqi; # Normal force w.r.t. inner ring
	No = F_elo + F_sqo; # Normal force w.r.t. outer ring
	
	mu_i = 2.*mu0/pi*atan(gamma*dthdoti*pi/2./mu0); # Coefficient of friction
	mu_o = 2.*mu0/pi*atan(gamma*dthdoto*pi/2./mu0); # Coefficient of friction   
	
	F_sli = abs(Ni)*mu_i;
	F_slo = abs(No)*mu_o;

	M_sli = ai*F_sli;
	M_slo = ao*F_slo;
	
	M_rmi = abs(Ni)*c_ro1*atan(si/c_ro2)*2./pi;
	M_rmo = abs(No)*c_ro1*atan(so/c_ro2)*2./pi;

# Calculate relevant forces and transform to global coordinate system #

	F_ri = F_eli + F_sqi;
	F_ro = -F_elo - F_sqo;

	F_thi = -F_sli;
	F_tho = -F_slo;
	
	Mi = M_sli - M_rmi;
	Mo = -M_slo + M_rmo;

# Transform from outer ring coordinate system to global
	F_ro_t = F_ro*cos(th-tho) + F_tho*sin(th-tho);
	F_tho_t = -F_ro*sin(th-tho) + F_tho*cos(th-tho);
	F_ro = F_ro_t;
	F_tho = F_tho_t;
	
	F_r = F_ri + F_ro;
	F_th = F_thi + F_tho;
	M = Mi + Mo;

# Transform the forces acting on the outer ring to ones in x,y-directions
	F_x = -F_ro*cos(tho) + F_tho*sin(tho)
	F_y = -F_ro*sin(tho) - F_tho*cos(tho)   
	
# Finally, find the right-hand side for the ball
	rhs = array(zeros(3)); # Right hand side of equations of motion
	rhs[0] = (F_r - m*g*sin(th) + m*r*thdot**2)/m;
	rhs[1] = (F_th - m*g*cos(th) - 2.*m*rdot*thdot)/r/m;
	rhs[2] = M/I;
	
	return rhs, F_x, F_y

# This function calculates the right hand side of the equations of motion
def oderhs(z,t):
	n = (len(z)-4)/5
	dz = array(zeros(len(z)))
	ring_state = r_[z[0:2], z[2+2*n:4+2*n]]
	F_x = 0
	F_y = F_a

	for j in xrange(0,n):
		rhs, F_xj, F_yj = rhs_ball(r_[z[2+2*j:4+2*j], z[4+2*n+3*j:7+2*n+3*j]], ring_state)
		F_x += F_xj
		F_y += F_yj
		dz[2+3*j] = z[2+2*j]
		dz[3+3*j] = z[3+2*j]

		dz[4+2*n+3*j:7+2*n+3*j] = rhs
	
	dz[0:2] = z[2+2*n:4+2*n]
	dz[2+2*n] = F_x/m0
	dz[3+2*n] = F_y/m0
	return dz
from math import *
from scipy import *
from pylab import *
import lsodar

#    ## Deltas and deformations ##
#    di=0; # Intersection between RE and inner ring
#    do=0; # Intersection between RE and outer ring
#    ai=0; # Deformed  radius of RE w.r.t. inner ring
#    ao=0; # Deformed  radius of RE w.r.t. outer ring
#    b=0;  # Deformed  radius of outer ring
#    c=0;  # Deformed  radius of inner ring

#    drdoti=0;  # Derivative of di
#    drdoto=0;  # Derivative of do
#    dthdoti=0; # Relative velocity in tangential direction between RE and inner ring
#    dthdoto=0; # Relative velocity in tangential direction between RE and outer ring

#    si=0; # Sum of surface velocities in tangential direction between RE and inner ring
#    so=0; # Sum of surface velocities in tangential direction between RE and outer ring

#    ## Forces and torques
#    F_eli=0; # Elastic Hertz contact force w.r.t. inner ring
#    F_elo=0; # Elastic Hertz contact force w.r.t. outer ring

#    F_sqi=0; # Squeeze material damping force w.r.t. inner ring
#    F_sqo=0; # Squeeze material damping force w.r.t. outer ring

#    F_sli=0; # Slip friction force w.r.t. inner ring
#    F_slo=0; # Slip friction force w.r.t. outer ring

#    M_sli=0; # Torque from slip force w.r.t. inner ring
#    M_slo=0; # Torque from slip force w.r.t. outer ring

#    M_rmi=0;  # Rolling material damping torque w.r.t. inner ring
#    M_rmo=0;  # Rolling material damping torque w.r.t. outer ring

## Physical constants
g = 9.82*0; # Gravity constant
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
d0 = 0*4.e-5; # Offset between the origins of the ring coordinate systems

rho = 7.82e3; # Density of the RE (iron)
m = 4.*pi/3.*a0**3*rho; # Mass of the RE (it's a sphere)
I = 2.*m*a0**2/5.; # Inertia of the RE
w0 = 0*52.36; # Angular speed of the outer ring (~500 rpm)


# This function calculates the right hand side of the equations of motion:
# [d2r/dt2; d2theta/dt2; d2phi/dt2] = rhs(r, dr/dt, theta, dtheta/dt, phi,
# dphi/dt)
def update_rhs(r, th, rdot, thdot, phidot):

# Update the deltas
	ro = sqrt(r**2 + d0**2 -2.*d0*r*sin(th));
	tho = atan2(r*sin(th)-d0, r*cos(th));
	rdoto = ((r - d0*sin(th))*rdot - d0*r*thdot*cos(th))/ro;
	thdoto = (r**2*thdot + d0*rdot*cos(th) - d0*r*thdot*sin(th))/ro**2;

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
	
# Finally, find the right-hand side
	rhs = array(zeros(3)); # Right hand side of equations of motion
	rhs[0] = (F_r - m*g*sin(th) + m*r*thdot**2)/m;
	rhs[1] = (F_th - m*g*cos(th) - 2.*m*rdot*thdot)/r/m;
	rhs[2] = M/I;
	
	return rhs

def oderhs(y,t):
	r = y[0]
	th = y[1]
	rdot = y[2]
	thdot = y[3]
	phidot = y[4]
	rhs = update_rhs(r, th, rdot, thdot, phidot)

	dy = zeros(5)

	dy[0] = y[2];
	dy[1] = y[3];

	dy[2] = rhs[0];
	dy[3] = rhs[1];
	dy[4] = rhs[2];
	return dy
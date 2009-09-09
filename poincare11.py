from scipy import *
from scipy.integrate import *
from numpy.linalg import norm, lstsq
from exceptions import *
import lsodar


def poincare(f, t0, y0, dist = None, transform = None, inv_transform=None, filename = None):

    if dist == None:
        dist = lambda y1,y2: norm(y1-y2)
    if transform == None:
        transform = lambda y,t: y
    if inv_transform == None:
        inv_transform = lambda y,t: y

    # Point to fix the Poincare section in space
    y0t = transform(y0,t0)

    # Normal vector of the Poincare section:
    #n = f(y0t, t0)
    #n = n/norm(n)
    n = zeros(len(y0t))
    n[1] = 1
    print "Hyperplane normal vector:", n

    # Solution is in the Poincare section if (y(t), n) = C
    C = dot(y0t,n)

    # "x = y" if abs(x-y)< tol
    tol = 1e-2
    rtol=1e-4
    atol=1e-4

    # Maximum number of iterations
    maxit = 1000
    # Will be the difference between two iterates x_j and x_{j-1}
    delta = 1000000

    # Will hold the iterates
    x = [copy(y0)]
    xt = [copy(transform(y0,t0))] # Transformed iterates
    T = [copy(t0)] # Iterate return times
    d  = [delta]   # Deltas

    def event(y,t):
        ym = copy(transform(y,t))
        return dot(ym,n) - C

    its = 0                 # Number of iterations.
    MPE_count=0             # Number of MPE approximations.
    max_t = 1             # Maximum time to integrate before stopping.
    t_small = 1e-3          # Time to integrate without checking if we're in the.
                            # Poincare section. If =0 the integrator might stop at the.
                            # ~same point we started at.
    p_it = 1                # Number of Poincare steps to take for each iterate.
    max_p_it = 1            # Max number of above steps to try.
    pre_it = 0              # Number of Poincare steps to do before anything else.
                            # These iterates will be discarded and the true cycling will.
                            # begin with the last iterate as y0.
    k = 2                   # Number of iterates to compute before doing the MPE approximation.
    kmax = min(6, len(y0))  # Max value of k to try.
    np = 100000
    np2 = 1000
    T_wasted = 0


    def poincare_step(z, t):
        " Integrate until the poincare section is reached again, in the right direction"
        # Integrate just a tiny step to make sure we don't stop at the
        # same point

        (y, tout, t_root, y_root, i_root, info_dict) = lsodar.odeintr(f, copy(z), linspace(t, t+t_small, np2), rtol=rtol, atol=atol, full_output = 1, root_func=event, root_term=array([0]))
        ynew = copy(y[-1,:])
        tnew = copy(tout[-1])
        print tnew

        # Main integration
        (y, tout, t_root, y_root, i_root, info_dict) = lsodar.odeintr(f, copy(ynew), linspace(tnew, tnew+max_t, np), rtol=rtol, atol=atol, full_output = 1, root_func=event, root_term=array([1]))
        ynew = copy(y[-1,:])
        if tout[-1] == tnew + max_t:
            raise RuntimeError, ("Never reached the Poincare section again. Oh noes!", ynew)
        tnew = copy(tout[-1])
        print tnew

        # Must check that the trajectory is going the right way
        if dot(f(ynew, tnew), n) < 0:
            print "Reached wrong side of the hyperplane, continuing..."
            # The solution cannot pass the hyperplane in the same direction two
            # times in a row, so just integrate one more time
            (y, tout, t_root, y_root, i_root, info_dict) = lsodar.odeintr(f, copy(ynew), linspace(tnew, tnew+max_t, np), rtol=rtol, atol=atol, full_output = 1, root_func=event, root_term=array([1]))
            ynew = copy(y[-1,:])
            tnew = copy(tout[-1])
            print tnew

        if len(t_root) == 0:
            raise RuntimeError, ("Never reached the Poincare section again. Oh noes!", ynew)
        else:
            print "Reached Poincare section again at", ynew

        return (ynew, tnew)


    def MPE(X):
        " Compute the minimum polynomial extrapolation (MPE) approximation to the sequence of vectors in X"
        uend = X[:,-1]-X[:,-2]
        U = X[:,1:-1]-X[:,0:-2]

        (c, res, rank, sing_vals) = lstsq(U,-uend)
        c = concatenate((c,array([1])))
        s = dot(X[:,0:-1],c)/sum(c)

        return s


    def main_it(i, p_it, its):
        # Take p_it Poincare steps
        print "Integrating... (", p_it, " times)"

        for j in xrange(0,p_it):
            try:
                (ynew, tnew) = poincare_step(x[i-1], T[i-1])
            except RuntimeError, e:
                print e[0]
                return e[1],0

            its = its+1
            # Check delta, adjust p_it if necessary
            delta = dist(ynew,x[i-1])
            if delta < d[i-1] and j < p_it-1:
                p_it = j+1
                print "Decreasing p_it to", p_it
                break


        if delta > d[i-1]:
            print "Delta:", delta
            print "Delta not decreasing with current p_it. (", p_it, ")"
            print "Iterating until delta has decreased..."
            j = 1
            while delta > d[i-1] and j < max_p_it:
                print "New iterate nr.", j

                try:
                    (ynew, tnew) = poincare_step(copy(ynew), copy(tnew))
                except RuntimeError, e:
                    print e[0]
                    return e[1],0

                its = its+1
                delta = dist(ynew,x[i-1])
                print "New delta nr.", j, " :", delta
                j = j+1

            if j == max_p_it:   # Delta not decreasing seems to be not due
                                # to multiple period of Poincare map.
                p_it = 1        # Reset to one step and hope for the best.
                print "Max. nr. of poincare steps per iterate reached. Resetting p_it to 1."
            else:
                p_it = p_it + j
                print "Found new p_it:", p_it

        d.append(delta)
        print "Delta:", delta

        # Store the iterate
        x.append(copy(ynew))
        xt.append(copy(transform(ynew,tnew)))
        T.append(copy(tnew))

        if delta < tol:  # We're done
            period = T[i]-T[i-1]
            return [True, delta, p_it, its, period]
        else:
            return [False, delta, p_it, its, 0]



    def end_processing(x, period, T, its, MPE_count):
        x = array(x)
        if filename == None:
            print "x: ", x
            print "Last y: ", x[-2]
##            print "Square sum: ", sum(x[-2]**2)
            print "Period: ", period
            print "Iterations: ", its
            print "MPE count: ", MPE_count
            print "Total time integrated: ", T
        else:
            f = open(filename, 'w')
            f.write( "x: " + array2string(x, separator=',') + '\n')
            f.write( "Last y: " + array2string(x[-2], separator=',') + '\n')
##            f.write( "Square sum: " + array2string(sum(x[-2]**2), separator=',') + '\n')
##            f.write( "Radius: " + array2string( sqrt(sum(x[-2]**2))*2/pi, separator=',') + ' pi/2\n')
            f.write( "Period: " + str(period) + '\n')
            f.write( "Iterations: " + str(its + pre_it) + '\n')
            f.write( "MPE count: " + str(MPE_count) + '\n')
            f.write( "Total time integrated: " + str(T) + '\n')
            f.write('\n\n')
            savetxt(f, x, delimiter=',')

        return x[-2], period





    # Do the pre-iteration
    xpre = [copy(y0)]
    print "Taking", pre_it, "Poincare steps before starting the main cycling."
    pre_T = 0
    for i in xrange(1,pre_it+1):
        (ynew, tnew) = poincare_step(xpre[i-1], i-1)
        tnew = tnew -(i-1)
        pre_T = pre_T + tnew
        xpre.append(ynew)
        if dist(xpre[i],xpre[i-1]) < tol:
            return end_processing(xpre, period = tnew, T = pre_T, its = i+1, MPE_count = 0)

    x = [xpre[-1]]


    # Cycle 0 - find k
    # Do temporary cycles with k = 2,3,...
    # and MPE-approximate fix-point s with each k
    # stop when ||f(s)-s|| < 1.e-1*||s||
    print "Finding k."
    k = 2
##    s = copy(x[0])
##    if dist(s, zeros(len(s))) == 0:
##        fs = ones(len(s))
##    else:
##        fs = copy(2*s)

    (done, delta, p_it, its, period) = main_it(1,p_it, its)
    if done: return end_processing(x, period, T[-1]+pre_T+T_wasted, its+pre_it, 0 )
    while k < kmax and not done:
        (done, delta, p_it, its, period) = main_it(k,p_it, its)
        if done: return end_processing(x, period, T[-1]+pre_T+T_wasted, its+pre_it, 0 )

        print "Computing MPE approximation..."
        X = array(xt).T     # We have generated x1, x2,...,x_{k}
        X = X[:,-(k+1):]    # Include the starting point x0
        s = MPE(X)

##        print "Taking Poincare step to evaluate how suitable k is."
        print "Doing a main iteration to evaluate how suitable k is."

        # Add the MPE approximation to the list
        xt.append(s)
        x.append(copy(inv_transform(s,T[-1])))
        T.append(T[-1])
        d.append(d[-1])#dist(x[-1],x[-2]))

        n_comp = len(x)
        # Find fs
        (done, delta, p_it, its, period) = main_it(k+2,p_it, its)
        if done: return end_processing(x, period, T[-1]+pre_T+T_wasted, its+pre_it, 0 )
        n_comp = len(x)-n_comp
        fs = x[-1]
        tnew = T[-1]

##        (fs, tnew) = poincare_step(copy(s), 0)

        print "Norm of diff.:", dist(fs,s)
        norm_s = dist(s, zeros(len(s)))
        print "Norm of s:", norm_s

        if dist(fs,s) < 1.e-1*norm_s:  # Consider k to be good enough
            MPE_count = MPE_count + 1
            break
        else:
            # Remove what we computed after the MPE approximation
            T_wasted += T[-1]-T[-2]
            del xt[-(1+n_comp):]
            del x[-(1+n_comp):]
            del T[-(1+n_comp):]
            del d[-(1+n_comp):]

        k = k+1
    else:
        k = kmax

    print "Found (estimate of) k:", k



    i = len(x)
    print i

    # Main cycling
    for j in xrange(1,maxit+1):
        print "Cycle", j

        # Compute the next k iterates (added to xt)
        for l in xrange(0,k): #xrange(j*k+MPE_count + 1,(j+1)*k+MPE_count + 1):
##            if not i == k + 2:
            (done, delta, p_it, its, period) = main_it(i,p_it, its)
            if done:
                return end_processing(x, period, T[-1]+pre_T+T_wasted, its + pre_it, MPE_count)
            i = i+1

        # Approximate the fix-point by minimum polynomial extrapolation
        # and restart from there

        print "Computing MPE approximation..."
        X = array(xt).T     # We have generated x1, x2,...,x_{k}
        X = X[:,-(k+1):]    # Include the starting point x0
        s = MPE(X)

        print "MPE:", s

        xt.append(s)
        x.append(copy(inv_transform(s,t0)))
        T.append(T[-1])

        d.append(delta)
        MPE_count = MPE_count + 1
        i = i + 1

    else:
        return end_processing(x, 0, T[-1]+pre_T+T_wasted, its + pre_it, MPE_count)

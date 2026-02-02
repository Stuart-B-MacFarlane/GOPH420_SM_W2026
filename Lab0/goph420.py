Python Code per Appendix D

bisect	p. 130

def bisect(func,xl,xu,es=1.e-7,maxit=30):
    """
    Uses the bisection method to estimate a root of func(x).
    The method is iterated until the relative error from
    one iteration to the next falls below the specified
    value or until the maximum number of iterations is
    reached first.
    Input:
        func = name of the function
        xl = lower guess
        xu = upper guess
        es = relative error specification  (default 1.e-7)
        maxit = maximum number of iterations allowed (default 30)
    Output:
        xm = root estimate
        fm = function value at the root estimate
        ea = actual relative error achieved
        i+1 = number of iterations required
        or
        error message if initial guesses do not bracket solution
    """
    if func(xl)*func(xu)>0:
        return 'initial estimates do not bracket solution'
    xmold = xl
    for i in range(maxit):
        xm = (xl+xu)/2
        ea = abs((xm-xmold)/xm)
        if ea < es:  break
        if func(xm)*func(xl)>0:
            xl = xm
        else:
            xu = xm
        xmold = xm
    return xm,func(xm),ea,i+1

import numpy as np

def f(m):
    g = 9.81
    cd = 0.25
    t = 4
    v = 36
    return np.sqrt(m*g/cd)*np.tanh(np.sqrt(g*cd/m)*t)-v

(m,fm,ea,iter) = bisect(f,50,200,maxit=20)
print('mass = {0:10.6f} kg'.format(m))
print('function value = {0:7.3g}'.format(fm))
print('relative error = {0:7.3g}'.format(ea))
print('iterations = {0:5d}'.format(iter))



bisect1		p. 129

def bisect1(func,xl,xu,maxit=20):
    """
    Uses the bisection method to estimate a root of func(x).
    The method is iterated maxit (default = 20) times.
    Input:
        func = name of the function
        xl = lower guess
        xu = upper guess
    Output:
        xm = root estimate
        or
        error message if initial guesses do not bracket solution
    """
    if func(xl)*func(xu)>0:
        return 'initial estimates do not bracket solution'
    for i in range(maxit):
        xm = (xl+xu)/2
        if func(xm)*func(xl)>0:
            xl = xm
        else:
            xu = xm
    return xm


brentsimp	p. 159

import numpy as np
eps = np.finfo(float).eps
def brentsimp(f,xl,xu):
    a = xl ; b = xu ; fa = f(a) ; fb = f(b)
    c = a ; fc = fa ; d = b - c ; e = d
    while True:
        if fb == 0: break
        if np.sign(fa) == np.sign(fb):  # rearrange points as req'd
            a = c ; fa = fc ; d = b - c ; e = d
        if abs(fa) < abs(fb):
            c = b ; b = a ; a = c
            fc = fb ; fb = fa ; fa = fc
        m = (a-b)/2  # termination test and possible exit
        tol = 2 * eps * max(abs(b),1)
        if abs(m) < tol or fb == 0: break
        # choose open methods or bisection
        if abs(e) >= tol and abs(fc) > abs(fb):
            s = fb/fc
            if a == c: 
                # secant method here
                p = 2*m*s
                q = 1 - s
            else:
                # inverse quadratic interpolation here
                q = fc/fa ; r = fb/fa
                p = s * (2*m*q*(q-r)-(b-c)*(r-1))
                q = (q-1)*(r-1)*(s-1)
            if p > 0: 
                q = -q
            else:
                p = -p
            if 2*p < 3*m*q - abs(tol*q) and p < abs(0.5*e*q):
                e = d ; d = p/q
            else:
                d = m ; e = m
        else:
            # bisection here
            d = m; e = m
        c = b ; fc = fb
        if abs(d) > tol:
            b = b + d
        else:
            b = b - np.sign(b-a)*tol
        fb = f(b)
    return b

xsoln = brentsimp(lambda x: x**10-1,0.1,5)
print('Solution = ',xsoln)


cspline		p. 433

import numpy as np

def cspline(x,y,xx):
    """
    Cubic Spline Interpolation with Natural End Conditions
    input:
        x = array of independent variable values
        y = array of dependent variable values
        xx = input value for interpolation
    output:
        yy = interpolated value of y
    """
    n = len(x)
    if len(y) != n: return 'input arrays must be the same length'
    if xx < x[0] or xx > x[n-1]: return 'input value out of range of table'
    h = np.zeros((n-1))
    for i in range(n-1):
        h[i] = x[i+1] - x[i]
    df = np.zeros((n-1))
    for i in range(n-1):
        df[i] = (y[i+1]-y[i])/h[i]
    A = np.zeros((n,n))
    A[0,0] = 1 ; A[n-1,n-1] = 1
    for i in range(1,n-1):
        A[i,i-1] = h[i-1]
        A[i,i] = 2*(h[i-1] + h[i])
        A[i,i+1] = h[i]
    const = np.zeros((n))
    for i in range(1,n-1):
        const[i] = 3*(df[i]-df[i-1])
    c = np.linalg.solve(A,const)
    b = np.zeros((n-1))
    for i in range(n-1):
        b[i] = (y[i+1]-y[i])/h[i]-h[i]/3*(2*c[i]+c[i+1])
    d = np.zeros((n-1))
    for i in range(n-1):
        d[i] = (c[i+1]-c[i])/3/h[i]
    for i in range(n):
        if xx == x[i]:  # check for an exact match
            return y[i]
        elif x[i] > xx:  # check for upper interval
            i2 = i-1
            break
    yy = y[i2] + b[i2]*(xx-x[i2]) + c[i2]*(xx-x[i2])**2 + d[i2]*(xx-x[i2])**3
    return yy

T = np.array([-40.,   0.,  20., 50., 100., 150., 200., 250.,
              300., 400., 500.])
rho = np.array([1.52,  1.29,  1.20,  1.09,  0.946, 0.835,
                0.746, 0.675, 0.616, 0.525, 0.457])
    
Tx = 350.
rhox = cspline(T,rho,Tx)
print(rhox)

import pylab

Tplot = np.linspace(-40.,500.,100)
k = len(Tplot)
rhoplot = np.zeros((k))
for i in range(k):
    rhoplot[i] = cspline(T,rho,Tplot[i])
pylab.scatter(T,rho,c='k',marker='s')
pylab.plot(Tplot,rhoplot,c='k',ls = ':')
pylab.grid()


csplinesm	p. 442

import numpy as np

def csplinesm(x,y,xx,lam,sdest):
    mu = 2*(1-lam)/lam/3
    m = len(x)
    if len(y) != m: return 'x and y arrays must be the same length'
    n = m - 1
    h = np.zeros((n))  # compute x intervals
    for i in range(n):
        h[i] = x[i+1] - x[i]
    # compute elements for R and Q' matrices
    p = np.zeros((n))
    for i in range(1,n):
        p[i] = 2*(-h[i]+h[i-1])
    r = np.zeros((n))
    for i in range((n)):
        r[i] = 3/h[i]
    f = np.zeros((n))
    for i in range(1,n):
        f[i] = - 3*(1/h[i-1]+1/h[i])
    # compose R matrix
    R = np.zeros((n-1,n-1))
    for i in range(n-2):
        R[i,i] = p[i]
        R[i+1,i] = h[i]
        R[i,i+1] = h[i]
    R[n-2,n-2] = p[n-2]
    Qp = np.zeros((n-1,n+1))
    for i in range(n-2):
        Qp[i,i] = r[i]
        Qp[i,i+1] = f[i+1]
        Qp[i,i+2] = r[i+1]
    Qp[n-2,n-1] = f[n-1]
    Qp[n-2,n-2] = r[n-2]
    Qp[n-2,n] = r[n-1]
    #  Q from Q'
    Q = np.transpose(Qp)
    # diagonal matrix of sigma estimates
    SigMat = np.zeros((n+1,n+1))
    for i in range(n+1):
        SigMat[i,i] = sdest
    # set up linear equations to solve for c
    Qt = np.dot(Qp,SigMat)
    Qcoef = np.dot(Qt,Q)
    Qcoef = Qcoef*mu + R
    const = np.dot(Qp,y)
    # solve for c
    c = np.linalg.solve(Qcoef,const)
    # solve for d
    SigQ = np.dot(SigMat,Q)
    yc = np.dot(SigQ,c)
    yc = yc*mu
    a = y - yc
    # solve for d and b
    d = np.zeros((n))
    b = np.zeros((n))
    cx = np.zeros((n))
    for i in range(1,n-1):
        cx[i] = c[i-1]
    for i in range(0,n-1):
        d[i] = (cx[i+1]-cx[i])/3/h[i]
        b[i] = (a[i+1]-a[i]/h[i])
    # compute interpolation
    for i in range(n):   # calculate interpolation
        if xx == x[i]:  # check for an exact match
            return y[i],a,b,c,d 
        elif x[i] > xx:  # check for upper interval
            i2 = i-1
            break
    yy = a[i2] + b[i2]*(xx-x[i2]) + c[i2]*(xx-x[i2])**2 + d[i2]*(xx-x[i2])**3
    # return interpolated value and spline coefficient arrays
    return yy,a,b,c,d 

conc = np.loadtxt(fname='ConcentrationData.txt')
n = len(conc)
samp = np.zeros((n))
for i in range(n):
    samp[i] = float(i)
xx = 100.5
lam = 0.001
sconc = np.std(conc)

concint,a,b,c,d = csplinesm(samp,conc,xx,lam,sconc)
print(concint)

import pylab
pylab.scatter(samp,conc,c='k',marker='s')
pylab.plot(samp,a,c='k')
pylab.grid()
pylab.xlabel('sample number')
pylab.ylabel('concentration')


diffex	p. 108

import numpy as np
def diffex(func,dfunc,x,n):
    dftrue = dfunc(x)
    h = 1
    H = np.zeros(n)
    D = np.zeros(n)
    E = np.zeros(n)
    H[0] = h
    D[0] = (func(x+h)-func(x-h))/2/h
    E[0] = abs(dftrue-D[0])
    for i in range(1,n):
        h = h/10
        H[i] = h
        D[i] = (func(x+h)-func(x-h))/2/h
        E[i] = abs(dftrue-D[i])
    return H,D,E

ff = lambda x: -0.1*x**4 - 0.15*x**3 - 0.5*x**2 - 0.25*x + 1.2
df = lambda x: -0.4*x**3 - 0.45*x**2 - x - 0.25

H,D,E = diffex(ff,df,0.5,11)
print( '   step size    finite difference      true error')
for i in range(11):
    print('{0:14.10f}  {1:16.14f}   {2:16.13f}'.format(H[i],D[i],E[i]))


eulode		p. 538

import numpy as np

def eulode(dydt,tspan,y0,h,*args):
    """
    solve initial-value single ODEs with the Euler method
    input:
        dydt = function name that evaluates the derivative
        tspan = array of [ti,tf] where
            ti and tf are the initial and final values
            of the independent variable
        y0 = initial value of the dependent variable
        h = step size
        *args = additional argument to be passed to dydt
    output:
        t = an array of independent variable values
        y - an array of dependent variable values
    """
    ti = tspan[0] ; tf = tspan[1]
    if not(tf>ti+h): return 'upper limit must be greater than lower limit'
    t = []
    t.append(ti)  # start the t array with ti
    nsteps = int((tf-ti)/h)
    for i in range(nsteps):  # add the rest of the t values
        t.append(ti+(i+1)*h)
    n = len(t)    
    if t[n-1] < tf:  # check if t array is short of tf
        t.append(tf)
        n = n+1
    y = np.zeros((n)) ; y[0] = y0  # initialize y array
    for i in range(n-1):
        y[i+1] = y[i] + dydt(t[i],y[i],*args)*(t[i+1]-t[i])  # Euler step
    return t,y

dydt = lambda t,y,a,b,c: b*np.exp(a*t) - c*y
tspan = np.array([0.,4.])
y0 = 2.
h = 0.9
a = 0.8
b = 4.
c = 0.5
t,y = eulode(dydt,tspan,y0,h,a,b,c)
n = len(t)
for i in range(n):
    print('{0:4.1f}  {1:7.4f}'.format(t[i],y[i]))
    
    
fixpt	p. 147

def fixpt(g,x0,Ea=1.e-7,maxit=30):
    """
    This function solves x=g(x) using fixed-point iteration.
    The method is repeated until either the relative error
    falls below Ea (default 1.e-7) or reaches maxit (default 30).
    Input:
        g = name of the function for g(x)
        x0 = initial guess for x
        Ea = relative error threshold
        maxit = maximum number of iterations
    Output:
        x1 = solution estimate
        ea = relative error
        i+1 = number of iterations
    """    
    for i in range(maxit):
        x1 = g(x0)
        ea = abs((x1-x0)/x1)
        if ea < Ea:  break
        x0 = x1
    return x1,ea,i+1

import numpy as np
rho = 1.23  # kg/m3
mu = 1.79e-5  # Pa*s
D = 0.005  # m
V = 40  # m/s
eps = 1.5e-6  # m
Re = rho*V*D/mu  # Reynolds number

def g(f):
    fr = -4*np.log10(eps/3.7/D+1.26/Re/np.sqrt(f))
    return 1/fr**2

f0 = 0.001
fsoln,ea,n = fixpt(g,f0)
print('Fanning friction factor = ',fsoln)
print('Relative error = ',ea)
print('Number of iterations = ',n)


funcavg		p. 71

import numpy as np
def funcavg(a,b,n):
    """
    computes the average value of a function
    over a range
    requires import of numpy module
    input:
        a = lower bound of range
        b = upper bound of range
        n = number of intervals
    output:
        favg = average value of function
    """
    x = np.linspace(a,b,n)
    y = func(x)
    favg = np.average(y)
    return favg

def func(t):
    g = 9.81
    m = 68.1
    cd = 0.25
    f = np.sqrt(m*g/cd)*np.tanh(np.sqrt(cd*g/m)*t)
    return f



gaussnaive	p. 227

import numpy as np

def gaussnaive(A,b):
    """
    gaussnaive: naive Gauss elimination
    input:
    A = coefficient matrix
    b = constant vector
    output:
    x = solution vector
    """
    (n,m) = A.shape
    #n = nm[0]
    #m = nm[1]
    if n != m:
        return 'Coefficient matrix A must be square'
    nb = n+1
    # build augmented matrix
    Aug = np.hstack((A,b))
       # forward elimination
    for k in range(n-1):
        for i in range(k+1,n):
            factor = Aug[i,k]/Aug[k,k]
            Aug[i,k:nb]=Aug[i,k:nb]-factor*Aug[k,k:nb]
    # back substitution
    x = np.zeros([n,1])  # create empty x array
    x = np.matrix(x)  # convert to matrix type
    x[n-1]=Aug[n-1,nb-1]/Aug[n-1,n-1]
    for i in range(n-2,-1,-1):
        x[i]=(Aug[i,nb-1]-Aug[i,i+1:n]*x[i+1:n,0])/Aug[i,i]
    return x

A = np.matrix('0.55,0.25,0.25 ; 0.30,0.45,0.20 ; 0.15,0.30,0.55')
b = np.matrix('4800;5800,5700')
V = gaussnaive(A,b)
print('Volume from Pit 1:',V[0])
print('Volume from Pit 2:',V[1])
print('Volume from Pit 3:',V[2])


gausspivot	p. 233

import numpy as np

def gausspivot(A,b):
    """
    gausspivot: Gauss elimination with partial pivoting
    input:
    A = coefficient matrix
    b = constant vector
    output:
    x = solution vector
    """
    (n,m) = A.shape
    if n != m:
        return 'Coefficient matrix A must be square'
    nb = n+1
    # build augmented matrix
    Aug = np.hstack((A,b))
       # forward elimination
    for k in range(n-1):

        # partial pivoting
        imax = maxrow(Aug[k:n,k])
        ipr = imax + k
        if ipr != k:  # no row swap if pivot is max
            for j in range(k,nb):  # swap rows k and ipr
                temp = Aug[k,j]
                Aug[k,j] = Aug[ipr,j]
                Aug[ipr,j] = temp

        for i in range(k+1,n):
            factor = Aug[i,k]/Aug[k,k]
            Aug[i,k:nb]=Aug[i,k:nb]-factor*Aug[k,k:nb]
    # back substitution
    x = np.zeros([n,1])  # create empty x array
    x = np.matrix(x)  # convert to matrix type
    x[n-1]=Aug[n-1,nb-1]/Aug[n-1,n-1]
    for i in range(n-2,-1,-1):
        x[i]=(Aug[i,nb-1]-Aug[i,i+1:n]*x[i+1:n,0])/Aug[i,i]
    return x

def maxrow(avec):
    # function to determine the row index of the
    # maximum value in a vector
    maxrowind = 0
    n = len(avec)
    amax = abs(avec[0])
    for i in range(1,n):
        if abs(avec[i]) > amax:
            amax = avec[i]
            maxrowind = i
    return maxrowind
   
    
A = np.matrix('0 -0.1 -0.2 ; 0.1 7 -0.3 ; 0.3 -0.2 10')
b = np.matrix(' 7.85 ; -19.3 ; 71.4')
x = gausspivot(A,b)
print(x)
xtest = np.linalg.inv(A)*b
print(xtest)


GaussSeidel	p. 269

import numpy as np

def GaussSeidel(A,b,es=1.e-7,maxit=50):
    """
    Implements the Gauss-Seidel method
    to solve a set of linear algebraic equations
    without relaxation
    Input:
    A = coefficient matris
    b = constant vector
    es = stopping criterion (default = 1.e-7)
    maxit = maximum number of iterations (default=50)
    Output:
    x = solution vector
    """
    n,m = np.shape(A)
    if n != m :
        return 'Coefficient matrix must be square'
    C = np.zeros((n,n))
    x = np.zeros((n,1))
    for i in range(n):  # set up C matrix with zeros on the diagonal
        for j in range(n):
            if i != j:
                C[i,j] = A[i,j]
    d = np.zeros((n,1))
    for i in range(n):  # divide C elements by A pivots
        C[i,0:n] = C[i,0:n]/A[i,i]
        d[i] = b[i]/A[i,i]
    ea = np.zeros((n,1))
    xold = np.zeros((n,1))
    for it in range(maxit):  # Gauss-Seidel method
        for i in range(n):
            xold[i] = x[i]  # save the x's for convergence test
        for i in range(n):
            x[i] = d[i] - C[i,:].dot(x)  # update the x's 1-by-1
            if x[i] != 0:
                ea[i] = abs((x[i]-xold[i])/x[i])  # compute change error
        if np.max(ea) < es:  # exit for loop if stopping criterion met
            break
    if it == maxit:  # check for maximum iteration exit
        return 'maximum iterations reached'
    else:
        return x
    
A = np.matrix(' 3. -0.1 -0.2 ; 0.1 7. -0.3 ; 0.3 -0.2 10. ')
b = np.matrix(' 7.85 ; -19.3 ; 71.4 ')
x = GaussSeidel(A,b)
print('Solution is\n',x)

x2 = np.linalg.solve(A,b)
print('2nd solution is\n',x2)


goldmin	p. 183

import numpy as np

def goldmin(f,xl,xu,Ea=1.e-7,maxit=30):
    """
    use the golden-section search to find the minimum of f(x)
    input:
        f = name of the function
        xl = lower initial guess
        xu = upper initial guess
        Ea = absolute relative error criterion (default = 1.e-7)
        maxit = maximum number of iterations (default = 30)
    output:
        xopt = location of the minimum
        f(xopt) = function value at the minimum
        ea = absolute relative error achieved
        i+1 = number of iterations required
    """
    phi = (1+np.sqrt(5))/2
    d = (phi - 1)*(xu-xl)
    x1 = xl + d ; f1 = f(x1)
    x2 = xu - d ; f2 = f(x2)
    for i in range(maxit):
        xint = xu - xl
        if f1 < f2:
            xopt = x1
            xl = x2
            x2 = x1
            f2 = f1
            x1 = xl + (phi-1)*(xu-xl)
            f1 = f(x1)
        else:
            xopt = x2
            xu = x1
            x1 = x2
            f1 = f2
            x2 = xu - (phi-1)*(xu-xl)
            f2 = f(x2)
        if xopt != 0:
            ea = (2-phi)*abs(xint/xopt)
            if ea <= Ea: break
    return xopt,f(xopt),ea,i+1

g = 9.81  # m/s2
v0 = 55  # m/s
m = 80  # kg
c = 15  # kg/s
z0 = 100  # m

def f(t):
    return -(z0+m/c*(v0+m*g/c)*(1-np.exp(-t/(m/c)))-m*g/c*t)

tl = 0
tu = 8
tmin,fmin,ea,n = goldmin(f,tl,tu,Ea=1.e-5)

print('Time at maximum altitude = {0:5.2f} s'.format(tmin))
print('Function value = {0:6.2g} '.format(fmin))
print('Relative error = {0:7.2e} '.format(ea))
print('Iterations required = {0:4.0f} '.format(n))

zmax = z0 + m/c*(v0+m*g/c)*(1-np.exp(-tmin/(m/c)))-m*g/c*tmin
print('Maximum altitude = {0:6.2f} m'.format(zmax))


incsearch	p. 123

import numpy as np

def incsearch(func,xmin,xmax,ns=50):
    """
    incsearch: incremental search locator
        incsearch(func,xmin,xmax,ns)
        finds brackets of x that contain sign changes in
        a function of x on an interval
    input: 
        func = name of the function
        xmin, xmax = endpoints of the interval
        ns = number of subintervals, default value = 50
    output:  a tuple containing
        nb = number of bracket pairs found
        xb = list of bracket pair values
        or returns "no brackets found"
    """
    x = np.linspace(xmin,xmax,ns) # create array of x values
    f = []  # build array of corresponding function values
    for k in range(ns-1):  
        f.append(func(x[k]))
    nb = 0
    xb = []
    for k in range(ns-2):  # check adjacent pairs of function values
        if func(x[k])*func(x[k+1])<0:  # for sign change
            nb = nb + 1  # increment the bracket counter
            xb.append((x[k],x[k+1]))  # save the bracketing pair
    if nb==0:
        return 'no brackets found'
    else:
        return nb,xb


itermeth	p. 89

import math
def itermeth(x,es=1e-4,maxit=50):
    """
    Maclaurin series expansion of the exponential function
    requires math module
    input:
        x = value at which the series is evaluated
        es = stopping criterion (default = 1e-4)
        maxit = maximum number of iterations (default=10)
    output:
        fx = estimated function value
        ea = approximate relative error (%)
        iter = number of iterations
    """
    # initialization
    iter = 1 ; sol = 1 ; ea = 100
    # iterative calculation
    while True:
        solold = sol
        sol = sol + x**iter / math.factorial(iter)
        iter = iter + 1
        if sol != 0: ea = abs((sol-solold)/sol)*100
        if ea < es or iter == maxit: break
    fx = sol
    return fx,ea,iter


Lagrange	p. 407

import numpy as np

def Lagrange(x,y,xx):
    """
    Lagrange interpolating polynomial
    Uses an (n-1)th-order Lagrange interpolating polynomial
    based on n data pairs to return a value of the
    dependent variable, yint, at a given value of the
    independent variable, xx.
    Input:
        x = array of independent variable values
        y = array of dependent variable values
        xx = value of independent variable at which
             the interpolation is calculated
    Output:
        yint = interpolated value of the dependent variable
    """
    n = len(x)
    if len(y) != n:
        return 'x and y must be of same length'
    s = 0
    for i in range(n):
        product = y[i]
        for j in range(n):
            if i != j:
                product = product * (xx - x[j])/(x[i]-x[j])
        s = s + product
    yint = s
    return yint

T = np.array([-40., 0., 20., 50.])
rho = np.array([1.52, 1.29, 1.2, 1.09])
rhoint = Lagrange(T,rho,15.)
print(rhoint)


loess	p. 446

import numpy as np

def sortdist(dist,n):
    """
    function to sort the dist array
    into the sdist array
    using an insertion sort
    the original indices of the distances
    are retained in ind
    """
    pos = np.zeros((n),dtype=np.bool)
    ind = np.zeros((n),dtype=np.int16)
    sdist = np.zeros((n))
    for i in range(n):
        pos[i] = False
    maxdist = np.max(dist)
    for i in range(n):
        mindist = maxdist
        for j in range(n):
            if dist[j] <= mindist and not pos[j]:
                mindist = dist[j]
                minloc = j
        sdist[i] = mindist
        ind[i] = minloc
        pos[minloc] = True
    return sdist,ind

def loess(x0,x,y,alpha):
    """
    loess smoothing applied to the series {x,y}
    Input:
        x0 = independent variable value for interpolation
        x = independent variable series
        y = dependent variable series
        alpha = smoothing parameter between 0 and 1
    Output:
        ys = smoothed values of y
    """
    n = len(x)
    if n != len(y): return 'x and y series must be of same length'
    k = int(alpha*n)  # how many data for loess?
    # compute distance array
    dist = np.zeros((n))
    for i in range(n):
        dist[i] = abs(x0-x[i])
    # sort the distance array
    # with an accompanying index array
    sdist,ind = sortdist(dist,n)
    # extract nearest x and y data to x0
    Nx = np.zeros((k))
    Ny = np.zeros((k))
    for i in range(k):
        Nx[i] = x[ind[i]]
        Ny[i] = y[ind[i]]
    # set value for weight calculation
    # based on farthest distance
    delx0 = sdist[k]
    # zero out the rest of the sdist array
    for i in range(k,n):
        sdist[i] = 0
    # compute the weights and the diagonal
    # weight matrix
    z = np.zeros((k))
    w = np.zeros((k))
    for i in range(k):
        z[i] = sdist[i]/delx0
        w[i] = (1-z[i]**3)**3
    W = np.zeros((k,k))
    for i in range((k)):
        W[i,i] = w[i]
    # build the X matrix
    X = np.zeros((k,3))
    for i in range(k):
        X[i,0] = Nx[i]**2
        X[i,1] = Nx[i]
        X[i,2] = 1
    # formulate the coefficient matrix
    # and constant vector
    # for the normal equations
    Xt = np.transpose(X)
    XtW = np.dot(Xt,W)
    Xcoef = np.dot(XtW,X)
    const = np.dot(XtW,Ny)
    # solve the normal equations
    b = np.linalg.solve(Xcoef,const)
    # use the coefficients for the interpolation
    yhat0 = b[0]*x0**2 + b[1]*x0 + b[2]
    return yhat0

conc = np.loadtxt(fname='ConcentrationData.txt')
n = len(conc)
samp = np.zeros((n))
for i in range(n):
    samp[i] = float(i)
xx = 100.5
alpha = 0.35

concint = loess(xx,samp,conc,alpha)
print(concint)

concsm = np.zeros((n))
for i in range(n):
    concsm[i] = loess(samp[i],samp,conc,alpha)

import pylab
pylab.scatter(samp,conc,c='k',marker='s')
pylab.plot(samp,concsm,c='k')
pylab.grid()
pylab.xlabel('sample number')
pylab.ylabel('concentration')


Newtint	p. 405

import numpy as np

def Newtint(x,y,xx):
    """
    Newtint: Newton interpolating polynomial
    Uses an (n-1)th-order Newton interpolating polynomial
    based on n data pairs to return a value of the
    dependent variable, yint, at a given value of the
    independent variable, xx.
    Input:
        x = array of independent variable values
        y = array of dependent variable values
        xx = value of independent variable at which
             the interpolation is calculated
    Output:
        yint = interpolated value of the dependent variable
    """
    # compute the finite divided differences in the 
    # form of a difference table
    n = len(x)
    if len(y) != n:
        return 'x and y must be of same length'
    b = np.zeros((n,n))
    # assign the dependent variables to the first column of b
    b[:,0] = np.transpose(y)
    for j in range(1,n):
        for i in range(n-j):
            b[i,j] = (b[i+1,j-1]-b[i,j-1])/(x[i+j]-x[i])
    # use the finite divided differences to interpolate
    xt = 1
    yint = b[0,0]
    for j in range(n-1):
        xt = xt * (xx - x[j])
        yint = yint + b[0,j+1]*xt
    return yint

x = np.array([1., 4., 6., 5.])
y = np.log(x)
yi = Newtint(x,y,2.)
print(yi)


newtmult	276

import numpy as np

def newtmult(fandJ,x0,es=1.e-7,maxit=20):
    """
    Newton-Raphson solution of sets of nonlinear algebraic equations
    Input:
    fandJ = function name that supplies f and Jacobian values
    x0 = initial guesses for x
    es = convergence tolerance (default = 1.3-7)
    maxit = iteration limit (default = 20)
    Output:
    x = solution
    f = function values at the solution
    ea = relative error
    iter = number of iterations taken
    """
    n,m = np.shape(x0)  # get the number of equations in n
    x = np.zeros((n,m))
    for i in range(n):  # initialize x
        x[i] = x0[i]
    for i in range(maxit):
        f,J = fandJ(x)   # get the function values and the Jacobian
        dx = np.linalg.inv(J).dot(f)  # Newton-Raphson iteration
        x = x - dx
        ers = dx/x
        ea = max(abs(ers))
        if ea < es: break  # check for convergence
    if i == maxit:
        return 'iteration limit reached'
    else:
        return x,f,ea,i+1  # here if solution successful


newtraph	p. 153

def newtraph(f,fp,x0,Ea=1.e-7,maxit=30):
    """
    This function solves f(x)=0 using the Newton-Raphson method.
    The method is repeated until either the relative error
    falls below Ea (default 1.e-7) or reaches maxit (default 30).
    Input:
        f = name of the function for f(x)
        fp = name of the function for f'(x)
        x0 = initial guess for x
        Ea = relative error threshold
        maxit = maximum number of iterations
    Output:
        x1 = solution estimate
        f(x1) = equation error at solution estimate
        ea = relative error
        i+1 = number of iterations
    """    
    for i in range(maxit):
        x1 = x0 - f(x0)/fp(x0)
        ea = abs((x1-x0)/x1)
        if ea < Ea:  break
        x0 = x1
    return x1,f(x1),ea,i+1

def f(x):
    return x**2-9

def fp(x):
    return 2*x

x0 = 5
(xsoln,fxsoln,ea,n) = newtraph(f,fp,x0,Ea=1.e-5)
print('Solution = {0:8.5g}'.format(xsoln))
print('Function value at solution = {0:8.5e}'.format(fxsoln))
print('Relative error = {0:8.3e}'.format(ea))
print('Number of iterations = {0:5d}'.format(n))


odesimp	p. 76

def odesimp(dydt,dt,ti,tf,yi):
   t = ti ; y = yi ; h = dt
   while True:
       if t + dt > tf: h = tf - t
       y = y + dydt(y)*h
       t = t + h
       if t >= tf: break
   yend = y
   return yend


odesimp2	p. 76

def odesimp2(dydt,dt,ti,tf,yi,*args):
   t = ti ; y = yi ; h = dt
   while True:
       if t + dt > tf: h = tf - t
       y = y + dydt(y,*args)*h
       t = t + h
       if t >= tf: break
   yend = y
   return yend

def dvdt(v,cd,m):
    g=9.81
    return g - (cd/m)*v*abs(v)

print('{0:7.3f}'.format(odesimp2(dvdt,0.5,0,12,0,0.25,68.1)))


quadadapt	p. 498

def quadadapt(func,a,b,tol=1.e-8):
    """
    Evaluates the definite integral of f(x) from a to b
    """
    c = (a+b)/2
    fa = func(a) ; fb = func(b) ; fc = func(c)
    q = quadstep(func,a,b,tol,fa,fc,fb)
    return q
def quadstep(func,a,b,tol,fa,fc,fb):
    h = b - a ; c = (a+b)/2
    fd = func((a+c)/2) ; fe = func((c+b)/2)
    q1 = h/6 * (fa + 4*fc + fb)
    q2 = h/12 * (fa + 4*fd + 2*fc + 4*fe + fb)
    if abs(q1-q2) < tol:
        q = q2 + (q2-q1)/15
    else:
        qa = quadstep(func,a,c,tol,fa,fd,fc)
        qb = quadstep(func,c,b,tol,fc,fe,fb)
        q = qa + qb
    return q

def f(x):
    return 0.2 + 25.*x - 200.*x**2 + 675.*x**3 - 900*x**4 + 400*x**5

fint = quadadapt(f,0.,0.8)
print(fint)


quadroots	p. 67

def quadroots(a,b,c):
    """
    quadroots:  roots of the quadratic equation
        quadroots(a,b,c): real and complex roots
                          of a quadratic polynomial
    input:
        a = second-order coefficient
        b = first-order coefficient
        c = zero-order coefficient
    outpuit:
        r1: real part of the first root
        i1: imaginary part of first root
        r2: real part of the second root
        i2: imaginary part of second root
    """
    import math
    if a == 0:
        # special cases
        if b != 0:
            # single root
            r1 = -c/b
            return r1
        else:
            # trivial solution
            print('Trivial solution. Try again.')
    else:
        # quadratic formula
        d = b**2 - 4*a*c
        if d >= 0:
            # real roots
            r1 = (-b+math.sqrt(d))/(2*a)
            r2 = (-b-math.sqrt(d))/(2*a)
            i1 = 0
            i2 = 0
        else:
            # complex roots
            r1 = -b/(2*a)
            i1 = math.sqrt(abs(d))/(2*a)
            r2 = r1
            i2 = -i1
        return r1,i1,r2,i2
        


rk4sys		p. 551

import numpy as np

def rk4sys(dydt,tspan,y0,h,*args):
    """
    fourth-order Runge-Kutta method
    for solving a system of ODEs
    input:
        dydt = function name that evaluates the derivatives
        tspan = array of [ti,tf] where
            ti and tf are the initial and final values
            of the independent variable
        y0 = initial value of the dependent variable
        h = step size
        *args = additional argument to be passed to dydt
    output:
        t = array of independent variable values
        y = array of dependent variable values
    """
    ti = tspan[0] ; tf = tspan[1]
    if not(tf>ti+h): return 'upper limit must be greater than lower limit'
    t = []
    t.append(ti)  # start the t array with ti
    nsteps = int((tf-ti)/h)
    for i in range(nsteps):  # add the rest of the t values
        t.append((i+1)*h)
    n = len(t)    
    if t[n-1] < tf:  # check if t array is short of tf
        t.append(tf)
        n = n+1
    neq = len(y0)
    y = np.zeros((n,neq))  # set up 2-D array for dependent variables
    for j in range(neq):
        y[0,j] = y0[j]  #  set first elememts to initial conditions
    for i in range(n-1):  # 4th order RK
        hh = t[i+1] - t[i]
        k1 = dydt(t[i],y[i,:],*args)
        ymid = y[i,:] + k1*hh/2.
        k2 = dydt(t[i]+hh/2.,ymid,*args)
        ymid = y[i,:] + k2*hh/2.
        k3 = dydt(t[i]+hh/2.,ymid,*args)
        yend = y[i,:] + k3*hh
        k4 = dydt(t[i]+hh,yend,*args)
        phi = (k1 + 2.*(k2+k3) + k4)/6.
        y[i+1,:] = y[i,:] + phi*hh
    return t,y

def dydtsys(t,y):
    n = len(y)
    dy = np.zeros((n))
    dy[0] = -2.*y[0]**2 +2.*y[0] + y[1] - 1.
    dy[1] = -y[0] -3*y[1]**2 +2.*y[1] + 2.
    return dy

ti = 0. ; tf = 2.
tspan = np.array([ti,tf])
h = 0.01
y0 = np.array([2.,0.])
t,y = rk4sys(dydtsys,tspan,y0,h)             

import pylab
pylab.plot(t,y[:,0],c='k',label='y1')
pylab.plot(t,y[:,1],c='k',ls='--',label='y2')
pylab.grid()
pylab.xlabel('t')
pylab.ylabel('y')
pylab.legend()
        


rkf45		p. 564

import numpy as np

def rkf45(dydx,xspan,y0,es=1.e-6,maxit=50,hmin=1.e-15,*args):
    """
    Runge-Kutta-Fehlberg 4/5 algorithm for the 
    solution of one or more ODEs
    input:
        dydx = function name that evaluates the derivatives
        xspan = an array of independent variable values
            where the solution will be returned
        y0 = an array initial values of the dependent variable(s)
        es = local relative error tolerance (default = 1.e-10)
        hmin = minimum step size
        *args = additional arguments to be passed to dydt
    output:
        t = array of independent variable values
        y = array of dependent variable values
    """
    # compute all coefficients
    a2 = 0.25 ; a3 = 0.375 ; a4 = 12./13. ; a5 = 1.; a6 = 0.5
    b21 = 0.25 ; 
    b31 = 3./32. ; b32 = 9./32.
    b41 = 1932./2197. ; b42 = -7200./2197. ; b43 = 7296./2197.
    b51 = 439./216. ; b52 = -8. ; b53 = 3680./513. ; b54 = -845./4104.
    b61 = -8./27. ; b62 = 2. ; b63 = -3544./2565. ; b64 = 1859./4104.
    b65 = -11./40.
    c1 = 25./216. ; c3 = 1408./2565. ; c4 = 2197./4101. ; c5 = -0.2
    d1 = 16./135. ; d3 = 6656./12825. ; d4 = 28561./56430. 
    d5 = -9./50. ; d6 = 2./55.
        
    n = len(xspan)  # here if tspan contains step times
    x = xspan
    neq = np.size(y0)
    y = np.zeros((n,neq))  # set up 2-D array for dependent variables
    if neq > 1:
        y[0,:] = y0[:]
    else:
        y[0] = y0
        cd =-1
    hnew = x[1]-x[0]
    for i in range(n-1): # integrate steps given in tspan
        h = hnew
        xk = x[i] ; yk = y[i,:]
        while True:
            if xk >= x[i+1]: break
            for k in range(maxit):
                k1 = h*dydx(xk,yk,*args)  # compute k factors
                k2 = h*dydx(xk+a2*h,yk+b21*k1,*args)
                k3 = h*dydx(xk+a3*h,yk+b31*k1+b32*k2,*args)
                k4 = h*dydx(xk+a4*h,yk+b41*k1+b42*k2+b43*k3,*args)
                k5 = h*dydx(xk+a5*h,yk+b51*k1+b52*k2+b53*k3+b54*k4,*args)
                k6 = h*dydx(xk+a6*h,yk+b61*k1+b62*k2+b63*k3+b64*k4+b65*k5,*args)
                y4 = yk + c1*k1 + c3*k3 + c4*k4 + c5*k5  # 4th-order estimate
                y5 = yk + d1*k1 + d3*k3 + d4*k4 + d5*k6 + d6*k6  # 5th order
                yerr = abs((y5 - y4)/y4)  # error, perhaps array of errors
                yerrm = max(yerr)  # max error
                s = 1.
                if yerrm != 0:
                    s = 0.84 * (es/yerrm)**0.25  # s factor
                    if s < 0.25:  # clamp s factor, if necessary
                        s = 0.25
                    elif s > 4.0:
                        s = 4.0
                hnew = s*h
                if xk + hnew > x[i+1]:
                    hnew = x[i+1] - xk
                elif x[i] + 1.5*hnew > x[i+1]:
                    hnew = hnew/2.
                if hnew < hmin: return x,y,cd
                if yerrm < es: break
                h = hnew
            if k == maxit-1: 
                cd = 0
                return x,y,cd
            xk = xk + hnew ;  yk = y4 ; h = hnew
        y[i+1,:] = y4
    cd=1
    return x,y,cd

dydt = lambda t,y: 4.*np.exp(0.8*t)-0.5*y

tspan = np.linspace(0.,4.,11)
y0 = 2.
t,y,cd = rkf45(dydt,tspan,y0,es=1.e-8)
print(cd)
n = len(t)
print(t,y)

from scipy.integrate import solve_ivp

result = solve_ivp(dydt,[0.,4.],[y0],t_eval=tspan)
t1 = result.t
y1 = result.y
print(t1,y1)

import pylab
pylab.scatter(t,y,c='b')
pylab.scatter(t1,y1,c='g')
pylab.grid()


romberg	p. 490

import numpy as np

def trap(func,a,b,n=100):
    """
    Composite trapezoidal rule quadrature
    Input:
        func = name of function to be integrated
        a,b = integration limits
        n = number of segments (default = 100)
    Output:
        I = estimate of integral
    """
    if b <= a: return 'upper bound must be greater than lower bound'
    x = a
    h = (b-a)/n
    s = func(a)
    for i in range(n-1):
        x = x + h
        s = s + 2*func(x)
    s = s + func(b)
    I = (b-a)*s/2/n
    return I

def romberg(func,a,b,es=1.e-8,maxit=30):
    """
    Romberg integration quadrature
    input:
        func = name of function to be integrated
        a, b = integration limits
        es = desired relative error (default = 1.e-8)
        maxit = iteration limit (defaul = 30)
    output:
        q = integral estimate
        ea = approximate relative error achieved
        iter = iterations taken
    """
    n = 1
    I = np.zeros((2*maxit,maxit+1))
    I[0,0] = trap(func,a,b,n)
    for iter in range(1,maxit+1):
        n = 2**iter
        I[iter,0] = trap(func,a,b,n)
        for k in range(1,iter+1):
            j = iter-k
            I[j,k] = (4**(k)*I[j+1,k-1] - I[j,k-1])/(4**(k)-1)
        ea = abs((I[0,iter]-I[1,iter-1])/I[0,iter])
        if ea <= es: break
    q = I[0,iter]
    return q,ea,iter

T = 1.

def f(t):
    if t <= T/2:
        i = 8.*np.exp(-t/T)*np.sin(2.*np.pi*t/T)
    else:
        i = 0
    return  i**2

I2val,errel,iter = romberg(f,0.,T,es=0.001)
print(I2val)
print(errel)
print(iter)

Irms = np.sqrt(1/T*I2val)
print('Irms =',Irms)


strlinregr	p. 336

import numpy as np

def strlinregr(x,y):
    n = len(x)
    if len(y) != n: return 'x and y must be of same length'
    sumx = np.sum(x)
    xbar = sumx/n
    sumy = np.sum(y)
    ybar = sumy/n
    sumsqx = 0
    sumxy = 0
    for i in range(n):
        sumsqx = sumsqx + x[i]**2
        sumxy = sumxy + x[i]*y[i]
    a1 = (n*sumxy-sumx*sumy)/(n*sumsqx-sumx**2)
    a0 = ybar - a1*xbar
    e = np.zeros((n))
    SST = 0
    SSE = 0
    for i in range(n):
        e[i] = y[i] - (a0+a1*x[i])
        SST = SST + (y[i]-ybar)**2
        SSE = SSE + e[i]**2
    SSR = SST - SSE
    Rsq = SSR/SST
    SE = np.sqrt(SSE/(n-2))
    return a0,a1,Rsq,SE

x = np.array([10., 20., 30., 40., 50., 60., 70., 80.])
y = np.array([25.,70., 380., 550., 610., 1220., 830., 1450.])
a0,a1,Rsq,SE = strlinregr(x,y)
print('Intercept = {0:7.2f}'.format(a0))
print('Slope = {0:7.3f}'.format(a1))
print('R-squared = {0:5.3f}'.format(Rsq))
print('Standard error = {0:7.2f}'.format(SE))

xline = np.linspace(0,90,10)
yline = a0 + a1*xline
yhat = a0 + a1*x
e = y - yhat
import pylab
pylab.scatter(x,y,c='k',marker='s')
pylab.plot(xline,yline,c='k')
pylab.grid()
pylab.xlabel('x')
pylab.ylabel('y')
pylab.figure()
pylab.hist(e,bins=3,color='w',edgecolor='k',linewidth=2.)
pylab.grid()
pylab.xlabel('Residual')
pylab.figure()
pylab.plot(yhat,e,c='k',marker='o')
pylab.grid()
pylab.xlabel('Predicted y')
pylab.ylabel('Residual')
pylab.title('Residuals vs. Fits')

logx = np.log10(x)
logy = np.log10(y)
a0,a1,Rsq,SE = strlinregr(logx,logy)
print('Intercept = {0:7.2f}'.format(a0))
print('Slope = {0:7.3f}'.format(a1))
print('R-squared = {0:5.3f}'.format(Rsq))
print('Standard error = {0:7.2f}'.format(SE))
logxline = np.linspace(0.9,2,10)
logyline = a0 + a1*logxline
logyhat = a0 + a1*logx
loge = logy - logyhat
pylab.figure()
pylab.scatter(logx,logy,c='k',marker='s')
pylab.plot(logxline,logyline,c='k')
pylab.grid()
pylab.xlabel('log10(x)')
pylab.ylabel('log10(y)')


TableLookup	p. 424

import numpy as np

def TableLookup(x,y,xx):
    n = len(x)
    if n != len(y): return 'input arrays must be the same length'
    if xx < x[0] or xx > x[n-1]:
        return 'input value out of range of table'
    for i in range(n):
        if xx == x[i]:  # check for an exact match
            return y[i]
        elif x[i] > xx:  # check for upper interval
            i2 = i
            break
    xint = (xx-x[i2-1])/(x[i2]-x[i2-1])*(y[i2]-y[i2-1])+y[i2-1]
    return xint


T = np.array([-40, 0., 20., 50., 100., 150., 200., 250.,
              300., 400., 500.])
rho = np.array([1.52, 1.29, 1.20, 1.09, 0.946, 0.935,
                0.746, 0.675, 0.616, 0.525, 0.457])
    
Tx = -501.
rhox = TableLookup(T,rho,Tx)
print(rhox)


TableLookup2		p. 423

import numpy as np

Matl = np.array(['Wheat' , 'Rice' , 'Millet', 'Polyethylene', 
                 'Corn', 'Polystyrene', 'Barley', 'Flaxseeds'])
Prop = np.array(['Absolute Density', 'Bulk Density', 'Percent Void',
                 'Particle Diameter', 'Shape Factor'])
TableData = np.array([[1400., 865., 39.2, 3.61, 1.07],
                     [1457., 905., 37.9, 2.72, 1.04],
                     [1180., 727., 38.4, 1.99, 1.07],
                     [ 922., 592., 35.8, 3.43, 1.02],
                     [1342., 743., 44.6, 7.26, 1.50],
                     [1058., 641., 39.4, 1.56, 1.14],
                     [1279., 725., 43.4, 3.70, 1.14],
                     [1129., 703., 37.8, 2.09, 1.05]])

def TableLookup2(Row,Col,RowNames,ColNames,TableData):
    """
    Function for lookup in a two-dimensional table.
    Input:
        Row = name in row array
        Col = name in column array
        RowNames = row array
        ColNames = column array
        TableData = two-dimensional array of data in table
    Output:
        TableValue = value extracted from table
    """
    n = len(RowNames) ; m = len(ColNames)
    nt = np.size(TableData,0)
    mt = np.size(TableData,1)
    if n != nt or m != mt:
        return 'table information does not conform in size'
    ifind = False
    for i in range(n):
        if Row == RowNames[i]:
            isel = i
            ifind = True
            break
    if not ifind: return 'row name not found in table'
    jfind = False
    for j in range(m):
        if Col == ColNames[j]:
            jsel = j
            jfind = True
            break
    if not jfind: return 'column name not found in table'    
    return TableData[isel,jsel]

MatlName = 'Corn'
PropName = 'BulkDensity'
TableValue = TableLookup2(MatlName,PropName,Matl,Prop,TableData)
if type(TableValue) == type(str()):
    print(TableValue)
else:
    print('{0:1} of {1:11} is {2:7.5g}'.format(PropName,MatlName,TableValue))


TableLookupBin		p. 425

import numpy as np

def TableLookupBin(x,y,xx):
    n = len(x)
    if n != len(y): return 'input arrays must be the same length'
    if xx < x[0] or xx > x[n-1]:
        return 'input value out of range of table'
    iL = 0  ;  iU = n-1
    while True:
        if iU - iL <= 1: break
        iM = int((iL+iU)/2)
        if x[iM] == xx: 
            return y[iM]
        elif x[iM] < xx:
            iL = iM
        else:
            iU = iM
    xint = (xx-x[iL])/(x[iU]-x[iL])*(y[iU]-y[iL])+y[iL]
    return xint


T = np.array([-40, 0., 20., 50., 100., 150., 200., 250.,
              300., 400., 500.])
rho = np.array([1.52, 1.29, 1.20, 1.09, 0.946, 0.935,
                0.746, 0.675, 0.616, 0.525, 0.457])
    
Tx = 350.
rhox = TableLookupBin(T,rho,Tx)
print(rhox)


trap	p. 467

def trap(func,a,b,n=100):
    """
    Composite trapezoidal rule quadrature
    Input:
        func = name of function to be integrated
        a,b = integration limits
        n = number of segments (default = 100)
    Output:
        I = estimate of integral
    """
    if b <= a: return 'upper bound must be greater than lower bound'
    x = a
    h = (b-a)/n
    s = func(a)
    for i in range(n-1):
        x = x + h
        s = s + 2*func(x)
    s = s + func(b)
    I = (b-a)*s/2/n
    return I

g = 9.81  # m/s2
m = 68.1  # kg
cd = 0.25  # kg/m

import numpy as np

def zint(t):
    return np.sqrt(m*g/cd)*np.tanh(np.sqrt(g*cd/m)*t)

z = trap(zint,0.,3.,5)
print(z)


trap_cumulative	p. 474

def trap_cumulative(x,y):
    """
    trapezoidal rule for unequally spaced data
    returns an array of cumulative sums
    Input:
        x = array of independent variable values
        y = array of dependent variable values
        x and y arrays must be of equal length
            and in ascending order of x
    Output:
        s = array of sums
    """
    n = len(x)
    if len(y) != n: return 'x and y arrays must be of equal length'
    for i in range(n-1):
        if x[i+1] < x[i]: return 'x array not in ascending order'
    s = np.zeros((n))
    for k in range(1,n):
        s[k] = s[k-1] + (x[k]-x[k-1])*(y[k]+y[k-1])/2
    return s

import numpy as np

r = np.array([0., 1100., 1500., 2450., 3400., 3630.,
              4500., 5380., 6060., 6280., 6380.])
rho = np.array([13., 12.4, 12.,11.2, 9.7, 5.7,
                5.2, 4.7, 3.6, 3.4, 3.])
n = len(r)
rho2 = rho*1000.  # convert density to kg/m3
r2 = r*1000.  # convert radius to m

mint = 4.*np.pi*r2**2*rho2

m = trap_cumulative(r2,mint)
mt = m/1000.  # convert to tonnes
print('mass of the earth = {0:8.4e} tonnes'.format(mt[n-1]))

vol = 4./3.*np.pi*(r[n-1]*1000.)**3  # cubic meters
den = m[n-1]/vol  # kg/m3
print('average density = {0:7.2f} kg/m3'.format(den))


trapuneq	p. 474

def trapuneq(x,y):
    """
    trapezoidal rule for unequally spaced data
    returns an array of cumulative sums
    Input:
        x = array of independent variable values
        y = array of dependent variable values
        x and y arrays must be of equal length
            and in ascending order of x
    Output:
        s = array of sums
    """
    n = len(x)
    if len(y) != n: return 'x and y arrays must be of equal length'
    for i in range(n-1):
        if x[i+1] < x[i]: return 'x array not in ascending order'
    s = 0
    for k in range(0,n-1):
        s = s + (x[k+1]-x[k])*(y[k+1]+y[k])/2
    return s

import numpy as np

r = np.array([0., 0.12, 0.24, 0.36, 0.49, 0.62, 0.79, 0.86, 0.93, 1.])
rho = np.array([6., 5.81, 5.14, 4.29, 3.39, 2.7, 2.19, 2.1, 2.04, 2.])
r = r / 10.
mint = 4.*np.pi*r**2*rho

m = trapuneq(r,mint)
print('mass estimate = {0:6.4f} g'.format(m))

vol = 4./3./np.pi*(1./10.)**3
den = m/vol
print('average density = {0:6.2f} g/cm3'.format(den))


tridiag		p. 234

import numpy as np

def tridiag(e,f,g,r):
    """
    tridiag: solves a set of n linear algebraic equations
             with a tridiagonal-banded coefficient matris
    input:
    e = subdiagonal vector of length n, first element = 0
    f = diagonal vector of length n
    g = superdiagonal vector of length n, last element = 0
    r = constant vector of length n
    output:
    x = solution vector of length n
    """
    n = len(f)
    # forward elimination
    x = np.zeros([n])
    for k in range(1,n):
        factor = e[k]/f[k-1]
        f[k] = f[k] - factor*g[k-1]
        r[k] = r[k] - factor*r[k-1]
    # back substitution
    x[n-1] = r[n-1]/f[n-1] 
    for k in range(n-2,-1,-1):
        x[k] = ( r[k] - g[k]*x[k+1] )/f[k]
    return x

n = 4
e = np.zeros([n])
f = np.zeros([n])
g = np.zeros([n])
for i in range(n):
    f[i] = 2.04
    if i < n-1:
        g[i] = -1
    if i > 0 :
        e[i] = -1
r = np.array(([40.8],[0.8],[0.8],[200.8]))

T = tridiag(e,f,g,r)
print('Temperatures in degC are:')
for i in range(4):
    print('       {0:6.2f}'.format(T[i]))


velocity1	p. 74

def velocity1(dt,ti,tf,vi):
    """
    Solution of bungee jumper velocity
    by Euler's method
    Input:
        dt = time step (s)
        ti = initial time (s)
        tf = final time (s)
        vi = initial value of velocity (m/s)
    output:
        vf = velocity at tf (m/s)
    """
    t = ti
    v = vi
    n = int((tf-ti)/dt)
    for i in range(n):
        dvdt = deriv(v)
        v = v + dvdt*dt
        t = t + dt
    vf = v
    return vf

def deriv(v):
    g = 9.81
    m = 68.1
    cd = 0.25
    dv = g - cd/m * v * abs(v)
    return dv


wegstein	p. 148

def wegstein(g,x0,x1,Ea=1.e-7,maxit=30):
    """
    This function solves x=g(x) using the Wegstein method.
    The method is repeated until either the relative error
    falls below Ea (default 1.e-7) or reaches maxit (default 30).
    Input:
        g = name of the function for g(x)
        x0 = first initial guess for x
        x1 = second initial guess for x
        Ea = relative error threshold
        maxit = maximum number of iterations
    Output:
        x2 = solution estimate
        ea = relative error
        i+1 = number of iterations
    """    
    for i in range(maxit):
        x2 = (x1*g(x0)-x0*g(x1))/(x1-x0-g(x1)+g(x0))
        ea = abs((x1-x0)/x1)
        if ea < Ea:  break
        x0 = x1
        x1 = x2
    return x2,ea,i+1

import numpy as np

def g(x):
    return np.exp(-x)

x0 = 0
x1 = 0.25
(xsoln,ea,n) = wegstein(g,x0,x1,Ea=1.e-5)
print('Solution = {0:8.5g}'.format(xsoln))
print('Relative error = {0:8.3e}'.format(ea))
print('Number of iterations = {0:5d}'.format(n))


Case Study 15.6		p. 369

import numpy as np
import pylab
# create the data arrays
U = np.array([0.5, 2., 10., 0.5, 2., 10., 0.5, 2., 10.])
H = np.array([0.15, 0.15, 0.15, 0.3, 0.3, 0.3, 0.5, 0.5, 0.5])
KL = np.array([0.48, 3.9, 57., 0.85, 5., 77., 0.8, 9., 92.])
# log transformation
logU = np.log10(U)
logH = np.log10(H)
logKL = np.log10(KL)
n = len(U)
# calculate the X matrix
X = np.zeros((n,3))
for i in range(n):
    X[i,0] = 1
    X[i,1] = logU[i]
    X[i,2] = logH[i]
# formulate and solve the normal equations
Xt = np.transpose(X)
A = np.dot(Xt,X)
const = np.dot(Xt,np.transpose(logKL))
b = np.linalg.solve(A,const)
print('Estimated model parameters are:\n',b)
# back-transform the intercept
b0 = 10**b[0]
print('b0 = ',b0)
 # model predictions
logKLpred =  b[0] + b[1]*logU + b[2]*logH
# residuals and SSE
e = logKL - logKLpred
SSE = np.dot(e,e)
print('SSE = ',SSE)
# SST
SST = np.var(logKL)*(n-1)
# R2 and se
R2 = 1 - SSE/SST
se = np.sqrt(SSE/(n-3))
print('R-sqaured = ',R2)
print('Standard error = ',se)
# residuals versus fits plot
pylab.scatter(logKLpred,e,c='k',marker='s')
pylab.grid()
pylab.xlabel('Predicted logKL Values')
pylab.ylabel('Residual Error')
# performance plots
pylab.figure()
pylab.scatter(logKL,logKLpred,c='k',marker='o')
pylab.plot([-0.4,2.1],[-0.4,2.1],c='k',ls='--')
pylab.grid()
pylab.xlabel('Measured logKL')
pylab.ylabel('Predicted logKL')
pylab.figure()
KLpred = 10**logKLpred
pylab.scatter(KL,KLpred,c='k',marker='o')
pylab.plot([0.,100.],[0.,100.],c='k',ls='--')
pylab.grid()
pylab.xlabel('Measured KL')
pylab.ylabel('Predicted KL')


Case Study 16.7		p. 392

import numpy as np
import pylab

yr,numspots,sd,n1,n1 = np.loadtxt(fname='SN_y_tot_V2.0.csv',delimiter=';',unpack=True)
pylab.plot(yr,numspots,c='k')
pylab.grid()

coef = np.polyfit(yr,numspots,1)
pylab.plot(yr,np.polyval(coef,yr),c='k',ls='--')

y = numspots - np.polyval(coef,yr)
pylab.figure()
pylab.plot(yr,y,c='k')
pylab.grid()
pylab.xlabel('Year')
pylab.ylabel('Sunspot Number')

from scipy.fft import fft
Y = fft(y)
fs = 1  # 1/yr
n = len(yr)
f = np.arange(1,n)*fs/n
n2 = int(n/2)
f2 = f[0:n2]
Y2 = Y[1:n2+1]
Pyy = abs(Y2)**2
pylab.figure()
pylab.plot(f2,Pyy,c='k')
pylab.grid()
pylab.xlabel('Frequency - cycles/yr')
pylab.ylabel('Power')

pmax = np.max(Pyy)
for i in range(n2):
    if Pyy[i] >= pmax:
        imax = i
        fmax = f2[i]
        break

print('Frequency at max power = ',fmax,' 1/yr')
print('Period at max power = ',1/fmax,' years')   


Case Study 18.8		p. 448

z = np.array([0., 2.3, 4.9, 9.1, 13.7, 18.3, 22.9, 27.2])
T = np.array([22.8, 22.8, 22.8, 20.6, 13.9, 11.7, 11.1, 11.1])

zz = 10.
TT,b,c,d = cspline(z,T,zz)

zplot = np.linspace(0.,27.2)
n = len(zplot)
m = len(z)
Tplot = np.zeros((n))
dTplot = np.zeros((n))
d2Tplot = np.zeros((n))
for i in range(n):
    for j in range(m-1):
        if z[j] > zplot[i]:
            j2 = j-1
            break
    Tplot[i] = T[j2]+b[j2]*(zplot[i]-z[j2])+c[j2]*(zplot[i]-z[j2])**2 \
    + d[j2]*(zplot[i]-z[j2])**3
    dTplot[i] = b[j2]+2*c[j2]*(zplot[i]-z[j2])+3*d[j2]*(zplot[i]-z[j2])**2
    d2Tplot[i] = 2*c[j2]+6*d[j2]*(zplot[i]-z[j2])

import matplotlib.pyplot as plt
fig = plt.figure()
ax1 = fig.add_subplot(131)
fig.subplots_adjust(wspace=0.5)
ax1.plot(Tplot,zplot,c='k')
ax1.scatter(T,z,c='k',marker='s')
ax1.set_ylim(30.,0.)
ax1.set_title('T vs. z')
ax1.set_ylabel('Depth, z in m')
ax1.set_xlabel('degC')
ax1.grid()

ax2 = fig.add_subplot(132)
ax2.plot(dTplot,zplot,c='k')
ax2.set_ylim(30.,0.)
ax2.set_title('dT/dz vs. z')
ax2.set_xlabel('degC/m')
ax2.grid()

ax3 = fig.add_subplot(133)
ax3.plot(d2Tplot,zplot,c='k')
ax3.set_ylim(30.,0.)
ax3.set_title('d2T/dz2 vs. z')
ax3.set_xlabel('degC/m2')
ax3.grid()


Case Study 22.6		p. 553

import numpy as np
import matplotlib.pyplot as plt

def eulersys(dydt,tspan,y0,h=-1.,*args):
    """
    Euler method for solving a system of ODEs
    input:
        dydt = function name that evaluates the derivatives
        tspan = array of independent variable values where either
            ti and tf are the initial and final values
            of the independent variable when h is specified,
            or the array specifies the values of t for
            solution (h is not specified)
        y0 = initial value of the dependent variable
        h = step size, default = -1.
        *args = additional argument to be passed to dydt
    output:
        t = array of independent variable values
        y = array of dependent variable values
    """
    if np.any(np.diff(tspan) < 0): return 'tspan times must be ascending'
    # check if only ti and tf spec'd and no value for h
    if len(tspan) == 2 and h != -1.:  
            ti = tspan[0] ; tf = tspan[1]
            nsteps = int((tf-ti)/h)
            t = []
            t.append(ti)
            for i in range(nsteps):  # add the rest of the t values
                t.append((i+1)*h)
            n = len(t)    
            if t[n-1] < tf:  # check if t array is short of tf
                t.append(tf)
                n = n+1
    else:
        n = len(tspan)  # here if tspan contains step times
        t = tspan
    neq = len(y0)
    y = np.zeros((n,neq))  # set up 2-D array for dependent variables
    for j in range(neq):
        y[0,j] = y0[j]  #  set first elememts to initial conditions
    for i in range(n-1):  # 4th order RK
        hh = t[i+1] - t[i]
        dy = dydt(t[i],y[i,:],*args)
        y[i+1,:] = y[i,:] + dy*hh
    return t,y

def rk4sys(dydt,tspan,y0,h=-1.,*args):
    """
    fourth-order Runge-Kutta method
    for solving a system of ODEs
    input:
        dydt = function name that evaluates the derivatives
        tspan = array of independent variable values where either
            ti and tf are the initial and final values
            of the independent variable when h is specified,
            or the array specifies the values of t for
            solution (h is not specified)
        y0 = initial value of the dependent variable
        h = step size, default = -1.
        *args = additional argument to be passed to dydt
    output:
        t = array of independent variable values
        y = array of dependent variable values
    """
    if np.any(np.diff(tspan) < 0): return 'tspan times must be ascending'
    # check if only ti and tf spec'd and no value for h
    if len(tspan) == 2 and h != -1.:  
            ti = tspan[0] ; tf = tspan[1]
            nsteps = int((tf-ti)/h)
            t = []
            t.append(ti)
            for i in range(nsteps):  # add the rest of the t values
                t.append((i+1)*h)
            n = len(t)    
            if t[n-1] < tf:  # check if t array is short of tf
                t.append(tf)
                n = n+1
    else:
        n = len(tspan)  # here if tspan contains step times
        t = tspan
    neq = len(y0)
    y = np.zeros((n,neq))  # set up 2-D array for dependent variables
    for j in range(neq):
        y[0,j] = y0[j]  #  set first elememts to initial conditions
    for i in range(n-1):  # 4th order RK
        hh = t[i+1] - t[i]
        k1 = dydt(t[i],y[i,:],*args)
        ymid = y[i,:] + k1*hh/2.
        k2 = dydt(t[i]+hh/2.,ymid,*args)
        ymid = y[i,:] + k2*hh/2.
        k3 = dydt(t[i]+hh/2.,ymid,*args)
        yend = y[i,:] + k3*hh
        k4 = dydt(t[i]+hh,yend,*args)
        phi = (k1 + 2.*(k2+k3) + k4)/6.
        y[i+1,:] = y[i,:] + phi*hh
    return t,y

def predprey(t,y,a,b,c,d):
    dy = np.zeros((2))
    dy[0] = a*y[0] - b*y[0]*y[1]
    dy[1] = -c*y[1] + d*y[0]*y[1]
    return dy

h = 0.0625
tspan = np.array([0.,40.])
y0 = np.array([2.,1.])
a = 1.2 ; b = 0.6 ; c = 0.8 ; d = 0.3

t,y = eulersys(predprey,tspan,y0,h,a,b,c,d)
fig = plt.figure()
ax1 = fig.add_subplot(221)
ax1.plot(t,y[:,0],c='k',label='prey')
ax1.plot(t,y[:,1],c='k',ls='--',label='predator')
ax1.grid()
ax1.set_xlabel('a) Euler time plot')
ax2 = fig.add_subplot(222)
ax2.plot(y[:,0],y[:,1],c='k')
ax2.grid()
ax2.set_xlabel('b) Euler phase plot')

t,y = rk4sys(predprey,tspan,y0,h,a,b,c,d)
ax3 = fig.add_subplot(223)
ax3.plot(t,y[:,0],c='k',label='prey')
ax3.plot(t,y[:,1],c='k',ls='--',label='predator')
ax3.grid()
ax3.set_xlabel('c) RK4 time plot')
ax4 = fig.add_subplot(224)
ax4.plot(y[:,0],y[:,1],c='k')
ax4.grid()
ax4.set_xlabel('d) RK4 phase plot')

plt.subplots_adjust(hspace=0.5)


Example 14.4		p. 319

import numpy as np
import pylab

x = np.array([10.,20.,30.,40.,50.,60.,70.,80.])
y = np.array([25.,70.,380.,550.,610.,1220.,830.,1450.])

sumx = np.sum(x)
sumy = np.sum(y)

n = len(x)
sumxy = 0
sumsqx = 0
for i in range((n)):
    sumxy = sumxy + x[i]*y[i]
    sumsqx = sumsqx + x[i]**2

ybar = np.mean(y)
xbar = np.mean(x)
a1 = (n*sumxy - sumx*sumy)/(n*sumsqx -sumx**2)
a0 = ybar - a1*xbar

yhat = a0 + a1*x
e = y - yhat

pylab.plot(x,e,c='k',marker='o')
pylab.grid()
pylab.xlabel('x')
pylab.ylabel('e')


Example 16.3		p. 389

import numpy as np
import pylab
from scipy.fft import fft

n = 8 ; dt = 0.02 ; fs = 1/dt ; T = 0.16
tspan = np.arange(0,n)/fs
y = 5 + np.cos(2*np.pi*12.5*tspan) + np.sin(2*np.pi*18.75*tspan)
pylab.plot(tspan,y,c='k',marker='o')
pylab.grid()
pylab.xlabel('Time - s')
pylab.ylabel('f(t)')

np.set_printoptions(precision=3,suppress=True)
Y = fft(y)/n
print(Y)

import matplotlib.pyplot as plt
Ymag = abs(Y[1:])
freq = np.arange(1,n)/T
fig=plt.figure()
plt.stem(freq,Ymag,basefmt='k-',linefmt='k-',markerfmt='ko')
plt.xlabel('Frequency - Hz')
plt.ylabel('DFT magnitude')
plt.grid()


Example 16.4		p. 391

import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft

# compute the DFT
n = 8 ; dt = 0.02 ; fs = 1/dt ; T = 0.16
tspan = np.arange(0,n)/fs
y = 5 + np.cos(2*np.pi*12.5*tspan) + np.sin(2*np.pi*18.75*tspan)
Y = fft(y)/n
freq = np.arange(1,n)/T
# compute and plot the power spectrum
nyquist = fs/2
n2 = int(n/2)
fP = freq[0:n2]
Pyy = abs(Y[1:n2+1])**2
fig=plt.figure()
plt.stem(fP,Pyy,basefmt='k-',linefmt='k-',markerfmt='ko',
         use_line_collection=True)
plt.xlabel('Frequency (Hz)')
plt.title('Power Spectrum')
plt.grid()


Example 17.7		p. 412

import numpy as np
import pylab
# generate 5 equally-spaced points and function values
x = np.linspace(-1.,1.,5)
y = 1./(1.+25.*x**2)
# 50 interpolation poings
xx = np.linspace(-1.,1.)
# fit 4th-order polynomial
coef = np.polyfit(x,y,4)
# use polynomial to interpolate
y4 = np.polyval(coef,xx)
# Runge's function values
yr = 1./(1.+25.*xx**2)
# generate plot
pylab.scatter(x,y,c='k',marker='o')
pylab.plot(xx,y4,c='k',ls='--',label='4th-order')
pylab.plot(xx,yr,c='k',label='Runge')
pylab.grid()
pylab.xlabel('x')
pylab.ylabel('y')
pylab.legend()


Figure 15.5	p. 362

import numpy as np
import pylab
# set table date
x = np.array([0., 4.44, 10., 15.56, 21.11, 26.67, 32.22,
              37.78, 48.89, 60., 71.11, 82.22, 93.33])
y = np.array([1.794, 1.546, 1.31, 1.129, 0.982, 0.862,
              0.764, 0.682, 0.559, 0.47, 0.401, 0.347, 0.305])
# standardize x to z
xbar = np.mean(x)
sx = np.std(x)
z = (x-xbar)/sx
n = len(x)
m = 9
  # order of polynomial
# calculate the X matrix
X = np.zeros((n,m+1))
for i in range(n):
    for j in range(m+1):
        X[i,j] = z[i]**j
# formulate and solve the normal equations
Xt = np.transpose(X)
A = np.dot(Xt,X)
const = np.dot(Xt,np.transpose(y))
b = np.linalg.solve(A,const)
print('Estimated model parameters are:\n',b)
#compute model predictions and residuals
yhat = np.zeros((n))
e = np.zeros((n))
for i in range(n):
    for j in range(m+1):
        yhat[i] = yhat[i] + b[j]*z[i]**j
    e[i] = y[i] - yhat[i]
#compute the sum of squares
SSE = np.dot(e,e)  # residuals
SST = np.var(y)*(n-1)  # total corrected
SSR = SST - SSE   # regression
print('SSE = ',SSE)
print('SST = ',SST)
print('SSR = ',SSR)
# R2 and the standard error of the estimate
R2 = SSR/SST
se2 = SSE/(n-(m+1))
se = np.sqrt(se2)
print('R-squared = ',R2)
print('Standard error of the estimate = ',se)
# adjusted R2
R2adj = 1 - (SSE/(n-(m+1)))/(SST/(n-1))
print('Adjusted R-squared = ',R2adj)
# covariance matrix and parameter standard errors
XtXinv = np.linalg.inv(np.dot(Xt,X))
covb = XtXinv*se2
seb = np.zeros((m+1))
for j in range(m+1):
    seb[j] = np.sqrt(se2*covb[j,j])
print('Standard errors of the parameter estimates:\n',seb)
# hat matrix, PRESS statistic and predicted R2
H = np.dot(np.dot(X,XtXinv),Xt)
PRESS = 0
for i in range(n):
    PRESS = PRESS + (e[i]/(1-H[i,i]))**2
R2pred = 1 - PRESS/SST
print('PRESS = ',PRESS)
print('Predicted R-squared = ',R2pred)
# plots
def ypred(x,xavg,xstd,b,m):
    yp = 0
    zvar = (x-xavg)/xstd
    for j in range(m+1):
        yp = yp + b[j]*zvar**j
    return yp

xplot = np.linspace(0,95.)
npt = len(xplot)
yplot = np.zeros((npt))
for k in range(npt):
    yplot[k] = ypred(xplot[k],xbar,sx,b,m) 
pylab.scatter(x,y,c='k',marker='s')
pylab.plot(xplot,yplot,c='k')
pylab.grid()
pylab.xlabel('Temperature - degC')
pylab.ylabel('Viscosity - cP')
pylab.figure()
pylab.plot(yhat,e,c='k',ls='-',marker='s')
pylab.grid()
pylab.xlabel('Predicted Viscosity - cP')
pylab.ylabel('Residual - cP')


SSE plus code		p. 366

def SSE(params):
    A = params[0] ; B = params[1] ; C = params[2]
    logVPpred = A - B/(C+T)
    e = logVP - logVPpred
    return np.dot(e,e)

from scipy.optimize import minimize
result = minimize(SSE,(A,B,C))
print(result.x)

# plot data
pylab.figure()
pylab.scatter(T,VP,c='k',marker='.',label='data')
pylab.grid()
pylab.xlabel('Temperature - degC')
pylab.ylabel('Vapor Pressure - torr')
# regressed estimates for parameters
A = result.x[0]
B = result.x[1]
C = result.x[2]
# plot model curve
Tplot = np.linspace(30.,300.,100)
VPplot = 10**(A-B/(C+Tplot))
pylab.plot(Tplot,VPplot,c='k',ls=':',label='model')
pylab.legend()
pylab.title('Nonlinear Regression')

pylab.figure()
pylab.scatter(T,logVP,c='k',marker='.',label='data')
pylab.grid()
pylab.xlabel('Temperature - degC')
pylab.ylabel('log10(VP)')
pylab.title('Nonlinear Regression')

# regression statistics
SSe = SSE(result.x)
n = len(T)
SST = np.var(logVP)*(n-1)
R2 = 1 - SSe/SST
se = np.sqrt(SSe/(n-3))
print('R-squared = ',R2)
print('Standard error of the estimate = ',se)

logVPpred = A - B/(C+T)
e = logVP - logVPpred
pylab.figure()
pylab.plot(logVPpred,e,c='k')
pylab.grid()
pylab.xlabel('Predicted log10(VP)')
pylab.ylabel('Residual')
pylab.title('Residuals vs. Fits')

# formulate the X matrix
X = np.zeros((n,3))
for i in range(n):
    X[i,0] = 1
    X[i,1] = 1/(C+T[i])
    X[i,2] = - B/(C+T[i])**2
# compute the covariance matrix
XtXinv = np.linalg.inv(np.dot(np.transpose(X),X))
cov = se**2 * XtXinv
seA = cov[0,0]
seB = cov[1,1]
seC = cov[2,2]
print('Estimated parameter standard errors:')
print('A: ',seA)
print('B: ',seB)
print('C: ',seC)

# coefficients of variation
print('Coefficients of variation (%):')
print('A: {0:7.3g}'.format(seA/A*100))
print('B: {0:7.3g}'.format(seB/B*100))
print('C: {0:7.3g}'.format(seC/C*100))


from dgeqp3 import dgeqp3
from numpy import *
import scipy
import scipy.linalg


def qrp(a, lwork=None):
    """Compute QR decomposition of a matrix."""

    a1 = copy(a) # Avoid overwriting
    a1 = asarray_chkfinite(a1)
    if len(a1.shape) != 2:
        raise ValueError("expected 2D array")
    M, N = a1.shape
    jpvt = zeros(N)

    if lwork is None or lwork == -1:
        # Get optimal work array
        work = zeros(1)
        qr,jpvt,tau,work,info = dgeqp3(M,N,a1,jpvt,work,-1)
        lwork = work[0]

    work = zeros(lwork)
    qr,jpvt,tau,work,info = dgeqp3(M,N,a1,jpvt,work,lwork)

    if info<0:
        raise ValueError("illegal value in %-th argument of internal geqp3"
            % -info)

    # Setup Q from reflectors in qr
    Q,work,info = scipy.linalg.flapack.dorgqr(qr[:,0:M],tau)

    R = triu(qr)

    return Q,R,jpvt-1  # The -1 because indexing in FORTRAN begins at 1 instead of 0

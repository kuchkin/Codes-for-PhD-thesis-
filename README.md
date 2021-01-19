Solving the skyrmion equation  

	rT'' + T' = Sin(2T)/(2r) - 4\pi Sin(T)^2 + 4\pi^2 r (h Sin(T) + u Sin(2T)),
	T(r = 0) = \pi,
	T(r = L) = 0,

by quasilinearization (QL) method is given.

The initial guess for skyrmion profile is chosen to be
	T0 = 2ArcTan(Exp[-2\pi\sqrt(h+2u)r]/r),

The linearized equation has the form

	rT'' + T' + f1(r)T = f2(r),

Applying the second order finite difference scheme for derivatives leads to the system of linear equations:

	A_{i,j} * T_{i} = b_{j},

with nonzero elemnts A_{i,i}, A_{i,i+1}, A{i+1,i}.
The system is solved by Gauss elimination method with consequential back substitution.
The error is given by 

	Error = \Int (T-T0)^2 dr, r \in (0,L).

from __future__ import print_function
from orbits import *
# Compute the number of components of a curve on a genus 2 surface
# split along a separating curve bounding two one-holed tori.
# We assume there are three parallels bands, of weights p, q, r, on each
# one-holed torus, and a shift of n on the circle.

def PG(p,q,r,n):
   N = p+q+r
   n = n%(2*N)
   G= Pseudogroup([ Flip([1,p],[p+q+r+1, 2*p+q+r]),
		    Flip([p+1, p+q], [2*p+q+r+1, 2*p+2*q+r]),
		    Flip([p+q+1, p+q+r], [2*p+2*q+r+1, 2*p+2*q+2*r]),
		    Flip([2*N+1,2*N+r],[2*N+p+q+r+1, 2*N+p+q+2*r]),
		    Flip([2*N+r+1, 2*N+r+q], [p+q+2*r+1, p+2*q+2*r]),
		    Flip([2*N+r+q+1, 2*N+p+q+r], [2*N+p+2*q+2*r+1, 2*N+2*p+2*q+2*r]),
		    Shift([1, 2*N - n],[2*N+n+1,4*N]),
                    Shift([2*N-n+1,2*N],[2*N+1, 2*N+n])], 
                Interval(1,4*N))
   return G.reduce()

if __name__ == '__main__':
   print('p = q = r = n = 1:')
   print(PG(1,1,1,1), 'components')

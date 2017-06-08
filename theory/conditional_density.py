import sympy as sp

from sympy.abc import a,b, x, y
from sympy import exp

p1= sp.integrate(  1, (x, b-y -a, b-y + a))

print p1
print p1/sp.integrate(  p1, (y,b-a,b+a))

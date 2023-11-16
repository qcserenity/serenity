#!/usr/bin/python3
#@file   sharmonics.py
#
#@date   May 29, 2017
#@author Jan Unsleber
#@copyright \n
# This file is part of the program Serenity.\n\n
# Serenity is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.\n\n
# Serenity is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.\n\n
# You should have received a copy of the GNU Lesser General
# Public License along with Serenity.
# If not, see <http://www.gnu.org/licenses/>.\n

# ======================================================================
#  This piece of python code generates c++ or latex code for
#   the evaluation of spherical basisfunctions at points in 
#   space in a cartesian representation.
#  It also gives the code for first and second derivatives with
#   rescpect to Cartesian coordinates.
#  The code uses the Herglotz generating function including separation
#   into z+r and x+y factors for real harmonics in Cartesian space.
#  Then medium efforts to reduce the complexity of the expressions is
#   made.
#  Most major optimizations for computational time should have been 
#   included.
#     - Jan Unsleber
# ======================================================================

# ===========
#   Imports
# ===========
from functools import reduce
import fractions
from math import *
import sys

# =============
#   Functions
# =============
def pascal(rows):
    """
     A function for Pascal's triangle.
    """
    result = []
    for rownum in range (rows):
        newValue=1
        l = [newValue]
        for iteration in range (rownum):
            newValue = newValue * ( rownum-iteration ) * 1 / ( iteration + 1 )
            l.append(int(newValue))
        result.append(l)
    return result

def gcd(l):
  """
   A function for the 'greatest common divisor' of more than 3 integers
  """
  r = fractions.gcd(l[0],l[1])
  for i in range(2,len(l)):
    r = fractions.gcd(r,l[i])
  return r

# ==============================
#   Term : n * x^a * y^b * z^c
# ==============================
class Term:
  """
  A 'Term' is an experssion of there form n * x^a * y^b * z^c
    where only n,a,b,c are needed for all operatrions.
  """
  def __init__(self,n,a,b,c):  
    self._n = n
    self._a = a
    self._b = b
    self._c = c

  def isNumber(self):
    return (self._n==self._n+self._a+self._b+self._c)
 
  def __mul__(self,other):
    return Term(self._n*other._n,self._a+other._a,self._b+other._b,self._c+other._c)
  __rmul__ = __mul__

  def __add__(self,other):
    if not self.combinable(other):
      print("Error: tried to combine terms that can not be added.")
      print(self)
      print(other)
    return Term(self._n+other._n,self._a,self._b,self._c)

  def __str__(self):
    if self._n==0: return "0";
    string = "{:-.1f}".format(self._n)
    if (self._a>0): string += "x^{:d}".format(self._a)
    if (self._b>0): string += "y^{:d}".format(self._b)
    if (self._c>0): string += "z^{:d}".format(self._c)
    return string

  def short(self):
    if self._n==0: return 0;
    string = "{: .1f}".format(self._n)
    string += ":{:d}".format(self._a)
    string += ":{:d}".format(self._b)
    string += ":{:d}".format(self._c)
    return string
  
  def latex(self):
    return self.__str__()
 
  def cpp(self):
    if self._n==0: return 0;
    string = "{:.1f}".format(self._n)
    for i in range(self._a): string += "*x";
    for i in range(self._b): string += "*y";
    for i in range(self._c): string += "*z";
    return string

  def cppArray(self):
    if self._n==0: return "0";
    string = "{:.1f}".format(self._n)
    if self._a>0:string += "*x[{:d}]".format(self._a)
    if self._b>0:string += "*y[{:d}]".format(self._b)
    if self._c>0:string += "*z[{:d}]".format(self._c)
    return string

  def dx(self):
    return Term(self._n*self._a,self._a-1,self._b,self._c)

  def dy(self):
    return Term(self._n*self._b,self._a,self._b-1,self._c)

  def dz(self):
    return Term(self._n*self._c,self._a,self._b,self._c-1)

  def combinable(self,other):     
    return self._a==other._a and self._b==other._b and self._c==other._c

# ===========================
#   Harmonic Function Y^m_l
# ===========================

class CRSphHarmon:
  def __init__(self,l,m,calc=True):
    self._l = l
    if (abs(m)>l): 
      print("Error |m| > l")
      exit()
    self._m = m
    self._terms = []
    if calc:
      self._pt = pascal(l+1)
  
      # prefactor
      numers = []
      denoms = []
      for k in range(0,int((l-abs(m))/2)+1):
        n = 1 if k%2==0 else -1
        numer = int(factorial(2*l-2*k)/(factorial(k)*factorial(l-k)*factorial(l-2*k-abs(m))))
        denom = int(pow(2,l))
        while (numer%2==0 and denom%2==0):
          numer = int(numer/2)
          denom = int(denom/2)
        numers.append(numer)
        denoms.append(denom)
      maxdenom = max(denoms)
      for i in range(len(numers)):
        numers[i] = int(numers[i]*maxdenom/denoms[i])
     
      if len(numers)>1:
        gcf = gcd(numers)
      else:
        gcf = numers[0]
      for i in range(len(numers)):
        numers[i] = int(numers[i]/gcf)  
  
      # final prefactor sqrt
      if m!=0:
        self._numer = 2*factorial(l-abs(m))
        self._denom =   factorial(l+abs(m))
      else:
        self._denom = 1 
        self._numer = 1
      
      self._denom *= maxdenom*maxdenom  
      self._numer *= gcf*gcf  
  
      gcf = fractions.gcd(self._denom,self._numer)  
      self._denom = int(self._denom/gcf)
      self._numer = int(self._numer/gcf)
    
      # z-r parts
      zrs=[]
  
      for k in range(0,int((l-abs(m))/2)+1):
        z = Term(1.0,0,0,l-2*k-abs(m))
        rs = []
        pref = numers[k]
        if    2*k == 2:
          rs.append(Term(1.0*pref,2,0,0))
          rs.append(Term(1.0*pref,0,2,0))
          rs.append(Term(1.0*pref,0,0,2))
        elif  2*k == 4:
          rs.append(Term(1.0*pref,4,0,0))
          rs.append(Term(1.0*pref,0,4,0))
          rs.append(Term(1.0*pref,0,0,4))
          rs.append(Term(2.0*pref,2,2,0))
          rs.append(Term(2.0*pref,0,2,2))
          rs.append(Term(2.0*pref,2,0,2))
        elif  2*k == 6:
          rs.append(Term(1.0*pref,6,0,0))
          rs.append(Term(1.0*pref,0,6,0))
          rs.append(Term(1.0*pref,0,0,6))
          rs.append(Term(3.0*pref,4,2,0))
          rs.append(Term(3.0*pref,4,0,2))
          rs.append(Term(3.0*pref,2,4,0))
          rs.append(Term(3.0*pref,0,4,2))
          rs.append(Term(3.0*pref,0,2,4))
          rs.append(Term(3.0*pref,2,0,4))
          rs.append(Term(6.0*pref,2,2,2))
        elif  2*k == 8:
          rs.append(Term( 1.0*pref,8,0,0))
          rs.append(Term( 1.0*pref,0,8,0))
          rs.append(Term( 1.0*pref,0,0,8))
          rs.append(Term( 4.0*pref,6,2,0))
          rs.append(Term( 4.0*pref,6,0,2))
          rs.append(Term( 4.0*pref,2,6,0))
          rs.append(Term( 4.0*pref,0,6,2))
          rs.append(Term( 4.0*pref,0,2,6))
          rs.append(Term( 4.0*pref,2,0,6))
          rs.append(Term( 6.0*pref,4,4,0))
          rs.append(Term( 6.0*pref,4,0,4))
          rs.append(Term( 6.0*pref,0,4,4))
          rs.append(Term(12.0*pref,4,2,2))
          rs.append(Term(12.0*pref,2,4,2))
          rs.append(Term(12.0*pref,2,2,4))
        elif 2*k == 10:
          rs.append(Term( 1.0*pref,10,0,0))
          rs.append(Term( 1.0*pref,0,10,0))
          rs.append(Term( 1.0*pref,0,0,10))
          rs.append(Term( 5.0*pref,8,2,0))
          rs.append(Term( 5.0*pref,8,0,2))
          rs.append(Term( 5.0*pref,2,8,0))
          rs.append(Term( 5.0*pref,0,8,2))
          rs.append(Term( 5.0*pref,0,2,8))
          rs.append(Term( 5.0*pref,2,0,8))
          rs.append(Term(10.0*pref,6,4,0))
          rs.append(Term(10.0*pref,6,0,4))
          rs.append(Term(20.0*pref,6,2,2))
          rs.append(Term(10.0*pref,4,6,0))
          rs.append(Term(10.0*pref,0,6,4))
          rs.append(Term(20.0*pref,2,6,2))
          rs.append(Term(10.0*pref,4,0,6))
          rs.append(Term(10.0*pref,0,4,6))
          rs.append(Term(20.0*pref,2,2,6))
          rs.append(Term(30.0*pref,4,4,2))
          rs.append(Term(30.0*pref,2,4,4))
          rs.append(Term(30.0*pref,4,2,4))
        elif (2*k>10):
          print("Error r^x with x>10 required, this is not implemented.")
          exit()
        else:
          rs.append(Term(1.0*pref,0,0,0))
                 
        for r in rs:
          if k%2!=0: r._n = -1.0 * r._n
          zrs.append(z*r)
  
      # x-y parts
      xys = []
      if m>0:
        for k in range(0,m+1):
  
          # m-k=exp for x, k = exp for y
          if k%2==0:
            pref = self._pt[m][k] if k%4==0 else -self._pt[m][k]
            xys.append(Term(pref,m-k,k,0))
      elif m==0:
        xys.append(Term(1.0,0,0,0))
      elif m<0:
        for k in range(0,abs(m)+1):
          # m-k=exp for x, k = exp for y
          pref = self._pt[abs(m)][k] if (k-1)%4==0 else -self._pt[abs(m)][k]
          if k%2==1:
            xys.append(Term(pref,abs(m)-k,k,0))
      
      terms = []
      for xy in xys:
        for zr in zrs:
          terms.append(xy*zr)
      
      while len(terms)>0:    
        done = [0]
        current = terms[0] 
        for j in range(i+1,len(terms)):   
          if current.combinable(terms[j]):
            current = current + terms[j]
            done.append(j)
        if current._n != 0.0:
          self._terms.append(current)
        for d in done[::-1]:   
          terms.pop(d)
       

  def __str__(self):
    string =  "Y^{{{:2d}}}_{{{:2d}}}(x,y,z) &=& \sqrt{{\\frac{{{:3d}}}{{4\pi}}}} ".format(self._m,self._l,2*self._l+1)
    string += "\sqrt{{ \\frac{{{:-d}}}{{{:-d}}} }}".format(self._numer,self._denom)
    for t in range(len(self._terms)):
      if t!=0: string += " +"
      string += " ("+self._terms[t].latex()+")"
    return string

  def latex(self):
    return self.__str__()

  def cpp(self):
    string = "sqrt({:-.1f}/{:-.1f}) *".format(self._numer,self._denom)
    for t in range(len(self._terms)):
      if t!=0: string += " +"
      string += " "+self._terms[t].cpp()
    return string

  def cppArray(self):
    if self._numer==1.0 and self._denom == 1.0: string = "("
    elif self._denom==1.0: string = "sqrt({:-.1f}) * (".format(self._numer)
    else: string = "sqrt({:-.1f}/{:-.1f}) * (".format(self._numer,self._denom)
    copy = string
    for t in range(len(self._terms)):
      if t!=0 and self._terms[t]._n>0: string += "+"
      string += self._terms[t].cppArray()
    if string==copy:
      return "0.0";
    return string+")"
       
  def dx(self):
    new = CRSphHarmon(self._l,self._m,False)
    new._numer = self._numer
    new._denom = self._denom
    new._terms = [] 
    for t in self._terms:
      dx = t.dx()
      if dx._n != 0:
        new._terms.append(dx)
    return new

  def dy(self):
    new = CRSphHarmon(self._l,self._m,False)
    new._numer = self._numer
    new._denom = self._denom
    new._terms = [] 
    for t in self._terms:
      dy = t.dy()
      if dy._n != 0:
        new._terms.append(dy)
    return new

  def dz(self):
    new = CRSphHarmon(self._l,self._m,False)
    new._numer = self._numer
    new._denom = self._denom
    new._terms = [] 
    for t in self._terms:
      dz = t.dz()
      if dz._n != 0:
        new._terms.append(dz)
    return new

# ======= 
#   Run
# =======
if __name__ == "__main__":

  if len(sys.argv)<2:
    print("Arguments needed.") 
    print(" 'cpp' for C++     code")
    print(" 'doc' for Doxygen code")
    print(" 'tex' for Latex   code")
    exit()
  if sys.argv[1] == "cpp":
    print("/* ==================  ")
    print(" *   Initialization    ")
    print(" * ==================*/")
    print("double Y[2*AM_MAX+1];")
    print("double dYdx[2*AM_MAX+1];")
    print("double dYdy[2*AM_MAX+1];")
    print("double dYdz[2*AM_MAX+1];")
    print("double d2Ydxdx[2*AM_MAX+1];")
    print("double d2Ydxdy[2*AM_MAX+1];")
    print("double d2Ydxdz[2*AM_MAX+1];")
    print("double d2Ydydy[2*AM_MAX+1];")
    print("double d2Ydydz[2*AM_MAX+1];")
    print("double d2Ydzdz[2*AM_MAX+1];")
    print()
    print()
  
    print("/* =============  ")
    print(" *   Harmonics    ")
    print(" * =============*/")
    print("switch (angularMomentumOfMu) {")
    for l in range(7):
      print("  case {:d}:".format(l))
      s = []
      for m in range(-l,l+1):
        s.append(CRSphHarmon(l,m))
      for m in range(-l,l+1):
        print("    Y[{:2d}] = ".format(m+l)+s[m+l].cppArray()+";")
      print("    if (_highestDerivative >= 1) {")
      for m in range(-l,l+1):
        print("      dYdx[{:2d}] = ".format(m+l)+s[m+l].dx().cppArray()+";")
        print("      dYdy[{:2d}] = ".format(m+l)+s[m+l].dy().cppArray()+";")
        print("      dYdz[{:2d}] = ".format(m+l)+s[m+l].dz().cppArray()+";")
      print("      if (_highestDerivative >= 2) {")
      for m in range(-l,l+1):
        print("        d2Ydxdx[{:2d}] = ".format(m+l)+s[m+l].dx().dx().cppArray()+";")
        print("        d2Ydxdy[{:2d}] = ".format(m+l)+s[m+l].dx().dy().cppArray()+";")
        print("        d2Ydxdz[{:2d}] = ".format(m+l)+s[m+l].dx().dz().cppArray()+";")
        print("        d2Ydydy[{:2d}] = ".format(m+l)+s[m+l].dy().dy().cppArray()+";")
        print("        d2Ydydz[{:2d}] = ".format(m+l)+s[m+l].dy().dz().cppArray()+";")
        print("        d2Ydzdz[{:2d}] = ".format(m+l)+s[m+l].dz().dz().cppArray()+";")
      print("      }")
      print("    }")
      print("    break;")
  
  
    print("  default:")
    print("    cout << \"Angular momentum too high for transformation spherical harmonics.\" << endl;")
    print("    throw exception();")
    print("    break;")
    print("}")
    print()
    print()
    print("/* =======================  ")
    print(" *   Common Finalization    ")
    print(" * =======================*/")
    print("for (unsigned int m=0;m<2*angularMomentumOfMu+1;++m){")
    print("  valueStore[m*offset] = radial*Y[m];")
    print("}")
    print("if (_highestDerivative >= 1) {")
    print("  for (unsigned int m=0;m<2*angularMomentumOfMu+1;++m){")
    print("    derivativeStore.x[m*offset] = radial*dYdx[m] + dradial*x[1]*Y[m];")
    print("    derivativeStore.y[m*offset] = radial*dYdy[m] + dradial*y[1]*Y[m];")
    print("    derivativeStore.z[m*offset] = radial*dYdz[m] + dradial*z[1]*Y[m];")
    print("  }")
    print("}")
    print("if (_highestDerivative >= 2) {")
    print("  for (unsigned int m=0;m<2*angularMomentumOfMu+1;++m){")
    print("    hessianStore.xx[m*offset] = radial*d2Ydxdx[m] + 2.0*dradial*x[1]*dYdx[m] + ddradial*x[2]*Y[m] + dradial*Y[m];")
    print("    hessianStore.xy[m*offset] = radial*d2Ydxdy[m] + dradial*x[1]*dYdy[m] + dradial*y[1]*dYdx[m] + ddradial*x[1]*y[1]*Y[m];")
    print("    hessianStore.xz[m*offset] = radial*d2Ydxdz[m] + dradial*x[1]*dYdz[m] + dradial*z[1]*dYdx[m] + ddradial*x[1]*z[1]*Y[m];")
    print("    hessianStore.yy[m*offset] = radial*d2Ydydy[m] + 2.0*dradial*y[1]*dYdy[m] + ddradial*y[2]*Y[m] + dradial*Y[m];")
    print("    hessianStore.yz[m*offset] = radial*d2Ydydz[m] + dradial*y[1]*dYdz[m] + dradial*z[1]*dYdy[m] + ddradial*y[1]*z[1]*Y[m];")
    print("    hessianStore.zz[m*offset] = radial*d2Ydzdz[m] + 2.0*dradial*z[1]*dYdz[m] + ddradial*z[2]*Y[m] + dradial*Y[m];")
    print("  }")
    print("}")
  elif sys.argv[1] == "tex":
    for l in range(7):
      print("l={:d}\\n".format(l))
      print("\\f{eqnarray*}{")
      for m in range(-l,l+1):
        s = CRSphHarmon(l,m)
        print(s.latex()+" \\\\")
      print("\\f}")
  elif sys.argv[1] == "doc":
    for l in range(7):
      print("   * l={:d}\\n".format(l))
      print("   * \\f{eqnarray*}{")
      for m in range(-l,l+1):
        s = CRSphHarmon(l,m)
        print("   * "+s.latex()+" \\\\")
      print("   * \\f}")


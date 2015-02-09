import numpy
from scipy.optimize import fmin
from  numpy.random import uniform 
import pylab

d=0.5
n=50
m=500
tau=0.3
distances=numpy.linspace(0,d,n)
f=lambda x: (sum(numpy.exp(-numpy.linspace(0,d,m)/x))-n)**2
f2 =lambda x: sum(numpy.exp(-numpy.linspace(0,d,m)/x))
f3 =lambda x: sum(numpy.exp(-numpy.linspace(0,d,m)/x))
print uniform(0.5 ,0.5)
pos_MSN_D1=numpy.array([ numpy.sqrt(uniform(-0.5 ,0.5)**2+ uniform(-0.5 , 0.5 )**2) for j in xrange(m/(numpy.pi/4))])
print pos_MSN_D1
print len(pos_MSN_D1)
f3 =lambda x: sum(numpy.exp(-pos_MSN_D1/x))
f4=lambda x: (sum(numpy.exp(-pos_MSN_D1/x))-n)**2

xopt=fmin(f4,numpy.array([0.1]))
print xopt
print f3(xopt)
print f2(xopt)
         
         
pylab.plot(numpy.linspace(0,d,m),numpy.exp(-numpy.linspace(0,d,m)/xopt))
pylab.show()

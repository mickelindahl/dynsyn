# -*- coding: cp1252 -*-

from math import *
from numpy import *
class Particle(object):

    def __init__(self,point):
	'''
	Input
	   point - list with coordinates as [x,y,z]	
	'''        
	self.point = array(point)
        
    def __str__(self):
        return str(self.point)

    def move(self, direction, distance):
	'''
	Input
	   direction - list with coordinates for direction as [x,y,z]	
	   distance - distance to travel 		
	''' 
	print array(direction)**2
	d=sqrt(sum(array(direction)**2))
	print distance*self.point/d 
	self.point=distance*self.point/d       



	# d_ = sqrt(x**2+y**2+z**2)
        #x,y,z=d*z/d_,d*y/d_,d*z/d_
        
        #self.x+=x
        #self.y+=y
        #self.z+=z
        
    def distance(self, p2):
        v=self.point-p2.point
        return sqrt(sum(v**2))
    
class Charge(Particle):
	
   def __init__(self, point, q):
      '''
	Input
	   point - list with coordinates as [x,y,z]	
               q    - electrical charge
      '''  
      super(Charge, self).__init__(point)	
      self.q=q

   def __str__(self):
        return (str(self.point)+ 'q='+
		str(self.q))

   def W(self, c2):
	''' Calculates the electrostatic energy'''
        e0=8.85*10**(-12) # Fm^-1
        r=self.distance(c2)/2.0 # radie between charges      
        return 1/(4*pi*e0)*self.q*c2.q/r 
    	
	        
if __name__ == "__main__":
    p1=Particle([1,1,1])
    p2=Particle([1,1,1])
    
    print 'p1:%s, p2:%s'%(p1, p2)
    d=2*sqrt(2)
    p1.move([-1, 0, -1], d)
    print 'p1 moved in direction -1,0,-1 distance',d,'to',p1
    print 'Distance between p1 and p2 is now', p1.distance(p2)

    c1=Charge([1,1,1],-1)
    c2=Charge([1,0,1],-1)

    print 'c1:%s, c2:%s'%(c1, c2)	
    print 'Electrostatic energy between c1 and c2', c1.W(c2)

#Körexempel
#p1:[1 1 1], p2:[1 1 1]
#[1 0 1]
#[ 2.  2.  2.]
#p1 moved in direction -1,0,-1 distance 2.82842712475 to [ 2.  2.  2.]
#Distance between p1 and p2 is now 1.73205080757
#c1:[1 1 1]q=-1, c2:[1 0 1]q=-1
#Electrostatic energy between c1 and c2 17983609388.9



    

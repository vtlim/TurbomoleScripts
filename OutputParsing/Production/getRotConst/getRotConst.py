#! /usr/bin/python

import numpy as np
import sys
import argparse

# Matt Agee 7/31/15
#
# This program is meant to read a coord file output by a Turbomole geometry 
#optimization and then  return the rotational constants associated with the 
#system specified therein.
# 
# It can be called by either supplying the coord file as an argument or by 
# running the script in a directory containing the file called 'coord'. 
# 
# It is capable of taking multiple files as argument. 


#getAtomWeight() acts as a look-up table for atomic weights given atom labels.
#
#Input: atom - The atom label [String]
#
#Output: dict[atom] - The atomic weight of the atom. [float]
def getAtomWeight( atom ):

   #atomic masses ( in 'g/mol' aka. 'u' - often written 'amu' ) taken from: 
   #<http://www.chemicalelements.com/show/mass.html> 7/31/15 by Matt Agee


   dict={'h'  : 1.00794,   'he' : 4.002602,'li' : 6.941,    'be' : 9.012182,   \
         'b'  : 10.811,    'c'  : 12.0107, 'n'  : 14.00674, 'o'  : 15.9994,    \
         'f'  : 18.9984032,'ne' : 20.1797, 'na' : 22.989770,'mg' : 24.3050,    \
         'al' : 26.981538, 'si' : 28.0855, 'p'  : 30.973761,'s'  : 32.066,     \
         'cl' : 35.4527,   'ar' : 39.948,  'k'  : 39.0983,  'ca' : 40.078,     \
         'sc' : 44.955910, 'ti' : 47.867,  'v'  : 50.9415,  'cr' : 51.9961,    \
         'mn' : 54.938049, 'fe' : 55.845,  'co' : 58.933200,'ni' : 58.9634,    \
         'cu' : 63.546,    'zn' : 65.39,   'ga' : 69.723,   'ge' : 72.61,      \
         'as' : 74.92160,  'se' : 78.96,   'br' : 79.904,   'kr' : 83.80}
   
   try:
      if not atom in dict: 
         raise KeyError(atom)

      return dict[atom]

   except KeyError:
      print "'" + atom + "' isn't in the dictionary yet. Why not add it?\n"
      raise


#makeInertiaTensor() returns the inertia tensor for a system given the
#system's atom coordinates and atom identities. 
#
#Input: xyz_at - 2D array with format [ [x1,y1,z1,at1] , [x2,y2,z2,at2], ... ].
#                Note at's should be strings identifiable by getAtomWeight.
#                [list]
#
#Output: iTens - The inertia tensor, output as a 2D array. [list]

def makeInertiaTensor( xyz_at ):
   try:
      if(not hasattr(xyz_at, '__len__') or isinstance(xyz_at, str) ):
         raise TypeError(xyz_at)

      try:
         if(len(xyz_at) == 0):
            raise TypeError(xyz_at)
      except TypeError:
         print 'Error: empty list passed to makeInertiaTensor, check your coord file'
         raise


      ixx = ixy = ixz = iyy = iyz = izz = 0

      totWeight=0
      xCOM = yCOM = zCOM = 0

      #With some algebra each inertia tensor entry can be put in a form that 
      #only requires a single loop to build.  
      for entry in xyz_at:
         if( not hasattr(entry, '__len__') or isinstance(entry, str) ):
            raise TypeError(xyz_at)

         if len(entry) != 4 or not isinstance(entry[3], str):
            raise IndexError(entry)

         atWeight = getAtomWeight(entry[3])
         totWeight += atWeight
        
         del entry[3]
         entry = map(float, entry)

         ixx += (entry[1]**2 + entry[2]**2)*atWeight
         ixy -= (entry[0]*entry[1])*atWeight
         ixz -= (entry[0]*entry[2])*atWeight

         iyy += (entry[0]**2 + entry[2]**2)*atWeight
         iyz -= (entry[1]*entry[2])*atWeight

         izz += (entry[0]**2 + entry[1]**2)*atWeight

         xCOM += entry[0]*atWeight
         yCOM += entry[1]*atWeight
         zCOM += entry[2]*atWeight


      #Note that ?COM still need to be divided by totWeight to be
      #center of mass coordinates here. 
      ixx = ixx - (yCOM**2 + zCOM**2)/totWeight
      iyy = iyy - (xCOM**2 + zCOM**2)/totWeight
      izz = izz - (xCOM**2 + yCOM**2)/totWeight

      ixy = ixy + (xCOM*yCOM)/totWeight
      ixz = ixz + (xCOM*zCOM)/totWeight
      iyz = iyz + (yCOM*zCOM)/totWeight

      iTens = [ [ixx, ixy, ixz],
                [ixy, iyy, iyz],
                [ixz, iyz, izz] ]

      return iTens

   except TypeError:
      print 'List must be 2-D\n'
      raise

   except IndexError:
      print 'Every list entry must have format [x,y,z,at]\n'
      raise



#getRotConst() calculates the 3 rorational constants for a molecule
#
#Input: coord - an opened Turbomole format coord file. 
#               Assumes Bohr radii units [file]
#
#Output: vals - A list of rotational constants in MHz sorted from largest
#               to smallest. [list]
def getRotConst( coord ):

   co = 1.804741074*10**6 #Conversion factors and coeffs to 
                          #get hbar^2/2I in units of MHz. 

   atArray = []

   started=False
   for line in coord:
      if line.find('$coord') != -1 :
         started = True

      elif line.find('$') != -1 :
         started = False

      elif started :
         line = line.strip('\n')
         entry = line.split()

         atArray.append(entry)

   iTens = makeInertiaTensor(atArray)

   vals, vects = np.linalg.eig(iTens)

   vals[0] = co/vals[0]
   vals[1] = co/vals[1]
   vals[2] = co/vals[2]

   return sorted(vals, reverse=True)



# MAIN function starts here.

if __name__ == '__main__':
   parser=argparse.ArgumentParser(description='Takes as argument Turbomole '+
   'format coord files and returns the rotational constants associated with '+
   'the molecules specified therein in MHz.')
   parser.add_argument('coords', nargs='*', default=['coord'], help='The '+
         'Turbomole format coord files you want rotation constants for. '+
         '(Default: ./coord)')
   args = parser.parse_args()

   for file in args.coords : 
      coord = open(file, 'r')

      rots = getRotConst(coord)

      print file + ' - A: ', rots[0], ' MHz , B: ', rots[1],                   \
            ' MHz , C: ', rots[2], ' MHz\n'

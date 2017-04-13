#! /usr/bin/env python

"""
Calculate the rotational constants of a molecule.

Reads Turbomole format coord files and returns the three rotational
constants of the system specified therein. Can take multiple coord
files as argument.

@author Matt Agee
"""

import numpy as np
import sys
import argparse



def get_atom_weight( atom ):
    """
    Look-up table for atomic weights given atomic labels

    Input: 
    atom - The atomic label (type: String)

    Output:
    dict[atom] - Atomic weight of atom. (type: float)
    """

    # atomic masses ( in 'g/mol' aka. 'u' - often written 'amu' ) taken from: 
    # <http://www.chemicalelements.com/show/mass.html> 7/31/15 by Matt Agee
    dict={'h'  : 1.00794,   'he' : 4.002602,'li' : 6.941,    'be' : 9.012182,
          'b'  : 10.811,    'c'  : 12.0107, 'n'  : 14.00674, 'o'  : 15.9994, 
          'f'  : 18.9984032,'ne' : 20.1797, 'na' : 22.989770,'mg' : 24.3050, 
          'al' : 26.981538, 'si' : 28.0855, 'p'  : 30.973761,'s'  : 32.066,  
          'cl' : 35.4527,   'ar' : 39.948,  'k'  : 39.0983,  'ca' : 40.078,  
          'sc' : 44.955910, 'ti' : 47.867,  'v'  : 50.9415,  'cr' : 51.9961, 
          'mn' : 54.938049, 'fe' : 55.845,  'co' : 58.933200,'ni' : 58.9634, 
          'cu' : 63.546,    'zn' : 65.39,   'ga' : 69.723,   'ge' : 72.61,   
          'as' : 74.92160,  'se' : 78.96,   'br' : 79.904,   'kr' : 83.80}
    
    try:
        return dict[atom]
    except KeyError:
        print("'" + atom + "' isn't in the dictionary yet. Why not add it?\n")
        raise


def make_inertia_tensor( xyz_at ):
    """
    Given a system's atomic coordinates and atomic labels return its inertia
    tensor

    Input:
    xyz_at - 2D list with format [ [x1,y1,z1,at1] , [x2,y2,z2,at2], ... ].
             Note at's should be strings identifiable by get_atom_weight, other
             entries are floats. (type: N x 4 iterable)

    Output:
        itens - The inertia tensor, output as a 2D array. (type: list)
    """
    # Check iterable
    try:
        assert(hasattr(xyz_at,'__iter__'))
    except AssertionError:
        print('Error: Argument to make_inertia_tensor must be iterable')
        raise
    
    # Check argument has elements
    try:
        assert(len(xyz_at) != 0)
    except TypeError: 
        print('Error: Iterable passed to make_inertia_tensor'
            + ' must implement __len__')
        raise
    except AssertionError:
        print('Error: Empty list passed to make_inertia_tensor, check'
            + ' your coord file')
        raise

    ixx = ixy = ixz = iyy = iyz = izz = 0
    tot_weight = 0
    xcom = ycom = zcom = 0
    # With some algebra each inertia tensor entry can be put in a form
    # requiring a single loop to build.  
    for entry in xyz_at:
        try:
            assert(len(entry) == 4)
        except TypeError: 
            print('Error: Entries of xyz_at must implement __len__')
            raise
        except AssertionError:
            print('Error: Entries of xyz_at must have 4 entries')
            raise

        try:
            entry[0:3] = list(map(float, entry[0:3]))
            assert(isinstance(entry[3],str))
        except ValueError:
            print('Error: Entries of xyz_at must have floats as their first'
                + ' 3 entries')
            raise
        except TypeError:
            print('Error: Entries of xyz_at must be indexable')
            raise
        except AssertionError:
            print('Error: Entries of xyz_at must have a string as their'
                + ' fourth entry')
            raise

        at_weight = get_atom_weight(entry[3])
        tot_weight += at_weight
      
        ixx += (entry[1]**2 + entry[2]**2)*at_weight
        ixy -= (entry[0]*entry[1])*at_weight
        ixz -= (entry[0]*entry[2])*at_weight

        iyy += (entry[0]**2 + entry[2]**2)*at_weight
        iyz -= (entry[1]*entry[2])*at_weight

        izz += (entry[0]**2 + entry[1]**2)*at_weight

        xcom += entry[0]*at_weight
        ycom += entry[1]*at_weight
        zcom += entry[2]*at_weight

    # Note ?com must be divided by tot_weight to be
    # center of mass coordinates here. 
    ixx = ixx - (ycom**2 + zcom**2)/tot_weight
    iyy = iyy - (xcom**2 + zcom**2)/tot_weight
    izz = izz - (xcom**2 + ycom**2)/tot_weight

    ixy = ixy + (xcom*ycom)/tot_weight
    ixz = ixz + (xcom*zcom)/tot_weight
    iyz = iyz + (ycom*zcom)/tot_weight

    itens = [ [ixx, ixy, ixz],
              [ixy, iyy, iyz],
              [ixz, iyz, izz] ]

    return itens


def get_rot_const( coord ):
    """
    Calculates the 3 rotational constants for a molecule

    Input:
    coord - Opened Turbomole format coord file. 
            Assumes units of Bohr radii (file)

    Output:
    vals - List of rotational constants in MHz, sorted 
           largest to smallest (type: list)
    """
    # Conversion factors and coeffs to get hbar^2/2I in units of MHz. 
    co = 1.804741074*10**6 
    at_array = []
    started=False
    for line in coord:
        if line.find('$coord') != -1 :
            started = True

        elif line.find('$') != -1 :
            started = False

        elif started :
            line = line.strip('\n')
            entry = line.split()

            at_array.append(entry)

    itens = make_inertia_tensor(at_array)
    vals, vects = np.linalg.eig(itens)

    vals[0] = co/vals[0]
    vals[1] = co/vals[1]
    vals[2] = co/vals[2]

    return sorted(vals, reverse=True)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description='Takes as argument a list of'
            + ' Turbomole format coord files and returns the rotational'
            + ' constants associated with the molecules specified therein.'
            + ' (in MHz)')
    parser.add_argument('coords', nargs='*', default=['coord'], help='A'
            + ' list of paths to Turbomole format coord files.'
            + ' (Default: ./coord)')
    args = parser.parse_args()

    for fil in args.coords: 
        try:
            coord = open(fil, 'r')
            rots = get_rot_const(coord)
            print(fil + ' - A: ', rots[0], ' MHz, B: ', rots[1], 
                  ' MHz, C: ', rots[2], ' MHz')
        except IOError:
            print('Error: Could not open file '+fil+' for reading.'
                + ' Processing remaining files')

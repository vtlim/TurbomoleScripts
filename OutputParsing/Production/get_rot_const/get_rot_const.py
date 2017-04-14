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


class KnownError(Exception):
    """Used to report anticipated errors"""
    pass


def get_atom_weight(atom):
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
        raise KnownError("'" + atom + "' isn't in the dictionary yet. Why not"
                + " add it?")


def make_inertia_tensor(xyz_at):
    """
    Given a system's atomic coordinates and atomic labels return its inertia
    tensor. 

    Input:
    xyz_at - 2D list with format [ [x1,y1,z1,at1] , [x2,y2,z2,at2], ... ].
             Note at's should be strings identifiable by get_atom_weight, other
             entries are floats. (type: N x 4 iterable)

    Output:
    itens - The inertia tensor, output as a 2D array. (type: list)
    """
    # Check iterable
    assert hasattr(xyz_at,'__iter__'), \
        'Argument to make_inertia_tensor must be iterable'
    
    # Check xyz_at has elements
    try:
        assert len(xyz_at) != 0, \
            'Empty list passed to make_inertia_tensor, check coord file'
    except TypeError: 
        raise KnownError('Iterable passed to make_inertia_tensor must'
                + ' implement __len__')

    ixx = ixy = ixz = iyy = iyz = izz = 0
    tot_weight = 0
    com = np.array([0, 0, 0],dtype=np.float64)
    # With some algebra each inertia tensor entry can be put in a form
    # requiring a single loop to build.  
    for entry in xyz_at:
        try:
            assert len(entry) == 4,'Entries of xyz_at must have 4 entries'
        except TypeError: 
            raise KnownError('Entries of xyz_at must implement __len__')

        try:
            assert isinstance(entry[3],str), \
                'Entries of xyz_at must have a string as their fourth entry'
            at_weight = get_atom_weight(entry[3])
            entry = np.array(list(map(float, entry[0:3])))
        except ValueError:
            raise KnownError('Entries of xyz_at must have floats as their'
                    + ' first 3 entries')
        except TypeError:
            raise KnownError('Entries of xyz_at must be indexable')

        tot_weight += at_weight

        ixx += (entry[1]**2 + entry[2]**2)*at_weight
        iyy += (entry[0]**2 + entry[2]**2)*at_weight
        izz += (entry[0]**2 + entry[1]**2)*at_weight

        ixy -= (entry[0]*entry[1])*at_weight
        ixz -= (entry[0]*entry[2])*at_weight
        iyz -= (entry[1]*entry[2])*at_weight

        com += entry*at_weight

    # Note ?com must be divided by tot_weight to be
    # center of mass coordinates here. 
    ixx -= (com[1]**2 + com[2]**2)/tot_weight
    iyy -= (com[0]**2 + com[2]**2)/tot_weight
    izz -= (com[0]**2 + com[1]**2)/tot_weight

    ixy += (com[0]*com[1])/tot_weight
    ixz += (com[0]*com[2])/tot_weight
    iyz += (com[1]*com[2])/tot_weight

    itens = [ [ixx, ixy, ixz],
              [ixy, iyy, iyz],
              [ixz, iyz, izz] ]
    return itens


def read_coord(fil):
    """
    Read the entries in a Turbomole format coord file. 

    Input:
    fil - String with path to Turbomole format coord file. (Type: String)

    Output:
    at_array - List of coordinates (float) and atomic labels (String) in 
               coord file. Each entry has format [x, y, z, at].
               (Type: N x 4 list) 
    """
    try:
        coord = open(fil, 'r')
    except IOError:
        raise KnownError('Could not open file '+fil+' for reading.')
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
            assert len(entry) == 4, \
                'The line for each atom in coord file must have 4 entries'

            try:
                entry[0:3] = list(map(float, entry[0:3]))
            except ValueError:
                raise KnownError('First 3 entries for each atom must be'
                        + ' floats')

            at_array.append(entry)
    return at_array


def get_rot_const(fil):
    """
    Calculates the 3 rotational constants for a molecule. (in MHz)

    Formula: hbar^2/(2I), where I is the moment of inertia about
    a principle axis of rotation

    Input:
    fil - String with path to Turbomole format coord file. 
          Assumes units of Bohr radii. (Type: String)

    Output:
    vals - List of rotational constants in MHz, sorted 
           largest to smallest (type: list)
    """
    # Conversion factors and coeffs to get hbar^2/2I in units of MHz. 
    co = 1.804741074*10**6 
    tol = 10**(-8)

    at_array = read_coord(fil)
    itens = make_inertia_tensor(at_array)
    vals, vects = np.linalg.eig(itens)

    # Atoms and linear molecules have less than 3 moments of inertia. 
    # The excess moments of inertia have 0-valued eigenvalues
    for i in range(len(vals)):
        if vals[i] < tol:
            vals[i] = -1

    vals = co/vals

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
        # Assume error local to a single file
        try:
            rots = get_rot_const(fil)
            axes = ['A: ','B: ','C: ']
            print(fil + ' - ',end='')
            for i in range(len(axes)):
                if rots[i] > 0:
                    print(axes[i], rots[i], ' MHz ',end='')
            print()
        except (AssertionError, KnownError) as err:
            print('Error:',err,file=sys.stderr)
            print('get_rot_const failed for file: '+ fil + '\n'
                + 'Processing remaining files',file=sys.stderr)

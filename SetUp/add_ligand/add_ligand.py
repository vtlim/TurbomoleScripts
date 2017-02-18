#! /usr/bin/python

"""
This script combines two molecules together. The user specifies the
point they want the molecules to be joined at by an x in the 
coordinate file (One in each molecule's file). The user specifies
the distance apart these two x atoms must be. The script then joins
the two fragments in a physically reasonable way, treating each molecule
as a rigid body. 

The current algorithm for determining optimal rotation is very simplistic.
It could be improved by interfacing this script with a more advanced 
molecular mechanics package, or just by implementing it here. 

Points of improvement:
    - Determine convex hull only for molecules close to connection site
    - Get the atoms contributing to the solvent available surface instead
      of the convex hull (hull might be a good start point though)
    - Better implement force field. This could be interfacing with an existing
      package or doing it right here. Don't overcomplicate things though, this 
      script is supposed to give a fit that is just "close enough"
    - Implement collision detection for molecules. If the first guess is bad
      there could be slow optimization and the final minima might not be 
      very good. 
    - Fix for host-guest complexes. I expect this script will fail badly when 
      putting one molecule into another. Instant improvement would be seen
      by going to solvent available surface from convex hull. 

@author Matt Agee
"""

from __future__ import division, print_function

import numpy as np
import math
import elements
import argparse
import collections

from scipy.spatial import ConvexHull
from scipy.special import expit
from scipy.optimize import minimize


try:
    from functools import lru_cache
except ImportError:
    print('Use Python 3 for memoization support')
    pass


class Fragment(object):
    def __init__(self,coord_file,f_type):
        self.precision = 1e-5

        # store coordinates and labels in array
        self.file_name = coord_file
        coords = open(self.file_name,'r')
        self.read(coords,f_type)

        # Find the point other fragments will connect to this one
        print(self.at_labels)
        connect = np.where(self.at_labels == 'X')
        try:
            assert len(connect[0]) == 1
        except AssertionError:
            print('Error: Fragment file ' + coord_file + ' must have ONE '
                 + "entry with label 'x' specifying the connection point")
            raise
        print(connect[0],len(connect))
        self.connect = self.coord[connect[0][0]]
        self.coord = np.delete(self.coord,connect[0][0],axis=0)
        self.at_labels = np.delete(self.at_labels,connect[0][0],axis=0)


        # Need to have some atoms
        try:
            assert len(self.coord) != 0
        except AssertionError:
            print('Error: No coordinates found in '+coord_file)
            raise

        # Orient molecule in a convenient way
        self.orient_frag()

        # Only convex hull will be used when optimizing
        # to final orientation
        self.hull, self.hull_ats = self.build_hull()


    def __len__(self):
        """Returns number of elements in convex hull of molecule"""
        return len(self.hull)


    def __str__(self):
        """
        Returns name of the file this fragment got its coordinates from
        """
        return self.file_name

    
    def read(self,coords,f_type):
        try:
            if f_type == 'Turbomole':
                self.coord, self.at_labels = self.read_coord(coords)
            elif f_type == 'xyz':
                self.coord, self.at_labels = self.read_xyz(coords)
            else:
                raise NotImplemented
        except NotImplemented:
            print('Error: File ' + self.file_name + ' of file type ' + f_type 
                 +' could not be read. That file type not yet readable by the '
                 +'Fragment class.')
            raise

    def read_coord(self,coords):
        """
        Read a Turbomole format coord file, storing the atomic coordinates 
        and atomic labels in a numpy array.
    
        Input: coords - opened Turbomole format coord file
        
        Output: coord_array - numpy array of atomic coordinates (atoms x 3)
                at_array - numpy array of atomic coordinates (atoms x 3)
        """
        bohr_2_ang = 0.52917724900001 # We want our units in Angstroms
        # We want only the atomic label, this will be at most two characters
        # long. 
        at_array = np.array([],dtype='S5')
        coord_array = np.array([])
    
        started = False
        for line in coords:
            if line.find('$coord') != -1:
                started = True
            elif started and line[0] == '$':
                break
            elif started :
                line = line.strip('\n')
                coord_array = np.append(coord_array,
                                 bohr_2_ang
                                *np.array(map(np.float64,line.split()[0:3])),
                                 axis=0)
                at_array = np.append(at_array,line.split()[3].capitalize())
    
        coord_array = np.reshape(coord_array,(len(coord_array)//3,3))
        return coord_array, at_array
    
    
    def read_xyz(self,coords):
        """
        Read a .xyz format coord file, storing the atomic coordinates 
        and atomic labels in a numpy array.
    
        Input: coords - opened .xyz format coord file
        
        Output: coord_array - numpy array of atomic coordinates (atoms x 3)
                at_array - numpy array of atomic labels (atoms)
        """
        # We want only the atomic label, this will be at most two characters
        # long. Could read atom count from xyz and allocate arrays at
        # start for performance, this more useable though
        at_array = np.array([],dtype='S5')
        coord_array = np.array([])
    
        count = 0
        for line in coords:
            # skips the file header
            if count > 1:
                line = line.strip('\n')
                coord_array = np.append(coord_array,
                                np.array(map(np.float64,line.split()[1:4])),
                                axis=0)
                at_array = np.append(at_array,line.split()[0].capitalize())
            count += 1

        coord_array = np.reshape(coord_array,(len(coord_array)//3,3))
    
        return coord_array, at_array

    
    def write(self,outfile,out_t,header=True):
        try:
            if out_t == 'xyz':
                self._write_xyz(outfile,header)
            elif out_t == 'Turbomole':
                self._write_coord(outfile,header)
            else:
                raise NotImplemented
        except NotImplemented:
            print('Error: File type ' + out_t + ' is not yet implemented '
                 +'for writing.')
            raise


    def _write_coord(self,outfile,header=True):
        """
        outfile should be an opened file for writing
        """
        ang_2_bohr = 1.889725989 # Need to convert back to bohr

        if header:
            outfile.write('$coord\n')
        for coords,label in zip(self.coord,self.at_labels):
            outfile.write(('{:20.14f}  '*3 + '{:>5s}\n').format(
                          coords[0],coords[1],coords[2],label.lower()))
        if header:
            outfile.write('$end')


    def _write_xyz(self,outfile,header=True):
        """
        outfile should be an opened file for writing
        """
        n = len(self.coord)
        if header:
            outfile.write(str(n)+'\n\n')

        count = 1
        for coords,label in zip(self.coord,self.at_labels):
            outfile.write(('{:4s}'+'{:15.9f} '*3).format(
                          label.upper(),coords[0],coords[1],coords[2]))
            if count != n:
                outfile.write('\n')
            count += 1


    # This routine based on one found in: 
    # http://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector
    def rotate_frag(self,rots,hull=False):
        """
        Return the fragment coordinates rotated counterclockwise about 
        a rotation axis, with both degree rotated and axis direction 
        specified by the Euler-Rodrigues formula vector parameters in rots. 
        """
        aa = 1 - np.dot(rots,rots)
        a = math.sqrt(aa)
        bb, cc, dd = rots[0]*rots[0], rots[1]*rots[1], rots[2]*rots[2]
        bc, ad, ac, ab, bd, cd = rots[0]*rots[1], a*rots[2], a*rots[1], \
                                 a*rots[0], rots[0]*rots[2], rots[1]*rots[2]

        r_mat = np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                          [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                          [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])
        if hull:
            return np.dot(self.hull,r_mat.T)
        else:
            return np.dot(self.coord,r_mat.T), np.dot(self.connect,r_mat.T)

    
    def rotate(self,rots):
        """Sets the fragment coordinates to newly rotated values"""
        self.coord, self.connect = self.rotate_frag(rots)


    def shift_origin(self,new_origin):
        """Shifts the fragment coordinates to new origin supplied"""
        self.coord = np.add(self.coord,-new_origin)
        self.connect = np.add(self.connect,-new_origin)


    # We orient the fragment to simplify the math later when docking
    # and (I expect) to get a better starting guess
    def orient_frag(self):
        """
        Shifts the origin to the center of volume and rotates the fragment
        about its center of volume so the connection point is aligned with 
        the z-axis
        """
        cov = np.mean(self.coord,axis=0)
        self.shift_origin(cov)

        z = np.array([0,0,1])
        con_mag = np.dot(self.connect,self.connect)
        if con_mag < self.precision:
            return

        # note math.sqrt is faster for scalars than np.sqrt
        # We use the half-angle formula to get sin(\phi/2)
        sin_rot = math.sqrt((1-np.dot(self.connect,z)/con_mag)/2)
        rots = np.cross(z,self.connect)/con_mag
        rots *= sin_rot
        self.rotate(rots)
        dr_fil = open('rotaft.xyz','w')


    # For now we build the entire hull. If performance matters we 
    # can improve this routine by building only partial hulls
    def build_hull(self):
        """
        Method for finding points in the convex hull we want to compare with
    
        Input: coord_array - Cartesian coordinates for molecule (atoms x 3)
               at_array - numpy array of atomic labels (atoms)
        
        Output: hull - The atomic coordinates making up the convex hull 
                       of our molecule
                hull_ats - The atomic labels of the hull atoms
        """
        # hull needs to be built once, efficiency not important
        hull = self.coord[ConvexHull(self.coord).vertices]
        hull_ats = []
        for at in hull:
            hull_ats.append(self.at_labels[np.where(self.coord == at)[0][0]])
    
        hull_ats = np.array(hull_ats,dtype='S5')
        return hull, hull_ats


    def align_frags(self,other_frag,dist):
        """
        Places the fragment with the smaller number of hull
        points on top of the fragment with more. The new origin
        is the connection point of the fragment with the most 
        members in its hull. Molecules are brought into alignment 
        in such a way that the old origins are as far apart as possible. 
        """
        try:
            if(not isinstance(other_frag, Fragment)):
                raise TypeError
        except TypeError:
            print('ERROR: '+str(self)+' attempted to align with something '
                 + 'not of the Fragment class.')
            raise

        if len(self) <= len(other_frag):
            rot = self
            anchor = other_frag
        else:
            rot = other_frag
            anchor = self

        
        # Might be able to come up with better axis of rotation than
        # just y-axis. Goal is to get rest of molecule out of the way
        flip = np.array([1,0,0])
        rot.shift_origin(rot.connect)
#        dr_fil = open('flipb4.xyz','w')
#        rot.write(dr_fil,'xyz')
        rot.rotate(flip)
#        dr_fil = open('flipaft.xyz','w')
#        rot.write(dr_fil,'xyz')


#        dr_fil = open('alignb4.xyz','w')
#        dr_fil.write(str(len(rot.coord)+len(anchor.coord)) + '\n\n')
#        self.write(dr_fil,'xyz',header=False)
#        dr_fil.write('\n')
#        anchor.write(dr_fil,'xyz',header=False)

        # The new origin is dist below this fragments
        # connection, where the other fragments connect 
        # should be
#        dr_fil = open('origrotb4.xyz','w')
#        rot.write(dr_fil,'xyz')
#        dr_fil = open('origanchb4.xyz','w')
#        rot.write(dr_fil,'xyz')

        rot.shift_origin(np.array([0,0,-dist]))
        anchor.shift_origin(anchor.connect)

#        dr_fil = open('origrotaft.xyz','w')
#        rot.write(dr_fil,'xyz')
#        dr_fil = open('origanchaft.xyz','w')
#        rot.write(dr_fil,'xyz')


#        dr_fil = open('alignaft.xyz','w')
#        dr_fil.write(str(len(rot.coord)+len(anchor.coord)) + '\n\n')
#        self.write(dr_fil,'xyz',header=False)
#        dr_fil.write('\n')
#        anchor.write(dr_fil,'xyz',header=False)


            

class ForceField(object):
    def __init__(self,force_field):
        """
        Sets the force field related functions to the requested 
        force field type.
        """
        try:
            if force_field == 'MMFF94':
                self.name = force_field
                self.sep = self._sep_MMFF94
                self.depth = self._depth_MMFF94
                self.inter = self._inter_MMFF94
                self._init_MMFF94()
            else:
                raise NotImplemented
        except NotImplemented:
            print('Error: the force field ' + force_field + ' is not yet '
                 +'implemented')
            raise

    def __str__(self):
        return self.name

    
    def _init_MMFF94(self):
        """
        These Van der Waal parameters pulled from the MMFFVDW.PAR
        file associated with MMFF94. The original file has much more
        nuance. I take one instance of each atom type and say that 
        defines the atom. 
        """
        atom = collections.namedtuple('atom',
        'polar, eff_ele, sep_fact, dep_fact' )
        self.at_dict = {}
        self.at_dict['H']  = atom(0.250, 0.800, 4.200, 1.209)
        self.at_dict['Li'] = atom(0.15,  2, 4, 1.3)
        self.at_dict['C']  = atom(1.050, 2.490, 3.890, 1.282)
        self.at_dict['N']  = atom(.15, 2.820, 3.890, 1.282)
        self.at_dict['O']  = atom(0.70, 3.150, 3.890, 1.282)
        self.at_dict['F']  = atom(0.35, 3.480, 3.890, 1.282)
        self.at_dict['Na'] = atom(0.4, 3.5, 4, 1.3)
        self.at_dict['Mg'] = atom(0.35, 3.5, 4, 1.3)
        self.at_dict['P']  = atom(3.600, 4.500, 3.320, 1.345)
        self.at_dict['S']  = atom(3.00, 4.800, 3.320, 1.345)
        self.at_dict['Cl'] = atom(2.300, 5.100, 3.320, 1.345)
        self.at_dict['K']  = atom(1.0, 5, 4, 1.3)
        self.at_dict['Ca'] = atom(0.9, 5, 4, 1.4)
        self.at_dict['Fe'] = atom(0.45, 6, 4, 1.4)
        self.at_dict['Cu'] = atom(0.35, 6, 4, 1.4)
        self.at_dict['Zn'] = atom(0.43, 6, 4, 1.4)
        self.at_dict['Br'] = atom(3.400, 6.000, 3.190, 1.359)
        self.at_dict['I']  = atom(5.500, 6.950, 3.080, 1.404)



    # Could probably just save this info for all possible pairs if
    # we wanted to, then no memoization needed
    # Decorator only defined in Python 3
    #@lru_cache
    def _sep_MMFF94(self,label):
        """
        Calculates the minimum energy separation of the two atoms as 
        defined in the Merck Molecular Force Field II (MMFF94). There 
        are some environment specific modifications that should be made 
        (polar hydrogen, donor-accpetor pair) that aren't. Revisit later.
        DOI:10.1002/(SICI)1096-987X(199604)17:5/6<520::AID-JCC2>3.0.CO;2-W
        """
        try:
            r1 = (self.at_dict[label[0]].sep_fact
                 *self.at_dict[label[0]].polar)
        except KeyError:
            r1 = (self.at_dict['Fe'].sep_fact
                 *self.at_dict['Fe'].polar)
        if label[0] == label[1]:
            return r1

        try:
            r2 = (self.at_dict[label[1]].sep_fact
                 *self.at_dict[label[1]].polar)
        except KeyError:
            r2 = (self.at_dict['Fe'].sep_fact
                 *self.at_dict['Fe'].polar)

        b = .2; beta = 12
        gam12 = (r1 - r2)/(r1 + r2)
        sep = .5*(r1 + r2)*(1 + b*(1 - math.exp(-beta*gam12**2)))

        return sep
    
    
    # Decorator only defined in Python 3
    #@lru_cache
    def _depth_MMFF94(self,label,sep):
        """
        Calculates the energy well-depth as defined in the Merck 
        Molecular Force Field II (MMFF94). There are some 
        environment specific modifications that should be made (polar 
        hydrogen, donor-accpetor pair) that aren't. Revisit later.
        DOI:10.1002/(SICI)1096-987X(199604)17:5/6<520::AID-JCC2>3.0.CO;2-W
        """
        try:
            g1 = self.at_dict[label[0]].dep_fact
            p1 = self.at_dict[label[0]].polar
            n1 = self.at_dict[label[0]].eff_ele
        except KeyError:
            # Every element not defined in some way defaults to 
            # iron because they would likely be larger elements
            g1 = self.at_dict['Fe'].dep_fact
            p1 = self.at_dict['Fe'].polar
            n1 = self.at_dict['Fe'].eff_ele


        try:
            g2 = self.at_dict[label[1]].dep_fact
            p2 = self.at_dict[label[1]].polar
            n2 = self.at_dict[label[1]].eff_ele
        except KeyError:
            # Every element not defined in some way defaults to 
            # iron because they would likely be larger elements
            g2 = self.at_dict['Fe'].dep_fact
            p2 = self.at_dict['Fe'].polar
            n2 = self.at_dict['Fe'].eff_ele
            pass
    
        depth = (181.16*g1*g2*p1*p2
              / (math.sqrt(p1/n1) + math.sqrt(p2/n2)
              * sep**6))

        return depth
    
    
    def _inter_MMFF94(self,dist,label):
        """
        Returns the interatomic interaction energy. Uses Slater's rules to 
        calculate the partial charges
        """
        # Van der Waal interaction energy
        sep = self._sep_MMFF94(label)
        w_depth = self._depth_MMFF94(label,sep)
        e_vdw = ( w_depth
                * (1.07*sep/(dist + .07*sep))**7
                * (1.12*sep**7/(dist**7 + .12*sep**7) - 2))

        ## electrostatic interaction energy
        ## Neglect for now. Need routine to calculate formal charge
        #d = .05 # units in Angstroms
        #q1, q2 = self.partial[label[0]], self.partial[label[1]]
        #e_inter = 332.0716*q1*q2/(dist + d)
        return e_vdw




# May want to rewrite to take VDW radii into account with cut
def calc_dists(s_frag,f_frag,cut=3.0):
    """
    Calculate non-negligible dists and return them and with a tuple
    storing the atom labels of the corresponding elements.
    """
    # We allocate the max number of elements possible 
    # so we don't need to reallocate 
    dists = np.empty([len(s_frag)*len(f_frag)],dtype=np.float64)
    labels = np.empty([len(s_frag)*len(f_frag)],dtype=('S5',2))
    count = 0
    cut2 = cut**2
    for i in range(len(s_frag)):
        i_hull = s_frag.hull[i]
        for j in range(len(f_frag)):
            dist = np.dot(np.add(i_hull, - f_frag.hull[j]),
                          np.add(i_hull, - f_frag.hull[j]))
            if dist > cut2:
                dists[count] = math.sqrt(dist)
                labels[count] = (s_frag.hull_ats[i],f_frag.hull_ats[j])
                count += 1

    # Want to use length later so shrink, could pass back 
    # count too I suppose. 
    dists = np.resize(dists,count)
    labels = np.resize(labels,(count,2))
    return dists, labels


def score_rotation(rots,s_frag,f_frag,f_field):
    """
    Calculates the value of the objective function (energy) for this
    optimization. This objective function uses stored partial charge
    and VDW radii to calculate an energy value. 
    """
    # I use unconstrained variables and map them to [0,1] 
    # with expit (the logistic function) when optimizing. The hope 
    # is that this will be faster. Need to verify.
    norm_rots = expit(rots)
    s_frag.guess_hull = s_frag.rotate_frag(norm_rots,hull=True)

    dists, labels = calc_dists(s_frag,f_frag)

    score = 0
    for i in range(len(dists)):
        score += f_field.inter(dists[i],labels[i])

    print(score)
    return score


def determine_file_type(file_name):
    try:
        if file_name.find('coord') != -1:
            file_t = 'Turbomole'
        elif file_name.find('.xyz') != -1:
            file_t = 'xyz'
        else:
            raise TypeError
    except TypeError:
        print('Error: could not determine file type of '
             + file_name + '\nPlease specify.')
        raise

    return file_t


def add_ligand(file1,file2,dist,file_t1='',file_t2='',out='combo',out_t='xyz'):
    if file_t1 == '':
        file_t1 = determine_file_type(file1)
    if file_t2 == '':
        file_t2 = determine_file_type(file2)

    frag1 = Fragment(file1,file_t1)
    frag2 = Fragment(file2,file_t2)
    if len(frag1) <= len(frag2):
        s_frag, f_frag = frag1, frag2
    else:
        s_frag, f_frag = frag2, frag1
    f_frag.align_frags(s_frag, dist)


    f_field = ForceField('MMFF94')

    # Initial guess corresponds to vector Euler-Rodrigues coefficients of 
    # roughly 0, ie. no rotation
    try:
        opt = minimize(score_rotation,[0,0,0], #[-1000,-1000,-1000],
                       (s_frag,f_frag,f_field),callback=print)
        if not opt.success:
            raise RuntimeError
    except RuntimeError:
        print('Error: Optimization of ligand alignment failed. '
             +"I'm sorry.")
        raise
    print(expit(opt.x))
    print('start',s_frag.coord)
    s_frag.rotate(expit(opt.x))
    print('aft',s_frag.coord)

    try:
        if out_t == 'xyz':
            if out.find('.xyz') == -1:
                out += '.xyz'
            outfile = open(out,'w')
            outfile.write(str(len(s_frag.coord)+len(f_frag.coord)) + '\n\n')
        elif out_t == 'Turbomole':
            if out.find('coord') == -1:
                out += '.coord'
            outfile = open(out,'w')
            outfile.write('$coord\n')
        else:
            raise TypeError
    except TypeError:
        print('Error: Output file type ' + out_t + ' not implemented yet.')
        raise

    s_frag.write(outfile,out_t,header=False)
    if out_t == 'xyz':
        outfile.write('\n')
    f_frag.write(outfile,out_t,header=False)
    if out_t == 'Turbomole':
        outfile.write('$end')

    print('Molecules joined succesfully')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Tool for combining two '
        + 'molecular fragments together. The user supplies two coordinate '
        + 'files defining the molecules. These coordinates have one entry '
        + 'with label x. This marks the connection point at which the '
        + 'molecules are to be joined. The user specifies the distance '
        + 'these two x should be apart and the program attempts to find '
        + 'a physically reasonable geometry satisfying this.')

    parser.add_argument('file1',nargs=1,
            help='File defining first molecular fragment to be joined')
    parser.add_argument('file2',nargs=1,
            help='File defining second molecular fragment to be joined')
    parser.add_argument('dist',nargs='?',default=1.5,
            help='Distance apart (in Angstroms) the two connection '
                +'points should be (Default: 1.5)')

    parser.add_argument('-t1','--type1', default='',
            help='Format for file1. Possibilities are "xyz" and '
                +'"Turbomole" (Default: automatically determined)')
    parser.add_argument('-t2','--type2', default='',
            help='Format for file2. Possibilities are "xyz" and '
                +'"Turbomole" (Default: automatically determined)')
    parser.add_argument('-o','--outfile', default='combo',
            help='Filename for output molecule (Default: combo)')
    parser.add_argument('-to','--type_out', default='xyz',
            help='Format for output. Possibilities are "xyz" and '
                +'"Turbomole" (Default: xyz)')

    args = parser.parse_args()

    
    add_ligand(args.file1[0], args.file2[0], args.dist, args.type1,
               args.type2, args.outfile, args.type_out)

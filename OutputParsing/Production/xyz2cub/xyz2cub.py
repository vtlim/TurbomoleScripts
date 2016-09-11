#! /usr/bin/python

from subprocess           import call, Popen, PIPE
import re
import sys
import argparse

#xyz2cub - This program takes as input a Turbomole xyz file (with grid 
#          point coordinates followed by property values) and converts 
#          it to a cub format file. This xyz file is NOT the 
#          same as a coordinate xyz files. 
#
#          It does this by looking at the number of entries and the 
#          distance apart they are. It relies on the fact values start
#          at the gridpoint with the smallest possible coordinates and 
#          the file iterates along x, y, z in that order. Each y/z
#          increment must be marked by a newline. The output file will 
#          be the input xyz file's basename followed by .cub
def xyz2cub (xyzname,coord):
    xyzname=xyzname[0]
    xyz_end=re.compile('\.[xX][yY][zZ]$')

    #coord file must be in the current directory
    coordfile=open(coord,'r')


    #Requires the file name to end in xyz, asks for name again if
    #not. 
#    if not xyz_end.search(xyzname) :
#        bad_input=True
#        while bad_input :
#            xyzname=raw_input("File name must end in '.xyz'"+ \
#                    " please re-enter (q to quit): ")
#            if xyzname == 'q':
#                sys.exit(1)
#            elif xyz_end.search(xyzname) :
#                bad_input=False

    #I don't put .cub in the sub field in case file
    #does not end in .xyz
    cubname=str(re.sub(xyz_end,'',xyzname))+'.cub'

    #opening the two files to prepare for conversion
    cubfile=open(cubname,'w')
    xyzfile=open(xyzname,'r')

    atNum = writeAtCoords(coordfile,cubfile)
    printCubVals(xyzfile,cubfile,atNum)




#getAtomNumber - Given a TURBOMOLE atom label, looks up the atomic number
#    input: atom - A string representing the atom label
#
#    return: dict[atom] - The atomic number associated with the atom (int)

def getAtomNumber( atom ):

   dict={'h'  : 1,  'he' : 2,  'li' : 3,  'be' : 4,  'b'  : 5,  'c'  : 6,      \
         'n'  : 7,  'o'  : 8,  'f'  : 9,  'ne' : 10, 'na' : 11, 'mg' : 12,     \
         'al' : 13, 'si' : 14, 'p'  : 15, 's'  : 16, 'cl' : 17, 'ar' : 18,     \
         'k'  : 19, 'ca' : 20, 'sc' : 21, 'ti' : 22, 'v'  : 23, 'cr' : 24,     \
         'mn' : 25, 'fe' : 26, 'co' : 27, 'ni' : 28, 'cu' : 29, 'zn' : 30,     \
         'ga' : 31, 'ge' : 32, 'as' : 33, 'se' : 34, 'br' : 35, 'kr' : 36 }
   
   try:
      if not atom in dict: 
         raise KeyError(atom)

      return dict[atom]

   except KeyError:
      print "'" + atom + "' isn't in the dictionary yet. Why not add it?\n"
      raise



#isClose - An equality check for floating point numbers that
#          may fail to be equal when they should be due to rounding 
#          errors.
#
#   input:   a,b - floating point numbers that need to be compared
#        rel_tol - The max percent a,b are allowed to differ and 
#                  still be considered equal. 
#        abs_tol - The absolute tolerance for a,b to be considered
#                  equal.

def isClose(a, b, rel_tol=1e-09, abs_tol=0.0):
        return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


#writeAtCoords -  searches the supplied file 'coord' for a section
#                 beginning with '$coord' and ending with a line 
#                 containing '$'.
#                 Within this section it goes line by line taking 
#                 the atomic label, determining the atomic number 
#                 and returning it and the coordinates formatted for
#                 a cube file. 
#
#   input:  coord - The opened coord file for reading
#           cube  - The opened cube file for writing
#
#   return: count - The total number of atoms in the coord file

def writeAtCoords(coord,cube):
    started=False
    count=0

    for line in coord:
        if( line.find('$coord') != -1 ):
            started=True
        elif( started ):
            if( line.find('$') != -1 ):
                started=False
            else :
                entry=line.split()
                atLabel=entry[3]
                loc=entry[0:3]
                atNum=getAtomNumber(atLabel)
                if count == 0 :
                    cube.write('%5s %12s %12s %12s %12s' % (atNum, '0.000000', loc[0], loc[1], loc[2]))
                else :
                    cube.write('\n%5s %12s %12s %12s %12s' % (atNum, '0.000000', loc[0], loc[1], loc[2]))
                count+=1

    return count


#printCubVals - Pulls the potential values from the xyz input file and 
#               appends them to the cube file. While doing this it determines
#               the number of points along each vector, the increment in each 
#               direction, and the origin for the cube file. It then writes 
#               these values with a header at the top of the file. Note that
#               I use 1 to refer to the fastest moving vector, 2 for the 
#               second fastest and 3 for the slowest. These are typically 
#               x, y and z respectively in the xyz file. 
#
#   input:  xyz - The opened xyz file to be read
#           cub - The opened cube file for writing to
#         atNum - The number of atoms in the associated coord file (int)

def printCubVals(xyz,cub,atNum):
    cub.seek(0,2) #go to end of file
    comment=re.compile('^\w*\#') #Allows full line comments with #

    #Need a temp file because some text must be inserted at the top
    #of the final cube file that we won't know until the end. 
    header=open('cubHeader.tmp','w')

    #initializing variables to count the number of points and other flags
    vec1points = vec2points = vec3points = blockCount = valCount = 0
    first=True
    vec2incremented = cycleFin = False

    for line in xyz: 
        if not comment.search(line) :
            if ( not line.strip() == '' ) :
                entry = [ float(x) for x in line.split() ]
                loc = entry[0:3]
                vstep = entry[3]

                #Counts the number of times 2 increments itself using the 
                #fact there is always a blank line in the file before it happens 
                if not ( 'x3inc' in locals() or 'x3inc' in globals() ) and \
                        cycleFin :
                    vec2points+=1


                #Responsible for figuring out the variables associated with the
                #grid. 
                if first :
                    origin=loc
                    lastLoc=loc
                    first=False
                
                #First increment will be in vector 1
                elif not ( 'x1inc' in locals() or 'x1inc' in globals() ):
                    x1inc=loc[0] - origin[0]
                    y1inc=loc[1] - origin[1]
                    z1inc=loc[2] - origin[2]

                #The increment after first block break will be in vector 2
                elif not ( 'x2inc' in locals() or 'x2inc' in globals() ) and \
                        blockCount == 1 :
                    x2inc=loc[0] - origin[0]
                    y2inc=loc[1] - origin[1]
                    z2inc=loc[2] - origin[2]

                #After a wraparound in vector 2, the increment will be in vector 3.
                #A block end marks an increment in vector 2 OR 3. 
                #To tell if it was vector 2 or 3 that was iterated through we 
                #see whether the increment is comparable to what we've already
                #determined vector 2's increment looks like.
                elif not ( 'x3inc' in locals() or 'x3inc' in globals() ) and \
                     cycleFin and \
                     not ( isClose(loc[0]-lastLoc[0],x2inc,1e-5,0) and \
                           isClose(loc[1]-lastLoc[1],y2inc,1e-5,0) and \
                           isClose(loc[2]-lastLoc[2],z2inc,1e-5,0) ):

                    x3inc=loc[0] - origin[0]
                    y3inc=loc[1] - origin[1]
                    z3inc=loc[2] - origin[2]


                #Counts number of points before 2 is incremented, so the number
                #of points in the first block of numbers. 
                if not ( 'x2inc' in locals() or 'x2inc' in globals() ):
                    vec1points+=1


                if cycleFin:
                    lastLoc=loc

                #After a full cycle of vector1 or after printing 
                #six values we need to go to a new line. 
                if cycleFin or valCount % 6 == 0: 
                    cub.write('\n')
                    cycleFin=False


                cub.write('%14.6e' % (vstep) )
                valCount += 1


            else:
                cycleFin=True
                valCount=0

                #This count will be off if there is no final empty line in the xyz
                blockCount+=1

    
    vec3points = blockCount/vec2points

    header.write('\n')
    header.write('INCREMENT FAST,MED,SLOW: X,Y,Z\n')
    header.write('%5d %12.8f %12.8f %12.8f \n' % ( atNum, origin[0], origin[1], origin[2] ) )
    header.write('%5d %12.8f %12.8f %12.8f \n' % ( vec3points, x3inc, y3inc, z3inc ) )
    header.write('%5d %12.8f %12.8f %12.8f \n' % ( vec2points, x2inc, y2inc, z2inc ) )
    header.write('%5d %12.8f %12.8f %12.8f \n' % ( vec1points, x1inc, y1inc, z1inc ) )
    header.close()
    cubname=cub.name
    cub.close()


        
    p=Popen('cat cubHeader.tmp '+cubname+' > finalCub.tmp',shell=True)
    p.wait()
    p=Popen('mv finalCub.tmp '+cubname,shell=True)
    p.wait()
    p=Popen('rm cubHeader.tmp',shell=True)
    p.wait()


#            In testing this code .xyz files have been found that do not increment xyz in the 
#            correct order. If the dimensions of the .cub file are incorrect check if this is 
#            the case. If not, the next most likely points of failure are lack of a new line
#            ending the file and the rel_tol variable for the isClose calls being a bad size. 

if __name__ == '__main__' :
    parser = argparse.ArgumentParser(description='This program takes as input a Turbomole xyz '+ 
            'file (with grid point coordinates followed by property values) and converts it '  +
            'to a cub format file. This xyz file is NOT the same as a coordinate xyz file.\n'  +

            'It does this by looking at the number of entries and the distance apart they '    +
            'are. It relies on the fact values start at the gridpoint with the smallest '      +
            'possible coordinates and the file iterates along x, y, z, in that order. Each '   +
            'y/z increment must be marked by a newline. The output file will be the input xyz '+
            "file's basename followed by .cub\n")

    parser.add_argument('xyzName',nargs=1, help='The name of the .xyz file that needs to be '+
                        'converted.')

    parser.add_argument('-c','--coord',default='./coord',help='Allows the user to specify the '+
                        'path to a Turbomole coordinate file. (Default: ./coord)') 

    args = parser.parse_args()


    xyz2cub(args.xyzName,args.coord)

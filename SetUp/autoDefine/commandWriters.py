import sys
import re



# This script complements autoDefine.py. 
# 
# Its member functions determine what input to print and
# then print that input to the define.input file
#
# FOR THOSE WHO MODIFY/ADD FUNCTIONALITY:
#
# All functions in commandWriters.py must have at least the following arguments:
# botSpecs, entries, defInp
# 
# Please only include extra variables if absolutely necessary, we want to avoid 
# variable creep in the main function. 
#
# When adding a new writing function, define both the defaults and key formats 
# the same way previous functions defined them. If this is done a template 
# input file may be generated by running the commands:
# 
# sed -n -e '/^#KeyFormat:/,/^[^#]/ { /^[^#]/b; p }' commandWriters.py > templateOptions
# sed -i -e 's|^#KeyFormat: ||' -e 's|^#||' templateOptions
# echo '$end' >> templateOptions
#
# And similarly a list of the Default values may be generated by:
#
# sed -n -e '/^#Default/,/^#KeyFormat:/ { /^#KeyFormat:/b; p }' commandWriters.py > defaultList


#getLine() finds the line with key on it and returns everything 
#          until the next '$', replacing newlines with spaces. 
#
#   Input:
#       opts - The list of specified options
#       key  - The marker for the section we want,
#              should have format '$sectionName'
#
#   Output:
#       target - The section specified by key with newlines
#                replaced by spaces
def getLine(opts, key):
   thisIter = iter(opts)
   target=None

    
   if opts:
       line=thisIter.next()
       noEnd=True
   else:
       noEnd=False
 
   while noEnd:
       #Find the section specified by key
       if line.split()[0] == key:
           target=line.split(key)[1]
           line=thisIter.next()

           #Get continued lines
           while '$' not in line:
               target=target.strip('\n')
               target=target+' '+line
               line=thisIter.next()
 
       #Push through until file end or key reached. 
       if line.split()[0] == '$end':
           noEnd = False
       elif line.split()[0] != key :
           line=thisIter.next()

   return target



#getEscapeChars() replaces the phrases used to represent escaped characters
#                 with those characters. Currently only needed/used for 
#                 basis set definitions. 
#
#   Input:
#       string - The string to be modified
#
#   Output:
#       string - The string with new substitutions
def getEscapeChars( string ):
    #For now replace SAVEUNDERSCORE and SAVEPERCENT with what they
    #would be in the options file. Need some way to specify 
    #_/% from directory though
    escapeDict={ 'SAVESPACE'      : ' ', 'SAVEEQUALS'  : '=', \
                 'SAVEUNDERSCORE' : ' ', 'SAVEPERCENT' : '='  }

    for phrase in escapeDict:
        string = string.replace( phrase, escapeDict[phrase] )

    return string
                 


#title() writes title
#
#Default: ''
#
#KeyFormat: $title [title]
def title(botSpecs, entries, defInp):
   key='$title'
    
   # The three sources for options are
   # in order: 
   # - The user specified/default options file
   # - The options file in directory we are currently
   #   iterating over.
   # - The options specified by the names of the supplied
   #   directories. 
   title=getLine(entries, key)
   titleBot=getLine(botSpecs, key)

   # precedence: title < titleBot < titleDir
   if titleBot:
      title = titleBot

   if not title:
     defInp.write('\n')
   else:
     title=title.strip()
     defInp.write(title+'\n')



#readCoord() tells the name of the coordinate file
#
#Default: only $coord present/$coord omitted - file=coord
#
#KeyFormat: $coord file=[fileName]
def readCoord(botSpecs, entries, defInp):
   key='$coord'

   coord=getLine(entries, key)
   coordBot=getLine(botSpecs, key)

   if coordBot:
      coord = coordBot
   
   if coord:
      if 'file=' in coord:
         coord=coord.split('file=')[1]
         coord = getEscapeChars( coord )
         defInp.write('a '+coord)
      else:
         defInp.write('a coord\n')
   else:
      defInp.write('a coord\n')



################################################################################
#                          Begin Coordinate Menu Writers                       #
################################################################################



#assignSym() tells symmetry to use, groups must have d6h
#symmetry or a group deriving from some subset of D6h 
#operations. (Is this correct?)
#'sym=auto' will attempt to automatically determine symmetry. 
#
#Default: c1
#
#KeyFormat: $sym sym=[auto/character] eps=[float]
def assignSym(botSpecs, entries, defInp):
   key='$sym'

   sym=getLine(entries, key)
   symBot=getLine(botSpecs, key)

   if symBot:
      sym = symBot

   if sym:
      if 'sym=' in sym:
         group=sym.split('sym=')[1]
         group=group.split()[0]
         if 'eps=' in sym:
            eps=sym.split('eps=')[1]
            eps=eps.split()[0]
            if group=='auto':
               defInp.write('desy '+eps+'\n')
            else:
               defInp.write('sy '+group+' '+eps+'\n')
         else:
            if group=='auto':
               defInp.write('desy\n')
            else:
               defInp.write('sy '+group+'\n')



#detInternals() tells whether internal redundant coordinates
#should be used. It returns True if they are. 
#
#Default: 'on'
#
#KeyFormat: $internal [on/off]
def detInternals(botSpecs, entries, defInp):
   key='$internal'

   internal=getLine(entries, key)
   internalBot=getLine(botSpecs, key)

   if internalBot:
      internal = internalBot

   if internal:
      internal=internal.strip()
      if internal == 'on':
         defInp.write('ired\n')
         return True
   else:
      defInp.write('ired\n')
      return True

   return False

#assignFrags() used to define fragments. It looks like the define section 
#for this has some issues so if $frag is used, you must define the 
#the fragment number for all atoms. Also, the symmetry assignments for
#individual fragments don't seem to work in define, so auto will always 
#be used. 
#
#Defaults: '$frag' omitted - No fragments
#          '$frag' used    - Must define atoms in fragments, no default
#                          - charge=0 for all fragments
#                          - auto symmetry determination for each fragment
#
#KeyFormat: $frag frag1=[1-4,6] [ frag2=[1-4,6] frag3=[1-4,6] ] 
#            chrg1=[int] chrg2=[int] chrg3=[int] 
#            sym=[auto/none/character] *currently disabled*
def assignFrags(botSpecs, entries, defInp):
   key='$frag'

   frag=getLine(entries, key)
   fragBot=getLine(botSpecs, key)

   if fragBot:
      frag = fragBot

   if frag:
      defInp.write('frag\non\nq\nq\n')
      if 'frag1=' in frag:
         frag1=frag.split('frag1=')[1]
         frag1=frag1.split()[0]
         ats=frag1.split(',')
         for at in ats:
            defInp.write('x\n')
            at.strip()
            rnge=at.split('-')
            if len(rnge) == 2:
               defInp.write(rnge[0]+'\n'+rnge[1]+'\n1\n')
            elif len(rnge) == 1:
               defInp.write(rnge[0]+'\n'+rnge[0]+'\n1\n')
            else:
               print 'Error: invalid entry in $frag'
               sys.exit()
      if 'frag2=' in frag:
         frag2=frag.split('frag2=')[1]
         frag2=frag2.split()[0]
         ats=frag2.split(',')
         for at in ats:
            defInp.write('x\n')
            at.strip()
            rnge=at.split('-')
            if len(rnge) == 2:
               defInp.write(rnge[0]+'\n'+rnge[1]+'\n2\n')
            elif len(rnge) == 1:
               defInp.write(rnge[0]+'\n'+rnge[0]+'\n2\n')
            else:
               print 'Error: invalid entry in $frag'
               sys.exit()
      if 'frag3=' in frag:
         frag3=frag.split('frag3=')[1]
         frag3=frag3.split()[0]
         ats=frag3.split(',')
         for at in ats:
            defInp.write('x\n')
            at.strip()
            rnge=at.split('-')
            if len(rnge) == 2:
               defInp.write(rnge[0]+'\n'+rnge[1]+'\n3\n')
            elif len(rnge) == 1:
               defInp.write(rnge[0]+'\n'+rnge[0]+'\n3\n')
            else:
               print 'Error: invalid entry in $frag'
               sys.exit()

#      if 'sym=' in frag:
#         sym=frag.split('sym=')[1]
#         sym=sym.split()[0]
#         if sym == '':

      if 'chrg' in frag:
         defInp.write('cha\n')
      if 'chrg1=' in frag:
         chrg1=frag.split('chrg1=')[1]
         chrg1=chrg1.split()[0]
         defInp.write(chrg1+'\n')
      if 'chrg2=' in frag:
         chrg2=frag.split('chrg2=')[1]
         chrg2=chrg2.split()[0]
         defInp.write(chrg2+'\n')
      if 'chrg3=' in frag:
         chrg3=frag.split('chrg3=')[1]
         chrg3=chrg2.split()[0]
         defInp.write(chrg3+'\n')

      defInp.write('\n\n\n')
      



################################################################################
#                          Begin Basis Set Menu Writers                        #
################################################################################


#defBasis() tells what basis set to use for each atom type. Note that
#the basis set names are case-sensitive. More options may be implemented
#in the future. 
#
#Default: '$basis' omitted - all def2-SV(P)
#         No basis specified for an atom - def2-SV(P)
#         No atoms specified for a basis - all atoms
#         atoms specified with no basis  - don't do this
#
#KeyFormat: $basis  [def2-QZVP]=[all/1-4,6/"c"] ...
def defBasis(botSpecs, entries, defInp):
   key='$basis'

   basis=getLine(entries, key)
   basisBot=getLine(botSpecs, key)

   if basisBot:
      basis = basisBot

   #Assigning the basis set information
   if basis:
      basis = basis.split()
      for entry in basis:
         defInp.write('b\n')
         if '=' in entry:
            entry = entry.split('=')
            if entry[1] != '':
                entry[0] = getEscapeChars(entry[0])
                entry[1] = getEscapeChars(entry[1])
                defInp.write(entry[1]+' '+entry[0]+'\n')
            else:
                print 'Error: empty basis set assignment for '+entry[0]
                sys.exit()
         else:
            entry = getEscapeChars(entry)
            defInp.write('all '+entry+'\n')


################################################################################
#                  Begin Molecular Orbital Calculation Writers                 #
################################################################################

#useHcore() turns hcore guess on or off.
#
#Default: off
#
#KeyFormat: $hcore [on/off]
def useHcore(botSpecs, entries, defInp):
   key='$hcore'

   hcore=getLine(entries, key)
   hcoreBot=getLine(botSpecs, key)

   if hcoreBot:
      hcore = hcoreBot

   # future input depends on hcore 
   if hcore:
      hcore = hcore.strip()
      if hcore == 'on': 
         defInp.write('hcore\n')
         return True

   return False



#eht() manage Hueckel guess options. EHT is always used
#unless hcore is used instead, it can't be turned off. 
#
#Default: global_constant=1.70 mod_Wolfsberg-Helmholz=off
#
#KeyFormat: $eht global=[float] modWH=[on/off] 
#               (atom index pair format -->)[int,int]=[float]
def eht(botSpecs, entries, defInp):
   key='$eht'

   eht=getLine(entries, key)
   ehtBot=getLine(botSpecs, key)

   if ehtBot:
      eht = ehtBot

   defInp.write('eht\n')

   modAtomEht=False
   #If key present
   if eht:
      eht = eht.split()
      #Should be 2-D
      eht = [ entry.split('=') for entry in eht ]
      #Just in case only $eht specified
      if len(eht) == 0:
         defInp.write('y\n')
         return
      else:
         defInp.write('n\n')
         #Set global Hueckel constant first
         if 'global' in [item for sub in eht for item in sub]:
            for entry in eht:
               if len(entry) != 2:
                  print 'Error: improper eht option supplied'
                  sys.exit()
               if entry[0] == 'global':
                  defInp.write('y\n'+entry[1]+'\n')
         else:
            defInp.write('n\n')

         #followed by deciding whether to use the 
         #modified Helmholz-Wolfsberg formula. 
         if 'modWH' in [item for item in sub for sub in eht]:
            for entry in eht:
               if len(entry) != 2:
                  print 'Error: improper eht option supplied'
                  sys.exit()
               if entry[0] == 'modWH' and entry[1] == 'on':
                  defInp.write('y\n')
               else:
                  defInp.write('n\n')
         else:
            defInp.write('n\n')

         #Then check for any entry with a comma in the first
         #range. It's assumed these are atomic indices. 
         for entry in eht:
            if len(entry) != 2:
               print 'Error: improper eht option supplied'
               sys.exit()
            if ',' in entry[0]:
               modAtomEht=True
               defInp.write('y\n')
               indices=entry[0].split(',')
               defInp.write(indices[0]+' '+indices[1]+'\n'+entry[1]+'\n')

         #Needed for when it asks about printing out
         #Hueckel coeffs. Default is no. 
         if modAtomEht:
            defInp.write('\n')

   #Will either skip the Hueckel modifying section or exit
   #the atom pair section depending on choices. 
   defInp.write('\n')




#molCharge() define the molecular charge
#
#Default: 0
#
#KeyFormat: $charge [int]
def molCharge(botSpecs, entries, defInp):
   key='$charge'

   charge=getLine(entries, key)
   chargeBot=getLine(botSpecs, key)

   if chargeBot:
      charge = chargeBot

   if charge:
      charge = charge.strip()
      defInp.write(charge)

   defInp.write('\n')



#setOcc() handles options involving electron occupation. 
#
#Default: lowest energy occupation
#
#KeyFormat: $occ [int]
def setOcc(botSpecs, entries, defInp):
   key='$occ'

   occ=getLine(entries, key)
   occBot=getLine(botSpecs, key)

   if occBot:
      occ = occBot

   defInp.write('\n')

################################################################################
#                          Begin Method Definition Writers                     #
################################################################################



#dft() handles specifications within the density 
#functional theory menu.
#
#Default: $dft omitted - DFT not used
#         $dft used    - functional b-p
#                        gridsize   m3
#
#KeyFormat: $dft func=[label] grid=[label]
def dft(botSpecs, entries, defInp):
   key='$dft'

   dft=getLine(entries, key)
   dftBot=getLine(botSpecs, key)

   if dftBot:
      dft = dftBot

   if dft:
      defInp.write('dft\non\n')
      if 'func=' in dft:
         func=dft.split('func=')[1]
         func=func.split()[0]
         defInp.write('func '+func+'\n')
      if 'grid=' in dft:
         grid=dft.split('grid=')[1]
         grid=grid.split()[0]
         defInp.write('grid '+grid+'\n')
      
      defInp.write('\n')



#ri() handles specifications within the resolution 
#of identity menu.
#
#Default: $ri omitted - ri not used
#         $ri used    - ricore 500 Mb
#                     - jbas matches basis type
#                     - file named auxbasis holds jbas
#                     - when specifying jbas, same defaults 
#                       as specifying basis used
#
#KeyFormat: $ri mem=[int in Mb] file=[name] 
#           jbas=[def2-QZVP]=[all/1-4,6/"c"] ...
def ri(botSpecs, entries, defInp):
   key='$ri'

   ri=getLine(entries, key)
   riBot=getLine(botSpecs, key)

   if riBot:
      ri = riBot

   if ri:
      defInp.write('ri\non\n')

      #specify $ricore
      if 'mem=' in ri:
         mem=ri.split('mem=')[1]
         mem=mem.split()[0]
         defInp.write('m '+mem+'\n')

      # modify file name for $jbas
      if 'file=' in ri:
         fil=ri.split('file=')[1]
         fil=fil.split()[0]
         defInp.write('f '+fil+'\n')

      # Set the $jbas basis type
      if 'jbas=' in ri:
         defInp.write('jbas\n')
         jbas = ri.split('jbas=')[1:]
         for entry in jbas:
            entry = entry.split()[0]
            if '=' in entry:
               ats  = entry.split('=')[1]
               entry = entry.split('=')[0]
            else:
               ats = 'all'
            entry = getEscapeChars( entry )
            defInp.write('b\n'+ats+' '+entry+'\n')
         defInp.write('*\n')

      defInp.write('\n')


   
#cc() handles specifications within the cc menu. 
#There are quite a few options, may want to break
#it up more. For now only those options relevant to
#rirpa will be implemented. I kind of think cc should
#be more than one menu in define itself. 
#
#A note on freeze=: pick one from num/energy/list exactly 
#as written. Options can't be combined. 
#input for num    is number of frozen core orbitals
#          energy is the orbital energy cutoff value in Hartree
#          list   is orbital specification via [1-4,6] style format
#
#Default: $cc used/omitted - maxcor=500mb
#                            denconv=1d-7
#   freeze with no options - energy=-3H
#     cbas with no options - cbas matches basis type
#
#KeyFormat: $cc freeze=[num/energy/list]=[input] mem=[int in Mb]
#               cbas=[def2-QZVP]=[all/1-4,6/"c"] ...
#               denconv=[float]
def cc(botSpecs, entries, defInp):
   key='$cc'

   cc=getLine(entries, key)
   ccBot=getLine(botSpecs, key)

   if ccBot:
      cc = ccBot

   if cc:
      defInp.write('cc\n')

      # Handles the freeze sub-menu
      if 'freeze' in cc:
         defInp.write('freeze\n')
         if 'freeze=' in cc:
            freeze = cc.split('freeze=')[1]
            freeze = freeze.split()[0]
            if 'num=' in freeze:
               num = freeze.split('num=')[1]
               defInp.write('core '+num+'\n')
            elif 'energy=' in freeze:
               energy = freeze.split('energy=')[1]
               defInp.write('fp '+energy+'\n')
            elif 'list=' in freeze:
               lis = freeze.split('list=')[1]
               defInp.write('core 0\nf '+lis+'\n')
         defInp.write('*\n')

      # Handles the cbas sub-menu
      if 'cbas' in cc:
         defInp.write('cbas\n')
         if 'cfail' in cc:
            defInp.write('\n')
         if 'cbas=' in cc:
            cbas = cc.split('cbas=')[1:]
            if not cbas[0].split()[0] == 'default':
                for entry in cbas: 
                   entry = entry.split()[0]
                   if '=' in entry:
                      ats = entry.split('=')[1]
                      bas = entry.split('=')[0]
                   else:
                      ats = 'all'
                      bas = entry
                   bas = getEscapeChars( bas )
                   defInp.write('b\n'+ats+' '+bas+'\n')
         defInp.write('*\n')

      # Sets $maxcor
      if 'mem=' in cc:
         mem = cc.split('mem=')[1]
         mem = mem.split()[0]
         defInp.write('memory '+mem+'\n')

      # Sets $denconv
      if 'denconv=' in cc:
         dens = cc.split('denconv=')[1]
         dens = dens.split()[0]
         defInp.write('denconv '+dens+'\n')
      
             
      defInp.write('*\n')




#rirpa() handles rirpa related specifications.
#
#Default: $rirpa omitted - rirpa not used
#         $rirpa used    - npoints 60, all other settings off
#
#KeyFormat: $rirpa npoints=[int] rpagrad nohxx rpaprof 
def rirpa(botSpecs, entries, defInp):
   key='$rirpa'

   rirpa=getLine(entries, key)
   rirpaBot=getLine(botSpecs, key)

   if rirpaBot:
      rirpa = rirpaBot

   if rirpa:
      defInp.write('rirpa\n')

      if 'npoints=' in rirpa:
         npoint = rirpa.split('npoints=')[1]
         npoint = npoint.split()[0]
         defInp.write('npoints '+npoint+'\n')
         
      if 'rpagrad' in rirpa:
         defInp.write('rpagrad\n')

      if 'nohxx' in rirpa:
         defInp.write('nohxx\n')

      if 'rpaprof' in rirpa:
         defInp.write('rpaprof\n')

      defInp.write('\n')
      



#scf() handles scf related specifications.
#
#Default: convergence criteria     : conv=6 --> 1d-6
#         integral threshold values: thize=1d-5, thime=5
#         2e- integral storage     : ints=None
#         scf iteration limit      : iter=60
#         
#
#KeyFormat: $scf  conv=[int] iter=[int] 
def scf(botSpecs, entries, defInp):
   key='$scf'

   scf=getLine(entries, key)
   scfBot=getLine(botSpecs, key)

   if scfBot:
      scf = scfBot

   if scf:
      defInp.write('scf\n')

      if 'conv=' in scf:
         conv=scf.split('conv=')[1]
         conv=conv.split()[0]
         defInp.write('conv\n'+conv+'\n')

      if 'iter=' in scf:
         iters=scf.split('iter=')[1]
         iters=iters.split()[0]
         defInp.write('iter\n'+iters+'\n')


      defInp.write('\n')

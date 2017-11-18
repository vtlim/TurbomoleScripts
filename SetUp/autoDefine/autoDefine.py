#! /usr/bin/python

import sys
import os
import re
import argparse
import subprocess as sp
from shutil import copyfile

import commandWriters as cw


#This script is meant to automate the set up of Turbomole 
#calculations using a short list of input options. For 
#information on specifying these options see the accompanying
#commandWriters.py file and/or run 'autoDefine.py -h'
#
#This script has not been too thoroughly tested and may fail
#with some versions of Turbomole and for certain edge cases.
#All output is saved in a def.out file which may be 
#examined for help in debugging. 
#
#primary author: Matt Agee (magee@uci.edu)
#
#The structure of this program is as follows:
#
#   - inputBuilder(): Acts as a controller. Its function
#     is to call other functions in the correct order. 
#     These other functions will either parse the input 
#     file, directory structure or write the necessary 
#     commands to the define input file. 
#
#   - optsParser(): Responsible for gathering the 
#     run specifications. It does this by reading the main
#     'options' file, by default located in the 
#     working directory and/or the directory specific 
#     'options' located in the directories supplied as 
#     argument
# 
#   - Command Printers: The options from each 'options' file
#     will be stored in lists.
#     These lists will be passed to command printers.
#     These command printers will find the section relevant to
#     them, parse the options present in their 
#     relevant section, and print the commands define needs to 
#     generate the control file. 
#     It is recommended there be one command printer function 
#     for each menu option. Consider making further functions 
#     for each menu suboption. 
#
#It is expected this will be a collaborative project. So the 
#protocol for adding a new command printer is as follows:
#
#   - Edit the optsParser() function so it looks for the new
#     section you define. 
#
#   - The section key should match as closely as possible the one
#     in the control file it concerns. The idea is to use 
#     a control file as a template for the options file. But this does not 
#     always make sense so deviate as you see fit. 
#
#   - Your function should be defined in commandWriters.py
#
#   - Find a sensible place to put your function call in inputBuilder().
#     Don't do something like try to define coord menu options in the 
#     middle of defining basis set options. 
#     An '&' should not appear in the define input file without
#     very good reason. 
#
#   - Your function should take at least three arguments:
#     - A list of strings containing options from the user-specified
#       options file.
#     - A list of strings containing options from the directory specific
#       options file.
#     - an opened, writable, input file. 
#
#   - When printing your commands, make sure to finish with a new line
#     character. You may also assume you will be starting on a new 
#     line when first writing your commands. 
#
#   - If no options match those you are looking for your function should
#     do nothing unless there is good reason to have some sort of default. 
#
#That's pretty much it. Hopefully following these rules will keep the code
#readable. 



#getEscapeChars() replaces escaped characters with phrases used to represent
#                 them. Used just before adding an entry to the list of 
#                 options. 
#
#   Input:
#       string - The string to be modified
#
#   Output:
#       string - The string with new substitutions
def escapeChars( string ):
    escapeDict={ 'SAVESPACE'      : '^ ', 'SAVEEQUALS'  : '^=', \
                 'SAVEUNDERSCORE' : '^_', 'SAVEPERCENT' : '^%'  }

    for phrase in escapeDict:
        string = string.replace( escapeDict[phrase], phrase )

    return string



#optsParser() Looks through a supplied options file, storing the valid section 
#             keywords and their corresponding subsections in a list. 
#
#   Input:
#       opts - The opened options file
#
#   Output:
#       entries - The list with sections and their subsections. 
def optsParser( opts ):
    #This is a list of every keyword this program is able to 
    #recognize. Put your keyword here if you wish to expand the
    #program. Key words must start with '$'. 
    #
    #Subsections for a recognized keyword are automatically stored. 
    #To add a new subsection, find that subsection's corresponding
    #function in commandWriters.py
    allowedKeys=[ '$title', '$coord', '$sym', '$internal', '$frag', '$basis', 
            '$hcore', '$eht', '$charge', '$occ', '$dft', '$ri', '$cc', 
            '$rirpa', '$scf', '$dsp', '$fix', '$cosmo' ] 
    entries=[]

    line=opts.next()
    noEnd=True
    while noEnd:
        #Add lines to array until new section is reached. 
        if line.split()[0] in allowedKeys:
            line = escapeChars( line )
            entries.append(line)
            line=opts.next()
            while '$' not in line:
                entries.append(line)
                line = opts.next()

        #Push through until file end
        if line.split()[0] == '$end':
            entries.append('$end')
            noEnd = False
        elif line.split()[0] not in allowedKeys:
            print('Warning: ',line.split()[0],' is not an allowed key.')
            line=opts.next()


    return entries



#inputBuilder() Manages creation of the define input file and runs define. 
#               It calls all the functions for parsing input and writing 
#               define commands. 
#
#   Input:
#       directories   - A list of strings containing the names of target 
#                       directories. 
#       options       - A string with the name of the user-defined options file. 
#       keep_going    - If true run define for all directories even if define fails 
#                       for one. 
#       save_intermed - If true do not remove files generated in setting up control
def inputBuilder( directories, options, keep_going, save_intermed ): 
    p=sp.Popen('pwd',stdout=sp.PIPE,shell=True)
    workDir=p.communicate()[0].strip()

    entries=''
    if os.path.exists(options):
        opts=open(options,'r')
        entries=optsParser(opts)
        opts.close()


    if entries == '':
        response = raw_input('Warning: no valid options in '+options+'\n'+
                             'Continue anyways? [y/n]\n')
        while response.strip() != 'y':
            if response.strip() == 'n':
                sys.exit()
            else:
                response = raw_input("Sorry I didn't understand that.\n"+
                                     'Continue? [y/n]\n')




    #Iterate through all argument directory paths. If 
    #no directories were supplied, do this for the working
    #directory. 
    for dirs in directories:
        #getting options specified in options file located at 
        #bottom of dirs
        botSpecs=''
        if os.path.exists(dirs+'/options'):
            botOpts=open(dirs+'/options','r')
            botSpecs=optsParser(botOpts)
            botOpts.close()

        #One of these will be made for each directory supplied as
        #argument, then moved to that directory. 
        defInp=open('def.input','w')


        #No pre-existing control files allowed, changes command
        #pattern too much. 
        defInp.write('\n')


        #Begin the command printers. These should be defined within
        #their own scripts to reduce clutter here. Be sure to mention
        #the proper syntax for providing them keywords there.
        cw.title(botSpecs, entries, defInp)


        #Coordinate definition menu
        cw.readCoord(botSpecs, entries, defInp)

        cw.assignSym(botSpecs, entries, defInp)

        cw.fix(botSpecs, entries, defInp)

        cw.detInternals(botSpecs, entries, defInp)

        cw.assignFrags(botSpecs, entries, defInp)

        defInp.write('*\n\n')


        #Basis set definition menu
        cw.defBasis(botSpecs, entries, defInp)

        defInp.write('*\n')
 

        #Molecular orbital calculation menu
        hcore=cw.useHcore(botSpecs, entries, defInp)

        #Only run through eht if hcore not used. 
        if not hcore:
            cw.eht(botSpecs, entries, defInp)

        cw.molCharge(botSpecs, entries, defInp)

        cw.setOcc(botSpecs, entries, defInp)

        #Just to make sure we're through the previous menu. 
        defInp.write('\n\n\n')
        

        #Method definition menu
        cw.dft(botSpecs, entries, defInp)

        cw.ri(botSpecs, entries, defInp)

        cw.cc(botSpecs, entries, defInp)

        cw.rirpa(botSpecs, entries, defInp)

        cw.scf(botSpecs, entries, defInp)

        defInp.write('*')
        defInp.close()


        #Now move the define input file over and run define with it. 
        os.rename('def.input',dirs+'/def.input')
        os.chdir(dirs)
        if os.path.exists('control'):
            os.remove('control')

        p=sp.Popen('define < def.input > def.out',shell=True)
        p.wait()

        if not save_intermed:
            os.remove('def.input')
            os.remove('def.out')
        if any("cosmo" in s for s in entries):
            cosInp=open('cosmoprep.input','w')
            cw.cosmo(botSpecs, entries, cosInp)
            cosInp.close()
            p=sp.Popen('cosmoprep < cosmoprep.input > cosmoprep.out',shell=True)
            p.wait()


        #Place setup files no longer needed in setup dir
        if not os.path.exists('setup'):
            os.makedirs('setup')
        for i in ['input.xyz','options','cosmoprep.input','cosmoprep.out','def.input','def.out']:
            try: os.rename(i,'setup/'+i)
            except OSError: pass

        #Check if define finished, if the user did not choose
        #to ignore failed set ups, exit program 
        if os.path.exists('tmp.input'):
            print('Error: Define failed in '+dirs+'\n')
            if not keep_going:
                sys.exit()
        else:
            #Post processing of control done here
            cw.dsp(botSpecs, entries, keep_going)

        os.chdir(workDir)




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Tool for setting up Turbomole'
           + ' calculations. Takes as argument a list of directories containing'
           + ' coord files. Will delete any control file already present, and '
           + ' run define based on options specified in the global options'
           + ' file (called ./options or specified with option -o), or through'
           + ' the local options file (called options and located in'
           + ' the argument directory).')
    parser.add_argument('dirs',nargs='*',default='.', help='The directories'
           + ' containing the coord files to be used with define.'
           + ' (ex: dir1 dir2 ...) (Default: current directory) ')
    parser.add_argument('-o','--options', default='options',help='Allows the'
           + ' user to specify the path to the global options file.'
           + ' (ex: -o global_options) (Default: ./options)')
    parser.add_argument('-k','--keep_going', action='store_true',
            help='Continue iterating through directories even if define fails.')
    parser.add_argument('-i','--save_intermed', action='store_true',
            help='Do not remove files generated while setting up control')
    args = parser.parse_args()

    
    inputBuilder(args.dirs, args.options, args.keep_going, args.save_intermed)

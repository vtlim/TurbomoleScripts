#! /bin/bash


if [ "${BASH_SOURCE[0]}" == "${0}" ]; then 

if [ "$1" = "-h" ]; then
  cat << HERE
USAGE: $0 [-h] [-r] [-d] [-n] [-x] [-c]  DIRECTORIES

This script takes as input a list of directories. 
It then sets up the slurm batch script for each directory,
copies that script to the directory and submits it. 
It names the job based off the supplied directory paths. 
You must have a slurm_template file set up and linked appropriately
within the bash script for this to work. 

The command to issue and most other settings must be specified
in the script. This script was made for running Turbomole
calculations. It should work for other calculation types
but may require adjustment. 

Memory is calculated by adding the control file's
\$ricore and \$maxcor specifications with user-specified buffer added
on top. It takes into account multiple processors being used. 

If TASK_PER_NODE is greater than 1 PARALLEL_TEMPLATE is used. 

Options:
-r              restart, only argument directories that have a file contained 
                in RESTART_FILES will be submitted. 

-d              debug, the batch script won't be submitted. The script just 
                creates the submit script and moves it to the target directory. 

-n              NumForce, the amount of memory requested won't be divided 
                by the number of cores you request. Intended for NumForce 
                calculations. 

-x              exclusive, --exclusive will be included in the list of SBATCH 
                options. In theory this will make the job run as the only job 
                on its node, but it may not work on green planet. 

-c              no control, the script will not require a control file to be in
                the argument directories. To use this option you must specify 
                MEM_PER_CPU at the top of the script. 
HERE

  exit 0 
fi

# Path information needed to find Slurm templates
# SLURM_BATCH_DIR is the directory holding files
# SLURM_TEMPLATE and PARALLEL_TEMPLATE
SLURM_BATCH_DIR="{BASH_SOURCE[0]}"
SLURM_TEMPLATE='slurm_template'
PARALLEL_TEMPLATE='slurm_template_parallel'

#The memory beyond that specified in $ricore and $maxcor in MB
#set this to make sure there is memory left over for stuff like 
#BLAS and LAPACK. 
BUFFER=1000

#The options to go in OPTIONS, names should be self-explanantory
THE_NODES='1'
TASK_PER_NODE='8'
CPU_PER_TASK='1'
MEM_PER_CPU='4000' #In MB, only used if '-c' option is set
TIME_HR_MN_SC='200:00:00'
THE_PARTITION='mf_nes2.8'
CHOOSE_MPI_OR_SMP='SMP' #Make sure you set this if you're running parallel


#COMMAND=( 'jobex -ri' )
#COMMAND=( 'jobex -ri -energy 7 -gcart 4' )
#COMMAND=( 'jobex -ri -md > md.out' )
#COMMAND=( 'jobex -rijk -level cc2 -c 50 -energy 7 -gcart 4' )
#COMMAND=( 'jobex -ri -level rirpa -l /export/home/magee/myRPAbranch/bin/em64t-unknown-linux-gnu 
#                 -energy 7 -gcart 4 > jobex.out' )
COMMAND=( 'jobex -ri -dscf -level rirpa -l /export/home/magee/rirpa_binaries
                 -energy 7 -gcart 4 > jobex.out' )


#COMMAND=( 'TURBODIR=~/V7-0-branch; TURBOIMG=~/V7-0-branch \n
#           export TURBODIR TURBOIMG \n
#           srun hostname > machines.tm\n
#           ~/bin/NumGrad > numGrad.out' )
#COMMAND=( 'ridft >> ridft.out' )
#COMMAND=( 'jobbsse -ri -setup' )
#COMMAND=( 'touch new_file VPQ_ISQR ; sleep 360' )


#COMMAND=( 'rirpa >> rirpa.out' )
#COMMAND=( '~/rirpa_binaries/rirpa_mkl > rirpa.out' )


#COMMAND=( "srun hostname > machines.turbomole\n
#           NumForce -ri -c -central -level rirpa -mfile machines.turbomole -scrpath /work/magee >> NumForce.out" )
#COMMAND=( "srun hostname > machines.turbomole\n
#           NumForce -ri -cosmo -central -mfile machines.turbomole -scrpath /work/magee >> NumForce.out" )
#COMMAND=( "srun hostname > machines.turbomole\n
#           NumForce -ri -c -central -mfile machines.turbomole -scrpath /work/magee >> NumForce.out" )

echo "COMMAND to be executed: $COMMAND"


######################
# The possible job submission partitions are listed here.
# Name:       CPU type:      Cores(Threads)/Node:  Speed:
# mf_m-c1.9   1.9GHz AMD       48                  0.5
# mf_ilg2.3   2.3GHz AMD       32(64)              0.7
# mf_i-b2.8   2.8GHz Intel     40                  1.0
# mf_nes2.8   2.8GHz Intel     8(16)-12(24)        1.0
#
# Slurm schedules each processor core to a single job only. Multithreaded jobs
# can take advantage of multiple threads per core. 
# 
# The 2.3GHz AMD nodes have 32 modules, each with 2 integer cores and 1 floating-
# point core. Each module is scheduled as a single core with two threads.
#
###############################################################################


#identifying info and other options likely to be iterated over

#The base job name to which the path name will be appended if desired
THE_JOB_NAME=''

#The number of directories up the argument directories in $@ needed to 
#uniquely specify your job
P_SAVED=6

###############################################################################
#Flag variables for the various options

#restart, set '-r' and only argument directories that have a 
#file contained in RESTART_FILES will be submitted. 
rstart=false

#debug, set '-d' and the batch script won't be submitted. The program 
#should just create a submit script and move it to the target directory. 
debug=false

#NumForce, set '-n' and the amount of memory requested won't be divided 
#by the number of cores you're going to get. 
numfor=false

#exclusive, set '-x' and --exclusive will be included in the 
#list of SBATCH options. In theory this will make the job run as the only 
#job on its node, but it may not work on green planet. 
exclusive=false

#no control, set '-c' and the script will not require a control file to be in
#the argument directories. To use this option you must specify MEM_PER_CPU at
#the top of the script. 
nocontrol=false

#NO LONGER AN OPTION
#Check to see if job is trying to use more than one processor. If so, it will
#then load up a different slurm batch script template. I may want to modify
#this later to just use the one. 
parallel=false

if (( $TASK_PER_NODE > 1 )); then
   parallel=true
fi

while getopts "rdnxc" OPT; do
   case $OPT in
      r) 
         rstart=true
         ;;

      d)
         debug=true
         ;;

      n)
         numfor=true
         ;;

      x)
         exclusive=true
         ;;

      c)
         nocontrol=true
         ;;

      *)
         eval BADOPT=\$$((OPTIND - 1))
         echo "Option $BADOPT is not yet implemented"
         exit 1
         ;;

   esac
done

shift $(( OPTIND - 1 ))

fi

###############################################################################

#defining functions to be used later in the main body. 

#this function checks if any of the RESTART_FILES exist in the argument
#directory. If they do it echoes 'true'. 
check_restart(){
   RESTART_FILES=( 'GEO_OPT_RUNNING' 'restarthess' 'GEO_OPT_FAILED' )

   for i in ${RESTART_FILES[@]}; do
      if [[ -a $1/$i ]]; then
         echo true
         return 0
      fi
   done
   echo false
}

###############################################################################

if [ "${BASH_SOURCE[0]}" == "${0}" ]; then 
    
OPTIONS=( "$THE_PARTITION" "$THE_NODES" "$TASK_PER_NODE" "$CPU_PER_TASK" \
          "$TIME_HR_MN_SC" "$CHOOSE_MPI_OR_SMP" "$COMMAND" )

OPTIONSREF=( 'THE_PARTITION' 'THE_NODES' 'TASK_PER_NODE' 'CPU_PER_TASK' \
          'TIME_HR_MN_SC' 'CHOOSE_MPI_OR_SMP' 'COMMAND' )

NUMOPTIONS=${#OPTIONSREF[@]}

LAST_OPT=''

WORKDIR=`pwd`

#Must have a SLURM_TEMPLATE batch script with the OPTIONSREF entries written where 
#the OPTIONS are to be placed
if [[ $parallel = true ]]; then 
   SLURM_TEMPLATE="$PARALLEL_TEMPLATE"
fi

cp "$SLURM_BATCH_DIR/$SLURM_TEMPLATE" temp_slurm_


#Preparing the SLURM submit script with all shared options
#
#We use printf instead of the commented out seds  
#to avoid interpretation of special characters in 
#the OPTIONS array
for i in `seq 0 $((NUMOPTIONS - 1))`; do
   while read -r line; do
      printf "${line/${OPTIONSREF[$i]}/${OPTIONS[$i]}}\n"
   done < "temp_slurm_${LAST_OPT}" > temp_slurm_${OPTIONSREF[$i]}
   rm temp_slurm_${LAST_OPT}
   LAST_OPT=${OPTIONSREF[$i]}
#   OPTIONS[$i]=`echo ${OPTIONS[$i]} | sed "s/|/SOME_REALLY_LONG_STRING_NOT_IN_COMMAND/g"`
#   sed "s|${OPTIONSREF[$i]}|${OPTIONS[$i]}|" temp_slurm_${LAST_OPT} > temp_slurm_swp
#   sed "s/SOME_REALLY_LONG_STRING_NOT_IN_COMMAND/|/g" temp_slurm_swp > temp_slurm_${OPTIONSREF[$i]}
#   rm temp_slurm_${LAST_OPT} temp_slurm_swp
#   LAST_OPT=${OPTIONSREF[$i]}
done

mv temp_slurm_$LAST_OPT temp_slurm_hold

if [[ $exclusive = true ]]; then
   sed -i 's|# SBATCH --exclusive|#SBATCH --exclusive|' temp_slurm_hold
fi

#The default values for the two forms of memory.
#ricore might not be used leading to excessive amounts, maybe examine
#COMMAND later to decide it. 
maxcor=500
ricore=500

for i in $@; do
   echo '////////////////////////////////////////////////////////////'

   if [[ ! -f $i/control && $nocontrol = false ]]; then
      echo "no control in $i"

   # If GEO_OPT_CONVERGED found assume calculation done even if some 
   # RESTART_FILES are present. 
   elif [[ ( $rstart = false ) || ( ($rstart = true) && ( `check_restart "$i"` = true ) && ! -f "$i"/GEO_OPT_CONVERGED ) ]]; then

      if [[ $nocontrol = false ]]; then
         #Setting the amount of memory per cpu so it matches the amount
         #of memory TURBOMOLE asks for. The 500 are default values used. 
         maxcor=`grep '$maxcor' "$i"/control | sed "s|\S*\s*\([0-9]*\)|\1|"`
         if [[ -z $maxcor ]]; then
            maxcor=500
         fi
         ricore=`grep '$ricore' "$i"/control | sed "s|\S*\s*\([0-9]*\)|\1|"`
         if [[ -z $ricore ]]; then
            ricore=500
         fi

         if [[ $numfor = false && $nocontrol = false ]]; then
            MEM_PER_CPU=$(( (maxcor + ricore + THE_NODES * TASK_PER_NODE * CPU_PER_TASK + $BUFFER - 1 ) \
                          / (THE_NODES * TASK_PER_NODE * CPU_PER_TASK) ))
         elif [[ $nocontrol = false ]]; then
            MEM_PER_CPU=$(( (maxcor + ricore + $BUFFER) ))
         fi
      fi

      MEM_PER_CPU=${MEM_PER_CPU}mb

      #Building the job name from the parent directory names
      PDIR_DOWN=`echo "$i" | sed -rn  "s|(/[^/]*){$P_SAVED}$||p"`
      PDIR_DOWN=`echo "${i#$PDIR_DOWN}"`
      PDIR_DOWN=`echo $PDIR_DOWN | sed -r 's|/|_|g'`
      echo $PDIR_DOWN

      #Writing the individualized values to the run_slurm batch script
      sed "s|DIR_HERE|${i}|" temp_slurm_hold > temp_slurm_dir
      sed "s|MEM_PER_CPU|$MEM_PER_CPU|" temp_slurm_dir > temp_slurm_dir2
      sed "s|THE_JOB_NAME|${THE_JOB_NAME}${PDIR_DOWN}|" temp_slurm_dir2 > run_slurm
      mv run_slurm "${i}"

      if [[ $debug = false ]]; then
         cd "${i}" && sbatch run_slurm
         cd "$WORKDIR"
      fi

      rm temp_slurm_dir*
      if [[ $nocontrol = false ]]; then
         MEM_PER_CPU=''
      fi
   fi

   echo '\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'
done

rm temp_slurm_hold

fi

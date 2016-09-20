# /bin/sh

if [ "$1" = "-h" ]; then
  cat << HERE
USAGE: $0 [-h] DIRECTORIES and/or FILES

The purpose of this script is to back up whatever files are given to
it as argument. These files may be on local or remote servers. 

This script works best in conjunction with the cron task scheduler. 
Simply create a textfile containing the full file pathnames (preceded
by the hostname if file is on a remote server) of whatever needs 
to be backed up. Then schedule a cron job that calls this script
with an argument of \$(cat filename).

The script is written such that filename expansion may be
used to specify many files at once ( *, {...}, etc.)

The script starts by rsyncing all files given to it as argument to BACKUP. 
It then checks to see if the file matches any other versions contained
in the ARCHIVE folder. If it does not, the script puts a copy of the 
file in ARCHIVE with the date that file was last modified as extension. 

For a less storage hungry implementation the script should be rewritten
to use git and just commit changes after every rsync. 
HERE

  exit 0 
fi


BACKUP=
ARCHIVE=

#Expected input is all files specified in ~/.cron_backup 
#if the files are on a remote host then the host name must be specified. 
for i in $@; do

    #Get just the names of the files copied over. The way I do it
    #allows for expansions to be used ( *, {...}, etc. )
    host=''
    entries=''
    nohost=`echo $i | sed -rn "s|[^@]*@[^:]*:(.*)|\1|p"`
    if [ "$nohost" = '' ]; then
        entries=( $(eval ls -d1 $i) )
    else
        host=`echo $i | sed -rn "s|([^@]*@[^:]*):.*|\1|p"`
        entries=( $(ssh $host "ls -d1  $(echo "$nohost")") )
        host=${host}:
    fi

    


    #This will return the name of $entry and all files it contains
    #to all levels
    for entry in ${entries[@]}; do
        path=$(echo $entry | sed -rn "s|(.*)/[^/]*/*$|\1|p")

        #Reproduce the given path to the file, for now there
        #is no host identification so directories may not be
        #unique. 
        mkdir -p $BACKUP/$path
        mkdir -p $ARCHIVE/$path

        #Retrieve this entry
        eval rsync -a $host$entry $BACKUP/$path

        files=`find $BACKUP/$entry  -name '*'`
 
        #for $entry and all files within
        for file in ${files[@]}; do
 
            #The below makes the archive directory if it needs making.
            #find lists the directory before any files it contains
            #so we'll never copy a file to a nonexistent directory. 
            if [ -d $file ]; then 
                if [ ! -d $ARCHIVE/${file#$BACKUP} ]; then
                    mkdir $ARCHIVE/${file#$BACKUP}
                fi
            else
                #for each file in $ARCHIVE with $file as basename and 
                #the last modified date appended, 
                #check if the file in $BACKUP is different. If no $ARCHIVE
                #files match, make a time stamped copy in $ARCHIVE. 
                newCopy=true
                for arch in `find $ARCHIVE -regextype sed \
                    -regex "$ARCHIVE/${file#$BACKUP}\.[0-9]\{4\}-[0-9]\{2\}-[0-9]\{2\}\.[0-9]\{2\}:[0-9]\{2\}:[0-9]\{2\}"`; do
                    check=`diff $file $ARCHIVE/${file#$BACKUP} 2>/dev/null`
 
                    if [ "$check" = '' ]; then
                        newCopy=false
                    fi
                done
                #make a copy timestamped with the last modification date
                if [ "$newCopy" = 'true' ]; then
                    cp $file $ARCHIVE/${file#$BACKUP}.`date +%Y-%m-%d.%H:%M:%S -r $file`
                fi
            fi
        done
    done
done

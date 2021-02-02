# Sets up these python libraries in user's python path

# BASH_SOURCE holds the name of the script being sourced
if [ $(uname) == "Linux" ]; then
    SCRIPTNAME=$(readlink -f $BASH_SOURCE)
elif [ $(uname) == "Darwin" ]; then
    SCRIPTNAME=$(greadlink -f $BASH_SOURCE)
else
    echo "***ERROR: unrecognized host OS"
fi

# % removes the shortest matching pattern
# /* expands to "/<and any other characters>"
#
LIBDIR=${SCRIPTNAME%/*}

KSPCALCDIR=$LIBDIR/kspcalc

NEWPYTHONPATH=$LIBDIR:$KSPCALCDIR

if ! echo $PYTHONPATH | grep -q "$NEWPYTHONPATH"; then
    
    if [ -z $PYTHONPATH ]; then
	export PYTHONPATH=$NEWPYTHONPATH
    else
	export PYTHONPATH=$PYTHONPATH:$NEWPYTHONPATH
    fi

    echo "PYTHONPATH HAS BEEN UPDATED TO $PYTHONPATH"
fi

if ! echo $PATH | grep -q "$NEWPYTHONPATH"; then
    
    if [ -z $PATH ]; then
	export PATH=$NEWPYTHONPATH
    else
	export PATH=$PATH:$NEWPYTHONPATH
    fi

    echo "PATH HAS BEEN UPDATED TO $PATH"
fi

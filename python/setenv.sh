# Sets up these python libraries in user's python path

# BASH_SOURCE holds the name of the script being sourced
SCRIPTNAME=$(readlink -f $BASH_SOURCE)

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
fi

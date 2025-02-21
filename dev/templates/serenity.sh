#!/bin/bash

###########################
# Serenity Main Variables #
###########################

# The main folder of Serenity

if [[ ${BASH_SOURCE[0]} == /* ]] 
then export SERENITY_HOME=$(dirname "${BASH_SOURCE[0]}") 
else export SERENITY_HOME=$PWD/$(dirname "${BASH_SOURCE[0]}")
fi

# Path to  files
export SERENITY_RESOURCES=$SERENITY_HOME/data/

# Path to the executables of Serenity
export SERENITY_BIN=$SERENITY_HOME/bin/

export PATH=$SERENITY_BIN:$PATH

# Add the shared serenity library to the library search path
export LD_LIBRARY_PATH=$SERENITY_HOME/lib/:$LD_LIBRARY_PATH

##########################################################
# The following variables are also respected by Serenity #
##########################################################

# Dynamic memory (only for the largest chunks of data) used by Serenity (in MB), if not set half of
# the total memory of the system is used.
# export SERENITY_MEMORY=4000

# Number of threads used for OpenMP parallelization (if compiled with OpenMP which is the default)
# export OMP_NUM_THREADS=4


##################
# Python Wrapper #
##################

# Add the Python module to the Python path
export PYTHONPATH=$SERENITY_HOME/:$PYTHONPATH

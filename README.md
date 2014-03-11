# GIP and gippy

**GIP** is a high-performance Geospatial Image Processing library written in
C++. 

**gippy** is a python interface to the GIP library and provides a framework for
maintaining spatial data assets.

## gippy Development Note
gippy development utilizes the distutils *setup.py* file for development work
and the generation of SWIG interfaces to the *GIP* library.

## GIP Development Note

For developing GIP, (until integrated into **setup.py**)it is recommended that
you use a python virtual environment and the _develop_ target in the
*Makefile*.  This allows multiple users on the same system to independently
develop the GIP library without collisions.   To maintain consistency with how
virtualenv works, consider adding the following chunks of code to your activate
file. 

1) In deactivate function definition:  

    if [ -n "$_OLD_LD_LIBRARY_PATH" ] ; then
        LD_LIBRARY_PATH="$_OLD_LD_LIBRARY_PATH"
        export LD_LIBRARY_PATH
        unset _OLD_LD_LIBRARY_PATH
    fi

2) And anywhere in the activate script but not in the deactivate function
definition:  

    if [ -n "${LD_LIBRARY_PATH}" ] ; then
       _OLD_LD_LIBRARY_PATH="${LD_LIBRARY_PATH}"
    fi
    LD_LIBRARY_PATH=${VIRTUAL_ENV}/lib/:${LD_LIBRARY_PATH}
    export LD_LIBRARY_PATH
    
    # # for MACOSX
    # if [ -n "${DYLD_LIBRARY_PATH}" ] ; then
    #    _OLD_DYLD_FALLBACK_LIBRARY_PATH="${DYLD_FALLBACK_LIBRARY_PATH}"
    # fi
    # DYLD_LIBRARY_PATH=${VIRTUAL_ENV}/lib/:${DYLD_LIBRARY_PATH}
    # export DYLD_LIBRARY_PATH


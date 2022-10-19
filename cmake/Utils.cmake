# List of likely locations to find header files and libraries

# clear symbols to prevent list bloating after multiple inclusions
set( USER_LIKELY_ROOT_DIR )
set( USER_LIKELY_INCLUDE_DIR )
set( USER_LIKELY_LIBRARY_DIR )
set( USER_ME $ENV{USER} )

# installation root directories
list( APPEND USER_LIKELY_ROOT_DIR 
    /usr
    /usr/local
    /opt
    /opt/local
    /sw
    /sw/local
    /scratch/${USER_ME}
    /scratch/${USER_ME}/local
    /scratch2/${USER_ME}
    /scratch2/${USER_ME}/local
    /scratch3/${USER_ME}
    /scratch3/${USER_ME}/local
    /scratch4/${USER_ME}
    /scratch4/${USER_ME}/local
    $ENV{HOME}/code
)

# corresponding include and lib directories
foreach( dir ${USER_LIKELY_ROOT_DIR} )
    list( APPEND USER_LIKELY_INCLUDE_DIR ${dir}/include )
    list( APPEND USER_LIKELY_LIBRARY_DIR ${dir}/lib ${dir}/bin )
endforeach()

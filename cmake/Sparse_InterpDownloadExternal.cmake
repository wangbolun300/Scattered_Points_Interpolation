include(DownloadProject)

# With CMake 3.8 and above, we can hide warnings about git being in a
# detached head by passing an extra GIT_CONFIG option
if(NOT (${CMAKE_VERSION} VERSION_LESS "3.8.0"))
    set(SPARSE_INTERP_EXTRA_OPTIONS "GIT_CONFIG advice.detachedHead=false")
else()
    set(SPARSE_INTERP_EXTRA_OPTIONS "")
endif()

function(sparse_interp_download_project name)
    download_project(
        PROJ         ${name}
        SOURCE_DIR   ${SPARSE_EXTERNAL}/${name}
        DOWNLOAD_DIR ${SPARSE_EXTERNAL}/.cache/${name}
        QUIET
        ${SPARSE_INTERP_EXTRA_OPTIONS}
        ${ARGN}
    )
endfunction()

################################################################################


# Eigen
function(sparse_interp_download_eigen)
  sparse_interp_download_project(eigen
  GIT_REPOSITORY           https://gitlab.com/libeigen/eigen.git
	GIT_TAG       3.3.7
  )
endfunction()



# libigl for timing and mesh processing
function(sparse_interp_download_libigl)
   sparse_interp_download_project(libigl
   GIT_REPOSITORY https://github.com/libigl/libigl.git
   GIT_TAG        aea868bd1fc64f71afecd2c51e51507a99d8e3e5
  )
endfunction()

# A modern string formatting library
function(sparse_interp_download_fmt)
  sparse_interp_download_project(fmt
    GIT_REPOSITORY https://github.com/fmtlib/fmt.git
    GIT_TAG        6.2.0
  )
endfunction()



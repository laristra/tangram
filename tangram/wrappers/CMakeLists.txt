#[[
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
]]

#-----------------------------------------------------------------------------~#


#-----------------------------------------------------------------------------~#

if (Jali_DIR)
    set(jali_wrapper_headers
        mesh/jali/jali_mesh_wrapper.h
        )
endif (Jali_DIR)

if (FLECSI_INSTALL_DIR)
    set(flecsi_wrapper_headers
        mesh/flecsi/flecsi_mesh_wrapper.h
        )
endif (FLECSI_INSTALL_DIR)

set(wrappers_HEADERS
    ${jali_wrapper_headers}
    ${flecsi_wrapper_headers}
    mesh/AuxMeshTopology.h
    mesh/simple_mesh/simple_mesh_wrapper.h
    PARENT_SCOPE)

if (ENABLE_UNIT_TESTS)
   if (Jali_DIR)
      cinch_add_unit(test_jali_mesh_wrapper
                    SOURCES mesh/jali/test/test_jali_mesh_wrapper.cc
                    LIBRARIES ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
                    INPUTS mesh/jali/test/shotshell.exo
                    mesh/jali/test/shotshell-v.exo
                    POLICY MPI
                    THREADS 1)
   endif (Jali_DIR)

   if (FLECSI_INSTALL_DIR)
      cinch_add_unit(test_flecsi_mesh_wrapper
                     SOURCES mesh/flecsi/test/test_flecsi_mesh_wrapper.cc
                     LIBRARIES ${FLECSI_LIBRARIES})
   endif (FLECSI_INSTALL_DIR)

   # Simple Mesh wrapper
   cinch_add_unit(test_simple_mesh_wrapper
     SOURCES mesh/simple_mesh/test/test_simple_mesh_wrapper.cc
#     LIBRARIES tangram
     POLICY MPI
     THREADS 1)
endif(ENABLE_UNIT_TESTS)


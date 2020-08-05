#[[
This file is part of the Ristra Tangram project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/tangram/blob/master/LICENSE
]]

add_executable(tangram_multidef_check main.cc dummy.cc)
target_link_libraries(tangram_multidef_check tangram)


set(lib lmfit)
set(${lib}_LIBRARY ${lib} PARENT_SCOPE)

set(src_files lmcurve.c lmcurve2.c lmmin.c lminvert.c)
set(inc_files lmcurve.h lmcurve2.h lmmin.h lmstruct.h lmdecls.h lmfit.hpp)

add_library(${lib} ${src_files})

target_compile_definitions(${lib} PRIVATE "LMFIT_EXPORT") # for Windows DLL export&import

set_target_properties(
    ${lib} PROPERTIES
    OUTPUT_NAME ${lib}
    VERSION ${lmfit_VERSION}
    SOVERSION ${lmfit_VERSION_MAJOR})

include(LinkLibMath)
link_libm(${lib})

install(
    TARGETS ${lib}
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION lib
    ARCHIVE DESTINATION lib
    COMPONENT Libraries)
install(
    FILES ${inc_files}
    DESTINATION include
    COMPONENT Headers)

if(WIN32)
    set(CPACK_PACKAGE_VENDOR "Forschungszentrum Jülich GmbH")
    set(CPACK_GENERATOR "ZIP")
    set(CPACK_PACKAGE_DIRECTORY ${CMAKE_SOURCE_DIR})
    include(CPack)
endif()

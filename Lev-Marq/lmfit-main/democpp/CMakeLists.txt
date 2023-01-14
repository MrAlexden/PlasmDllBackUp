set(cppdemos
curve1
nonlin1
surface1)

foreach(app ${cppdemos})
    add_executable(${app}cpp ${app}.cpp)
    target_include_directories(${app}cpp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${lmfit_SOURCE_DIR}/lib)
    target_link_libraries(${app}cpp ${lmfit_LIBRARY})
endforeach()

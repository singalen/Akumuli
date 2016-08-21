# Main executable
add_executable(akumuli_sample
		hello.cpp
)

target_link_libraries(akumuli_sample
    akumuli
    ${Boost_LIBRARIES}
)

include(CppcheckTargets)
add_cppcheck(akumuli_sample UNUSED_FUNCTIONS STYLE POSSIBLE_ERROR FORCE)

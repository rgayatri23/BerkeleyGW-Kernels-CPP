

if(GPP_USE_OPENMP)
    gpp_google_test(gpp-test-openmp
        DEFINITIONS
            GPP_FIXTURE_NAME=gpp_test_openmp
        SOURCES
            ${CMAKE_CURRENT_SOURCE_DIR}/gpp_test.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/gpp_test.cpp
        LINK_LIBRARIES
            gpp-openmp
            gpp-source
            gpp-std-complex
            gpp-compile-options
        )
endif()

if(GPP_USE_OPENMP_TARGET)
    gpp_google_test(gpp-test-openmp-target
        DEFINITIONS
            GPP_FIXTURE_NAME=gpp_test_openmp_target
        SOURCES
            ${CMAKE_CURRENT_SOURCE_DIR}/gpp_test.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/gpp_test.cpp
        LINK_LIBRARIES
            gpp-openmp-target
            gpp-source
            gpp-std-complex
            gpp-compile-options
        )
endif()

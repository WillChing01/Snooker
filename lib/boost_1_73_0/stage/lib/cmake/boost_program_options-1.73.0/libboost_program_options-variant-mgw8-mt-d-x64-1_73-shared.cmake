# Generated by Boost 1.73.0

# address-model=64

if(CMAKE_SIZEOF_VOID_P EQUAL 4)
  _BOOST_SKIPPED("libboost_program_options-mgw8-mt-d-x64-1_73.dll.a" "64 bit, need 32")
  return()
endif()

# layout=versioned

# toolset=mgw8

if(Boost_COMPILER)
  if(NOT Boost_COMPILER STREQUAL "mgw8" AND NOT Boost_COMPILER STREQUAL "-mgw8")
    _BOOST_SKIPPED("libboost_program_options-mgw8-mt-d-x64-1_73.dll.a" "mgw8, Boost_COMPILER=${Boost_COMPILER}")
    return()
  endif()
else()
  if(BOOST_DETECTED_TOOLSET AND NOT BOOST_DETECTED_TOOLSET STREQUAL "mgw8")
    _BOOST_SKIPPED("libboost_program_options-mgw8-mt-d-x64-1_73.dll.a" "mgw8, detected ${BOOST_DETECTED_TOOLSET}, set Boost_COMPILER to override")
    return()
  endif()
endif()

# link=shared

if(DEFINED Boost_USE_STATIC_LIBS)
  if(Boost_USE_STATIC_LIBS)
    _BOOST_SKIPPED("libboost_program_options-mgw8-mt-d-x64-1_73.dll.a" "shared, Boost_USE_STATIC_LIBS=${Boost_USE_STATIC_LIBS}")
    return()
  endif()
else()
  if(WIN32 AND NOT _BOOST_SINGLE_VARIANT)
    _BOOST_SKIPPED("libboost_program_options-mgw8-mt-d-x64-1_73.dll.a" "shared, default on Windows is static, set Boost_USE_STATIC_LIBS=OFF to override")
    return()
  endif()
endif()

# runtime-link=shared

if(Boost_USE_STATIC_RUNTIME)
  _BOOST_SKIPPED("libboost_program_options-mgw8-mt-d-x64-1_73.dll.a" "shared runtime, Boost_USE_STATIC_RUNTIME=${Boost_USE_STATIC_RUNTIME}")
  return()
endif()

# runtime-debugging=on

if(NOT "${Boost_USE_DEBUG_RUNTIME}" STREQUAL "" AND NOT Boost_USE_DEBUG_RUNTIME)
  _BOOST_SKIPPED("libboost_program_options-mgw8-mt-d-x64-1_73.dll.a" "debug runtime, Boost_USE_DEBUG_RUNTIME=${Boost_USE_DEBUG_RUNTIME}")
  return()
endif()

# threading=multi

if(DEFINED Boost_USE_MULTITHREADED AND NOT Boost_USE_MULTITHREADED)
  _BOOST_SKIPPED("libboost_program_options-mgw8-mt-d-x64-1_73.dll.a" "multithreaded, Boost_USE_MULTITHREADED=${Boost_USE_MULTITHREADED}")
  return()
endif()

# variant=debug

if(NOT "${Boost_USE_DEBUG_LIBS}" STREQUAL "" AND NOT Boost_USE_DEBUG_LIBS)
  _BOOST_SKIPPED("libboost_program_options-mgw8-mt-d-x64-1_73.dll.a" "debug, Boost_USE_DEBUG_LIBS=${Boost_USE_DEBUG_LIBS}")
  return()
endif()

if(Boost_VERBOSE OR Boost_DEBUG)
  message(STATUS "  [x] libboost_program_options-mgw8-mt-d-x64-1_73.dll.a")
endif()

# Create imported target Boost::program_options

if(NOT TARGET Boost::program_options)
  add_library(Boost::program_options SHARED IMPORTED)

  set_target_properties(Boost::program_options PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${_BOOST_INCLUDEDIR}"
    INTERFACE_COMPILE_DEFINITIONS "BOOST_ALL_NO_LIB"
  )
endif()

# Target file name: libboost_program_options-mgw8-mt-d-x64-1_73.dll.a

get_target_property(__boost_imploc Boost::program_options IMPORTED_IMPLIB_DEBUG)
if(__boost_imploc)
  message(SEND_ERROR "Target Boost::program_options already has an imported location '${__boost_imploc}', which is being overwritten with '${_BOOST_LIBDIR}/libboost_program_options-mgw8-mt-d-x64-1_73.dll.a'")
endif()
unset(__boost_imploc)

set_property(TARGET Boost::program_options APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)

set_target_properties(Boost::program_options PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG CXX
  IMPORTED_IMPLIB_DEBUG "${_BOOST_LIBDIR}/libboost_program_options-mgw8-mt-d-x64-1_73.dll.a"
  )

set_property(TARGET Boost::program_options APPEND
  PROPERTY INTERFACE_COMPILE_DEFINITIONS "BOOST_PROGRAM_OPTIONS_DYN_LINK"
  )

list(APPEND _BOOST_PROGRAM_OPTIONS_DEPS headers)
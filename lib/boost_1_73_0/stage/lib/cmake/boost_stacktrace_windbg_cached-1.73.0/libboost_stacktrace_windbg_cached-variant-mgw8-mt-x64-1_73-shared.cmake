# Generated by Boost 1.73.0

# address-model=64

if(CMAKE_SIZEOF_VOID_P EQUAL 4)
  _BOOST_SKIPPED("libboost_stacktrace_windbg_cached-mgw8-mt-x64-1_73.dll.a" "64 bit, need 32")
  return()
endif()

# layout=versioned

# toolset=mgw8

if(Boost_COMPILER)
  if(NOT Boost_COMPILER STREQUAL "mgw8" AND NOT Boost_COMPILER STREQUAL "-mgw8")
    _BOOST_SKIPPED("libboost_stacktrace_windbg_cached-mgw8-mt-x64-1_73.dll.a" "mgw8, Boost_COMPILER=${Boost_COMPILER}")
    return()
  endif()
else()
  if(BOOST_DETECTED_TOOLSET AND NOT BOOST_DETECTED_TOOLSET STREQUAL "mgw8")
    _BOOST_SKIPPED("libboost_stacktrace_windbg_cached-mgw8-mt-x64-1_73.dll.a" "mgw8, detected ${BOOST_DETECTED_TOOLSET}, set Boost_COMPILER to override")
    return()
  endif()
endif()

# link=shared

if(DEFINED Boost_USE_STATIC_LIBS)
  if(Boost_USE_STATIC_LIBS)
    _BOOST_SKIPPED("libboost_stacktrace_windbg_cached-mgw8-mt-x64-1_73.dll.a" "shared, Boost_USE_STATIC_LIBS=${Boost_USE_STATIC_LIBS}")
    return()
  endif()
else()
  if(WIN32 AND NOT _BOOST_SINGLE_VARIANT)
    _BOOST_SKIPPED("libboost_stacktrace_windbg_cached-mgw8-mt-x64-1_73.dll.a" "shared, default on Windows is static, set Boost_USE_STATIC_LIBS=OFF to override")
    return()
  endif()
endif()

# runtime-link=shared

if(Boost_USE_STATIC_RUNTIME)
  _BOOST_SKIPPED("libboost_stacktrace_windbg_cached-mgw8-mt-x64-1_73.dll.a" "shared runtime, Boost_USE_STATIC_RUNTIME=${Boost_USE_STATIC_RUNTIME}")
  return()
endif()

# runtime-debugging=off

if(Boost_USE_DEBUG_RUNTIME)
  _BOOST_SKIPPED("libboost_stacktrace_windbg_cached-mgw8-mt-x64-1_73.dll.a" "release runtime, Boost_USE_DEBUG_RUNTIME=${Boost_USE_DEBUG_RUNTIME}")
  return()
endif()

# threading=multi

if(DEFINED Boost_USE_MULTITHREADED AND NOT Boost_USE_MULTITHREADED)
  _BOOST_SKIPPED("libboost_stacktrace_windbg_cached-mgw8-mt-x64-1_73.dll.a" "multithreaded, Boost_USE_MULTITHREADED=${Boost_USE_MULTITHREADED}")
  return()
endif()

# variant=release

if(NOT "${Boost_USE_RELEASE_LIBS}" STREQUAL "" AND NOT Boost_USE_RELEASE_LIBS)
  _BOOST_SKIPPED("libboost_stacktrace_windbg_cached-mgw8-mt-x64-1_73.dll.a" "release, Boost_USE_RELEASE_LIBS=${Boost_USE_RELEASE_LIBS}")
  return()
endif()

if(Boost_VERBOSE OR Boost_DEBUG)
  message(STATUS "  [x] libboost_stacktrace_windbg_cached-mgw8-mt-x64-1_73.dll.a")
endif()

# Create imported target Boost::stacktrace_windbg_cached

if(NOT TARGET Boost::stacktrace_windbg_cached)
  add_library(Boost::stacktrace_windbg_cached SHARED IMPORTED)

  set_target_properties(Boost::stacktrace_windbg_cached PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${_BOOST_INCLUDEDIR}"
    INTERFACE_COMPILE_DEFINITIONS "BOOST_ALL_NO_LIB"
  )
endif()

# Target file name: libboost_stacktrace_windbg_cached-mgw8-mt-x64-1_73.dll.a

get_target_property(__boost_imploc Boost::stacktrace_windbg_cached IMPORTED_IMPLIB_RELEASE)
if(__boost_imploc)
  message(SEND_ERROR "Target Boost::stacktrace_windbg_cached already has an imported location '${__boost_imploc}', which is being overwritten with '${_BOOST_LIBDIR}/libboost_stacktrace_windbg_cached-mgw8-mt-x64-1_73.dll.a'")
endif()
unset(__boost_imploc)

set_property(TARGET Boost::stacktrace_windbg_cached APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)

set_target_properties(Boost::stacktrace_windbg_cached PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE CXX
  IMPORTED_IMPLIB_RELEASE "${_BOOST_LIBDIR}/libboost_stacktrace_windbg_cached-mgw8-mt-x64-1_73.dll.a"
  )

set_target_properties(Boost::stacktrace_windbg_cached PROPERTIES
  MAP_IMPORTED_CONFIG_MINSIZEREL Release
  MAP_IMPORTED_CONFIG_RELWITHDEBINFO Release
  )

set_property(TARGET Boost::stacktrace_windbg_cached APPEND
  PROPERTY INTERFACE_COMPILE_DEFINITIONS "BOOST_STACKTRACE_WINDBG_CACHED_DYN_LINK"
  )

list(APPEND _BOOST_STACKTRACE_WINDBG_CACHED_DEPS headers)

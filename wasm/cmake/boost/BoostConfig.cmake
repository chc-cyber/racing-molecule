# Our RDKit build needs Boost headers only. No host Boost binaries may be linked
# into a WebAssembly target. All non-header Boost components are disabled.
set(Boost_FOUND TRUE)
set(Boost_VERSION_STRING "1.85.0")
set(Boost_INCLUDE_DIRS "${MOLECULAR_CAR_BOOST_DIR}")
set(Boost_LIBRARIES Boost::headers)
if(Boost_FIND_COMPONENTS)
  message(FATAL_ERROR "This browser build permits header-only Boost only")
endif()
if(NOT TARGET Boost::headers)
  add_library(Boost::headers INTERFACE IMPORTED GLOBAL)
  set_target_properties(Boost::headers PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${Boost_INCLUDE_DIRS}")
  add_library(Boost::boost ALIAS Boost::headers)
endif()

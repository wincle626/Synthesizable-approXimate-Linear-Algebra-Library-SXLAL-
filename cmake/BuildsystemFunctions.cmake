# ------------------------------------------------------------------------------
#! @brief Search for a specific value in a list and return a boolean result
#!        indicating if the entry was found.
#! @param[in]  LIST A list of entries that should be searched.
#! @param[in]  ITEM The value to look for in the list 'LIST'
#! @param[out] OUTPUT_VAR Boolean value indicating whether ITEM was found in
#!             LIST.
function(mrg_list_contains LIST ITEM OUTPUT_VAR)
  list(FIND LIST "${ITEM}" INDEX)
  if(INDEX LESS 0)
    set(RESULT 0)
  else()
    set(RESULT 1)
  endif()
  set(${OUTPUT_VAR} "${RESULT}" PARENT_SCOPE)
endfunction()

# ----------------------------------------------------------------------------
#! @brief Initialize CMAKE_BUILD_TYPE cache variable and constrain its contents
#!
#! @param[in] DEFAULT_BUILD_TYPE [Optional]
#!            Build type to use if user hasn't specified. Defaults to 'Release'.
function(mrg_setup_cmake_build_type DEFAULT_BUILD_TYPE)
  if (NOT DEFAULT_BUILD_TYPE)
    set(DEFAULT_BUILD_TYPE Release)
  endif ()

  # The set of valid build types
  if(NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_CONFIGURATION_TYPES
      "Debug" "Release" "MinSizeRel" "RelWithDebInfo"
    )
  endif()

  # If user specifies the build type, use theirs, otherwise use
  # the value in DEFAULT_BUILD_TYPE.
  if(CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE STRING
      "Choose the type of build.")
  else()
    mrg_list_contains("${CMAKE_CONFIGURATION_TYPES}"
      "${DEFAULT_BUILD_TYPE}" DEFAULT_BUILD_TYPE_VALID)
    if(NOT DEFAULT_BUILD_TYPE_VALID)
      message(FATAL_ERROR "mrg_setup_cmake_build_type: DEFAULT_BUILD_TYPE is "
        "set to '${DEFAULT_BUILD_TYPE}', but it's not in the list of valid "
        "strings: '${CMAKE_CONFIGURATION_TYPES}'.")
    endif()
    set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE STRING
      "Choose the type of build." FORCE)
  endif()

  # Constrain build type to take only valid strings.
  set_property(CACHE CMAKE_BUILD_TYPE
    PROPERTY STRINGS "${CMAKE_CONFIGURATION_TYPES}"
  )
endfunction()


# ----------------------------------------------------------------------------
#! @brief Standard boilerplate for configuring executable RPATHs
#!
#! Causes installed executables to look for libs in the 'lib' dir adjacent to
#! the dir containing the executable. This should make installed packages
#! relocatable.
macro(mrg_setup_executable_rpath)

  # Use @rpath library install names on OSX
  # Also set by policy CMP00042
  set(CMAKE_MACOSX_RPATH TRUE)

  # Use (don't skip) RPATH for binaries in the build tree
  set(CMAKE_SKIP_BUILD_RPATH FALSE)

  # Don't append external paths to the RPATH
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)

  # Don't use the install RPATH for a normal build.
  # Only use it when the binaries are installed.
  set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

  # Set RPATH value for installed executables. This should point
  # at the directory containing the libraries. Relative paths used to allow
  # installed files to be relocated.
  if(APPLE)
    set(CMAKE_INSTALL_RPATH "@executable_path/../lib")
  elseif(UNIX)
    set(CMAKE_INSTALL_RPATH "\$ORIGIN/../lib")
  endif()
endmacro()

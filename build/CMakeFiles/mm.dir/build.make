# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/yunwu/gitworkspace/Xprecision_Linear_Algebra_Library

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/yunwu/gitworkspace/Xprecision_Linear_Algebra_Library/build

# Include any dependencies generated for this target.
include CMakeFiles/mm.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/mm.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/mm.dir/flags.make

CMakeFiles/mm.dir/examples/c/mm/mm.cpp.o: CMakeFiles/mm.dir/flags.make
CMakeFiles/mm.dir/examples/c/mm/mm.cpp.o: ../examples/c/mm/mm.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/yunwu/gitworkspace/Xprecision_Linear_Algebra_Library/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/mm.dir/examples/c/mm/mm.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mm.dir/examples/c/mm/mm.cpp.o -c /home/yunwu/gitworkspace/Xprecision_Linear_Algebra_Library/examples/c/mm/mm.cpp

CMakeFiles/mm.dir/examples/c/mm/mm.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mm.dir/examples/c/mm/mm.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/yunwu/gitworkspace/Xprecision_Linear_Algebra_Library/examples/c/mm/mm.cpp > CMakeFiles/mm.dir/examples/c/mm/mm.cpp.i

CMakeFiles/mm.dir/examples/c/mm/mm.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mm.dir/examples/c/mm/mm.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/yunwu/gitworkspace/Xprecision_Linear_Algebra_Library/examples/c/mm/mm.cpp -o CMakeFiles/mm.dir/examples/c/mm/mm.cpp.s

CMakeFiles/mm.dir/examples/c/mm/mm.cpp.o.requires:

.PHONY : CMakeFiles/mm.dir/examples/c/mm/mm.cpp.o.requires

CMakeFiles/mm.dir/examples/c/mm/mm.cpp.o.provides: CMakeFiles/mm.dir/examples/c/mm/mm.cpp.o.requires
	$(MAKE) -f CMakeFiles/mm.dir/build.make CMakeFiles/mm.dir/examples/c/mm/mm.cpp.o.provides.build
.PHONY : CMakeFiles/mm.dir/examples/c/mm/mm.cpp.o.provides

CMakeFiles/mm.dir/examples/c/mm/mm.cpp.o.provides.build: CMakeFiles/mm.dir/examples/c/mm/mm.cpp.o


# Object files for target mm
mm_OBJECTS = \
"CMakeFiles/mm.dir/examples/c/mm/mm.cpp.o"

# External object files for target mm
mm_EXTERNAL_OBJECTS =

mm: CMakeFiles/mm.dir/examples/c/mm/mm.cpp.o
mm: CMakeFiles/mm.dir/build.make
mm: CMakeFiles/mm.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/yunwu/gitworkspace/Xprecision_Linear_Algebra_Library/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable mm"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/mm.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/mm.dir/build: mm

.PHONY : CMakeFiles/mm.dir/build

CMakeFiles/mm.dir/requires: CMakeFiles/mm.dir/examples/c/mm/mm.cpp.o.requires

.PHONY : CMakeFiles/mm.dir/requires

CMakeFiles/mm.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/mm.dir/cmake_clean.cmake
.PHONY : CMakeFiles/mm.dir/clean

CMakeFiles/mm.dir/depend:
	cd /home/yunwu/gitworkspace/Xprecision_Linear_Algebra_Library/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/yunwu/gitworkspace/Xprecision_Linear_Algebra_Library /home/yunwu/gitworkspace/Xprecision_Linear_Algebra_Library /home/yunwu/gitworkspace/Xprecision_Linear_Algebra_Library/build /home/yunwu/gitworkspace/Xprecision_Linear_Algebra_Library/build /home/yunwu/gitworkspace/Xprecision_Linear_Algebra_Library/build/CMakeFiles/mm.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/mm.dir/depend

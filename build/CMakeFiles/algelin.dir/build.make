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
CMAKE_SOURCE_DIR = /home/caio/Documents/numerico

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/caio/Documents/numerico/build

# Include any dependencies generated for this target.
include CMakeFiles/algelin.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/algelin.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/algelin.dir/flags.make

CMakeFiles/algelin.dir/include/algelin.cpp.o: CMakeFiles/algelin.dir/flags.make
CMakeFiles/algelin.dir/include/algelin.cpp.o: ../include/algelin.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/caio/Documents/numerico/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/algelin.dir/include/algelin.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/algelin.dir/include/algelin.cpp.o -c /home/caio/Documents/numerico/include/algelin.cpp

CMakeFiles/algelin.dir/include/algelin.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/algelin.dir/include/algelin.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/caio/Documents/numerico/include/algelin.cpp > CMakeFiles/algelin.dir/include/algelin.cpp.i

CMakeFiles/algelin.dir/include/algelin.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/algelin.dir/include/algelin.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/caio/Documents/numerico/include/algelin.cpp -o CMakeFiles/algelin.dir/include/algelin.cpp.s

CMakeFiles/algelin.dir/include/algelin.cpp.o.requires:

.PHONY : CMakeFiles/algelin.dir/include/algelin.cpp.o.requires

CMakeFiles/algelin.dir/include/algelin.cpp.o.provides: CMakeFiles/algelin.dir/include/algelin.cpp.o.requires
	$(MAKE) -f CMakeFiles/algelin.dir/build.make CMakeFiles/algelin.dir/include/algelin.cpp.o.provides.build
.PHONY : CMakeFiles/algelin.dir/include/algelin.cpp.o.provides

CMakeFiles/algelin.dir/include/algelin.cpp.o.provides.build: CMakeFiles/algelin.dir/include/algelin.cpp.o


# Object files for target algelin
algelin_OBJECTS = \
"CMakeFiles/algelin.dir/include/algelin.cpp.o"

# External object files for target algelin
algelin_EXTERNAL_OBJECTS =

libalgelin.a: CMakeFiles/algelin.dir/include/algelin.cpp.o
libalgelin.a: CMakeFiles/algelin.dir/build.make
libalgelin.a: CMakeFiles/algelin.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/caio/Documents/numerico/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libalgelin.a"
	$(CMAKE_COMMAND) -P CMakeFiles/algelin.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/algelin.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/algelin.dir/build: libalgelin.a

.PHONY : CMakeFiles/algelin.dir/build

CMakeFiles/algelin.dir/requires: CMakeFiles/algelin.dir/include/algelin.cpp.o.requires

.PHONY : CMakeFiles/algelin.dir/requires

CMakeFiles/algelin.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/algelin.dir/cmake_clean.cmake
.PHONY : CMakeFiles/algelin.dir/clean

CMakeFiles/algelin.dir/depend:
	cd /home/caio/Documents/numerico/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/caio/Documents/numerico /home/caio/Documents/numerico /home/caio/Documents/numerico/build /home/caio/Documents/numerico/build /home/caio/Documents/numerico/build/CMakeFiles/algelin.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/algelin.dir/depend

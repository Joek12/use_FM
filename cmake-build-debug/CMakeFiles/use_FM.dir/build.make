# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.14

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
CMAKE_COMMAND = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake

# The command to remove a file.
RM = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/josephkang/CLionProjects/use_FM

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/josephkang/CLionProjects/use_FM/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/use_FM.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/use_FM.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/use_FM.dir/flags.make

CMakeFiles/use_FM.dir/src/main.cpp.o: CMakeFiles/use_FM.dir/flags.make
CMakeFiles/use_FM.dir/src/main.cpp.o: ../src/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/josephkang/CLionProjects/use_FM/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/use_FM.dir/src/main.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/use_FM.dir/src/main.cpp.o -c /Users/josephkang/CLionProjects/use_FM/src/main.cpp

CMakeFiles/use_FM.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/use_FM.dir/src/main.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/josephkang/CLionProjects/use_FM/src/main.cpp > CMakeFiles/use_FM.dir/src/main.cpp.i

CMakeFiles/use_FM.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/use_FM.dir/src/main.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/josephkang/CLionProjects/use_FM/src/main.cpp -o CMakeFiles/use_FM.dir/src/main.cpp.s

CMakeFiles/use_FM.dir/src/fastaReader.cpp.o: CMakeFiles/use_FM.dir/flags.make
CMakeFiles/use_FM.dir/src/fastaReader.cpp.o: ../src/fastaReader.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/josephkang/CLionProjects/use_FM/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/use_FM.dir/src/fastaReader.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/use_FM.dir/src/fastaReader.cpp.o -c /Users/josephkang/CLionProjects/use_FM/src/fastaReader.cpp

CMakeFiles/use_FM.dir/src/fastaReader.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/use_FM.dir/src/fastaReader.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/josephkang/CLionProjects/use_FM/src/fastaReader.cpp > CMakeFiles/use_FM.dir/src/fastaReader.cpp.i

CMakeFiles/use_FM.dir/src/fastaReader.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/use_FM.dir/src/fastaReader.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/josephkang/CLionProjects/use_FM/src/fastaReader.cpp -o CMakeFiles/use_FM.dir/src/fastaReader.cpp.s

CMakeFiles/use_FM.dir/src/spammodule.cpp.o: CMakeFiles/use_FM.dir/flags.make
CMakeFiles/use_FM.dir/src/spammodule.cpp.o: ../src/spammodule.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/josephkang/CLionProjects/use_FM/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/use_FM.dir/src/spammodule.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/use_FM.dir/src/spammodule.cpp.o -c /Users/josephkang/CLionProjects/use_FM/src/spammodule.cpp

CMakeFiles/use_FM.dir/src/spammodule.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/use_FM.dir/src/spammodule.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/josephkang/CLionProjects/use_FM/src/spammodule.cpp > CMakeFiles/use_FM.dir/src/spammodule.cpp.i

CMakeFiles/use_FM.dir/src/spammodule.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/use_FM.dir/src/spammodule.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/josephkang/CLionProjects/use_FM/src/spammodule.cpp -o CMakeFiles/use_FM.dir/src/spammodule.cpp.s

CMakeFiles/use_FM.dir/src/pbar.cpp.o: CMakeFiles/use_FM.dir/flags.make
CMakeFiles/use_FM.dir/src/pbar.cpp.o: ../src/pbar.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/josephkang/CLionProjects/use_FM/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/use_FM.dir/src/pbar.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/use_FM.dir/src/pbar.cpp.o -c /Users/josephkang/CLionProjects/use_FM/src/pbar.cpp

CMakeFiles/use_FM.dir/src/pbar.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/use_FM.dir/src/pbar.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/josephkang/CLionProjects/use_FM/src/pbar.cpp > CMakeFiles/use_FM.dir/src/pbar.cpp.i

CMakeFiles/use_FM.dir/src/pbar.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/use_FM.dir/src/pbar.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/josephkang/CLionProjects/use_FM/src/pbar.cpp -o CMakeFiles/use_FM.dir/src/pbar.cpp.s

CMakeFiles/use_FM.dir/src/primer_design.cpp.o: CMakeFiles/use_FM.dir/flags.make
CMakeFiles/use_FM.dir/src/primer_design.cpp.o: ../src/primer_design.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/josephkang/CLionProjects/use_FM/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/use_FM.dir/src/primer_design.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/use_FM.dir/src/primer_design.cpp.o -c /Users/josephkang/CLionProjects/use_FM/src/primer_design.cpp

CMakeFiles/use_FM.dir/src/primer_design.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/use_FM.dir/src/primer_design.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/josephkang/CLionProjects/use_FM/src/primer_design.cpp > CMakeFiles/use_FM.dir/src/primer_design.cpp.i

CMakeFiles/use_FM.dir/src/primer_design.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/use_FM.dir/src/primer_design.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/josephkang/CLionProjects/use_FM/src/primer_design.cpp -o CMakeFiles/use_FM.dir/src/primer_design.cpp.s

# Object files for target use_FM
use_FM_OBJECTS = \
"CMakeFiles/use_FM.dir/src/main.cpp.o" \
"CMakeFiles/use_FM.dir/src/fastaReader.cpp.o" \
"CMakeFiles/use_FM.dir/src/spammodule.cpp.o" \
"CMakeFiles/use_FM.dir/src/pbar.cpp.o" \
"CMakeFiles/use_FM.dir/src/primer_design.cpp.o"

# External object files for target use_FM
use_FM_EXTERNAL_OBJECTS =

use_FM: CMakeFiles/use_FM.dir/src/main.cpp.o
use_FM: CMakeFiles/use_FM.dir/src/fastaReader.cpp.o
use_FM: CMakeFiles/use_FM.dir/src/spammodule.cpp.o
use_FM: CMakeFiles/use_FM.dir/src/pbar.cpp.o
use_FM: CMakeFiles/use_FM.dir/src/primer_design.cpp.o
use_FM: CMakeFiles/use_FM.dir/build.make
use_FM: CMakeFiles/use_FM.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/josephkang/CLionProjects/use_FM/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking CXX executable use_FM"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/use_FM.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/use_FM.dir/build: use_FM

.PHONY : CMakeFiles/use_FM.dir/build

CMakeFiles/use_FM.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/use_FM.dir/cmake_clean.cmake
.PHONY : CMakeFiles/use_FM.dir/clean

CMakeFiles/use_FM.dir/depend:
	cd /Users/josephkang/CLionProjects/use_FM/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/josephkang/CLionProjects/use_FM /Users/josephkang/CLionProjects/use_FM /Users/josephkang/CLionProjects/use_FM/cmake-build-debug /Users/josephkang/CLionProjects/use_FM/cmake-build-debug /Users/josephkang/CLionProjects/use_FM/cmake-build-debug/CMakeFiles/use_FM.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/use_FM.dir/depend


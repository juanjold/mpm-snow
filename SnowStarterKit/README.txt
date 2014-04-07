Contained in this folder is everything you'll need to get started on the smoke simulation assignment.

A working project is included for Visual Studio 2010. If you're using a different platform, just copy the source files from the SourceFiles folder in the Visual Studio project. Almost everything is cross-platform without modifications. I originally modified it to be cross-platform in order to work on it in Xcode on a Mac (including removing Boost and DevIL dependencies). You might need to change the OpenGL includes in open_gl_headers.h depending on your setup. The only dependency that you might not already have on your computer is GLUT. In case you don't have it, I've included the header file and Windows binaries in the GLUT folder.

This fluid simulation is follows the semi-Lagrangian approach described in Robert Bridson and Matthias Muller-Fischer's fluid simulation course notes, which are included in this package (fluids_notes.pdf) and can also be found at http://www.cs.ubc.ca/~rbridson/fluidsimulation/fluids_notes.pdf . That document is key to completing this assignment.

You will find "TODO" statements throughout the code where you should add the required functionality. The heart of the fluid simulation is the MacGrid class. The bulk of your code will be written in mac_grid.cpp. If you want to extend the framework with additional features and make it clean and robust, feel free to add and modify code elsewhere as well.




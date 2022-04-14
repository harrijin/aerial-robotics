# navtoolbox

The navtoolbox library is an Eigen implementation of the GSS navtoolbox functions. 

To incorporate the navtoolbox library within your project, clone navtoolbox
into the same directory as your top-level CMakeLists.txt file:
```
git clone git@gitlab.com:radionavlab/core/navtoolbox.git
```
Then add this line to the top-level CMakeLists.txt file:
```
add_subdirectory(navtoolbox)
```
Finally, add the library where needed in the top-level or (more likely)
lower-level CMakeLists.txt files.  For example,
```
target_link_libraries(${TARGET} navtoolbox Eigen3::Eigen)
```

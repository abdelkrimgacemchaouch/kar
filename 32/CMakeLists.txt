cmake_minimum_required(VERSION 2.8)
enable_language (C Fortran)
file(GLOB g 
"src/*.f" "src/*.inc" )
add_executable(ben ${g})
target_link_libraries(ben m) #pour les bibiotheque mathematique

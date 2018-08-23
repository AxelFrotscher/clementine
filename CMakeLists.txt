#Set Project Name
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(clementine)

set(CMAKE_CXX_FLAGS "-DCMAKE_BUILD_TYPE=Release")

#Find ROOT installation
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#Locate ROOT package
find_package(ROOT REQUIRED COMPONENTS RIO Net)

#Define useful functions
include(${ROOT_USE_FILE})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include
        /home/afrotscher/soft/anaroot/anarootv4.5.36/include
        ${CMAKE_CURRENT_SOURCE_DIR}/include/cutclasses
        ${CMAKE_CURRENT_SOURCE_DIR}/include/corrections
        ${CMAKE_CURRENT_SOURCE_DIR}/include/helper)

add_executable(${PROJECT_NAME} main.cpp
      ./src/MakeAllTree_78Ni.cpp
      ./src/histograms.cpp
        src/helper/treereader.cpp
        src/helper/histogram_cuts.cpp
        src/cutclasses/triggercut.cpp
        src/cutclasses/ppaccut.cpp
        src/cutclasses/chargestatecut.cpp
        src/cutclasses/plasticcut.cpp
        src/corrections/higherorder.cpp
        src/cutclasses/ICcut.cpp
        src/cutclasses/targetcut.cpp
        src/PID/pid.cpp
        src/helper/progress.cpp)
target_link_libraries(${PROJECT_NAME} ${ROOT_LIBRARIES} -lThread
-L/home/afrotscher/soft/anaroot/anarootv4.5.36/lib
-lanacore -lanabrips -lanadali -lanaroot -lanaminos -lXMLParser
-L/home/afrotscher/soft/minos/minosv1.1.1/install/pentium-linux-release/lib -lminos-fem
-lmordicus-utl -lUtilities ${ROOT_LIBRARIES})

install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION bin)


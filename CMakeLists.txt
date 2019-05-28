#Set Project Name
cmake_minimum_required(VERSION 3.6 FATAL_ERROR)
project(clementine)

set(CMAKE_CXX_STANDARD 17)
set(ROOTSYS "/ldas-ikda/soft/src/cern/root/install/6.16.00/")
message(STATUS "ROOT: ${ROOTSYS}")
#Find ROOT installation
#set($ENV{ROOTSYS} /ldas-ikda/soft/src/cern/root/install/6.16.00/)
list(APPEND CMAKE_PREFIX_PATH ${ROOTSYS})

#Locate ROOT package
find_package(ROOT REQUIRED COMPONENTS RIO Net)

#Define useful functions
include(${ROOT_USE_FILE})
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    /home/afrotscher/soft/anaroot/4.5.36/include    # anaroot installation
    ${CMAKE_CURRENT_SOURCE_DIR}/include/cutclasses
    ${CMAKE_CURRENT_SOURCE_DIR}/include/corrections
    ${CMAKE_CURRENT_SOURCE_DIR}/include/helper)

add_executable(${PROJECT_NAME} main.cpp
    src/MakeAllTree_78Ni.cpp
    src/histograms.cpp
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
    src/helper/progress.cpp
    src/helper/txtwriter.cpp
    src/helper/zdssetting.cpp
    src/helper/minos.cpp src/corrections/minosdrift.cpp)

target_link_libraries(${PROJECT_NAME} ${ROOT_LIBRARIES} -lThread -lMinuit -lMinuit2
   -L/home/afrotscher/soft/anaroot/4.5.36/lib -lanaroot -lanacore -lanabrips -lanadali -lanaminos
        -lXMLParser
   -L/home/afrotscher/soft/minos/liboffline -lminos
   -L/home/afrotscher/soft/minos/minosv1.1.1/install/pentium-linux-release/lib
                                          -lmordicus-utl -lUtilities -lminos-fem
        -lstdc++fs
   ${ROOT_LIBRARIES})

if(CMAKE_BUILD_TYPE MATCHES Debug)
    target_compile_definitions(${PROJECT_NAME} PUBLIC DEBUG=1)
elseif(CMAKE_BUILD_TYPE MATCHES Release)
    target_compile_definitions(${PROJECT_NAME} PUBLIC DEBUG=0)
endif()

install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION bin)
cmake_minimum_required(VERSION 3.10)
project(PDS_Analysis)

# ROOT
find_package(ROOT REQUIRED COMPONENTS Core Tree RIO Hist)
include(${ROOT_USE_FILE})

set(CMAKE_CXX_STANDARD 17)

# Fontes do seu projeto
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)   # implementação da classe
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh) # headers

# Gera o dicionário ROOT (G__my_data.cxx)
ROOT_GENERATE_DICTIONARY(G__my_data ${headers} LINKDEF ${PROJECT_SOURCE_DIR}/include/LinkDef.h)

# Cria biblioteca compartilhada única: implementação + dicionário
add_library(my_data SHARED ${sources} G__my_data.cxx)
target_include_directories(my_data PRIVATE ${PROJECT_SOURCE_DIR}/include ${ROOT_INCLUDE_DIRS})
target_link_libraries(my_data PRIVATE ${ROOT_LIBRARIES})

#------------------------------------------------------------------------------

# Executável run
add_executable(run main.cpp)
target_include_directories(run PRIVATE ${PROJECT_SOURCE_DIR}/include ${ROOT_INCLUDE_DIRS})
target_link_libraries(run PRIVATE ${ROOT_LIBRARIES} my_data)

# Executável flux
add_executable(flux flux.cpp)
target_include_directories(flux PRIVATE ${PROJECT_SOURCE_DIR}/include ${ROOT_INCLUDE_DIRS})
target_link_libraries(flux PRIVATE ${ROOT_LIBRARIES} my_data)

# Executável flash
add_executable(flash flash.cpp)
target_include_directories(flash PRIVATE ${PROJECT_SOURCE_DIR}/include ${ROOT_INCLUDE_DIRS})
target_link_libraries(flash PRIVATE ${ROOT_LIBRARIES} my_data)

# Executável flash
add_executable(filter_flash filter_flash.cpp)
target_include_directories(filter_flash PRIVATE ${PROJECT_SOURCE_DIR}/include ${ROOT_INCLUDE_DIRS})
target_link_libraries(filter_flash PRIVATE ${ROOT_LIBRARIES} my_data)
#..........................................................................
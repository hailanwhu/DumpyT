cmake_minimum_required(VERSION 3.20)
project(Dumpy)

set(CMAKE_CXX_STANDARD 23)

set(CMAKE_CXX_FLAGS "-O0 -g")
#set(CMAKE_BUILD_TYPE Release)
#set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
find_package(Threads REQUIRED)

include_directories(/research/local/hai/boost/boost_install/include)

add_executable(Dumpy src/main.cpp include/DataStructures/TimeSeries.h include/DataStructures/PqItemSeries.h src/Utils/FileUtil.cpp include/Utils/FileUtil.h src/Utils/TimeSeriesUtil.cpp include/Utils/TimeSeriesUtil.h src/Utils/SaxUtil.cpp include/Utils/SaxUtil.h src/Utils/MathUtil.cpp include/Utils/MathUtil.h src/IndexConstruction/GraphConstruction.cpp include/DataStructures/GraphConstruction.h src/PqItemSeries.cpp include/Const.h src/IndexConstruction/FullAryTreeNode.cpp include/DataStructures/FullAryTreeNode.h src/SearchEngine/DumpyInMemorySearch.cpp include/Searchers/DumpyInMemorySearch.h src/Utils/INIReader.cpp include/Utils/INIReader.h src/Utils/ini.c src/Const.cpp src/IndexConstruction/DumpyNode.cpp include/DataStructures/DumpyNode.h src/SearchEngine/DumpySearcher.cpp include/Searchers/DumpySearcher.h src/IndexConstruction/DumpyFuzzy.cpp)

target_link_libraries(Dumpy Threads::Threads libboost_serialization.so.1.76.0)
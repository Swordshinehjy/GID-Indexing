cmake_minimum_required(VERSION 3.21)

project(GIDIndexing VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

if(MSVC)
  add_compile_options(/utf-8 /W3 /permissive- /Zc:__cplusplus /O2)
  add_compile_definitions(_USE_MATH_DEFINES NOMINMAX wxMSVC_VERSION_AUTO)
endif()

# Static runtime: the artifact is a single exe with no external DLL dependencies
set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")

set(VCPKG_STATIC_DIR "D:/vcpkg/installed/x64-windows-static")
# The static triplet must be used: falling back to x64-windows (dynamic) would make the artifact depend on zlib1.dll / tiff.dll,
# which contradicts the "single-exe distribution" goal, so fail with an error when it is missing rather than silently degrading.
if(NOT EXISTS "${VCPKG_STATIC_DIR}")
  message(FATAL_ERROR
    "vcpkg static triplet not found: ${VCPKG_STATIC_DIR}\n"
    "Please run vcpkg install ... --triplet x64-windows-static first, or modify VCPKG_STATIC_DIR in this file.")
endif()

# wxWidgets 3.2 static build (CONFIG mode)
list(APPEND CMAKE_PREFIX_PATH "D:/wxWidgets/install")
find_package(OpenGL REQUIRED COMPONENTS OpenGL)
find_package(wxWidgets CONFIG REQUIRED)

# ---------------------------------------------------------------- core (no wx dependency)
file(GLOB CORE_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/core/*.cpp)

add_library(gid_core STATIC ${CORE_SRC})
target_include_directories(gid_core PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/src
  ${CMAKE_CURRENT_SOURCE_DIR}/src/core
  ${VCPKG_STATIC_DIR}/include)

# xlsx reading (vcpkg static libraries)
find_library(XLSXIO_READ_LIBRARY NAMES xlsxio_read PATHS "${VCPKG_STATIC_DIR}/lib" NO_DEFAULT_PATH)
find_library(MINIZIP_LIBRARY NAMES minizips minizip PATHS "${VCPKG_STATIC_DIR}/lib" NO_DEFAULT_PATH)
find_library(EXPAT_LIBRARY NAMES libexpatMT libexpat PATHS "${VCPKG_STATIC_DIR}/lib" NO_DEFAULT_PATH)
find_library(ICONV_LIBRARY NAMES iconv PATHS "${VCPKG_STATIC_DIR}/lib" NO_DEFAULT_PATH)
find_library(CHARSET_LIBRARY NAMES charset PATHS "${VCPKG_STATIC_DIR}/lib" NO_DEFAULT_PATH)
find_library(ZLIB_LIBRARY NAMES zlib PATHS "${VCPKG_STATIC_DIR}/lib" NO_DEFAULT_PATH)
find_library(BZ2_LIBRARY NAMES bz2 PATHS "${VCPKG_STATIC_DIR}/lib" NO_DEFAULT_PATH)

target_link_libraries(gid_core PUBLIC
  "${XLSXIO_READ_LIBRARY}" "${MINIZIP_LIBRARY}" "${EXPAT_LIBRARY}"
  "${ICONV_LIBRARY}" "${CHARSET_LIBRARY}" "${ZLIB_LIBRARY}" "${BZ2_LIBRARY}")

# ---------------------------------------------------------------- ui (depends on wx)
file(GLOB UI_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/ui/*.cpp)

add_library(gid_ui STATIC ${UI_SRC})
target_include_directories(gid_ui PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/src
  ${CMAKE_CURRENT_SOURCE_DIR}/src/ui)
target_link_libraries(gid_ui PUBLIC gid_core wxWidgets::wxWidgets)

# ---------------------------------------------------------------- executables
if(WIN32)
  enable_language(RC)
  set(APP_RC ${CMAKE_CURRENT_SOURCE_DIR}/app.rc)
endif()

add_executable(GIDIndexing WIN32 src/App.cpp ${APP_RC})
target_link_libraries(GIDIndexing PRIVATE gid_ui)

# Merge the static libraries into a single exe
if(MSVC)
  target_link_options(GIDIndexing PRIVATE /WHOLEARCHIVE:gid_ui /WHOLEARCHIVE:gid_core)
else()
  target_link_options(GIDIndexing PRIVATE
    -Wl,--whole-archive $<TARGET_FILE:gid_ui> $<TARGET_FILE:gid_core> -Wl,--no-whole-archive)
endif()

set_target_properties(GIDIndexing PROPERTIES
  OUTPUT_NAME "GIDIndexing"
  WIN32_EXECUTABLE TRUE
  VS_DPI_AWARE "PerMonitor")

# ---------------------------------------------------------------- core algorithm self-test (optional)
option(GID_BUILD_SELFTEST "Build core self test" ON)
if(GID_BUILD_SELFTEST)
  add_executable(GIDSelfTest tools/SelfTest.cpp)
  target_link_libraries(GIDSelfTest PRIVATE gid_core)
  add_executable(GIDDipProbe tools/DipProbe.cpp)
  target_link_libraries(GIDDipProbe PRIVATE gid_core)
  add_executable(GIDFp tools/_Fp.cpp)
  target_link_libraries(GIDFp PRIVATE gid_core)
  add_executable(GIDFpCsv tools/FpCubicCsv.cpp)
  target_link_libraries(GIDFpCsv PRIVATE gid_core)
  add_executable(GIDNoSpec tools/NoSpecProbe.cpp)
  target_link_libraries(GIDNoSpec PRIVATE gid_core)
  add_executable(GIDNoSpec2 tools/NoSpecTwoStepsProbe.cpp)
  target_link_libraries(GIDNoSpec2 PRIVATE gid_core)
  add_executable(GIDSpecBench tools/SpecBench.cpp)
  target_link_libraries(GIDSpecBench PRIVATE gid_core)
  add_executable(GIDSimpleProbe tools/SimpleProbe.cpp)
  target_link_libraries(GIDSimpleProbe PRIVATE gid_core)
endif()

# ---------------------------------------------------------------- scrm self-check probe (off by default)
option(GID_BUILD_SCRM_PROBE "Build SCRM self-check probe (DIP + CSV)" OFF)
if(GID_BUILD_SCRM_PROBE)
  add_executable(GIDScrmProbe tools/ScrmProbe.cpp)
  target_link_libraries(GIDScrmProbe PRIVATE gid_core)
endif()

# CMake の最低バージョンを指定します。
cmake_minimum_required(VERSION 3.20)

# プロジェクト名と使用する言語を指定します。
project(mdspassproject CXX)

# C++ の標準を C++ 20 に設定します。
set(CMAKE_CXX_STANDARD 20)

# 必要なライブラリを検索します。
find_package(FLTK REQUIRED)
find_package(OpenGL REQUIRED)
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)
find_package(PNG REQUIRED)

# インクルード ディレクトリを指定します。
include_directories(${FLTK_INCLUDE_DIR})
include_directories(${OPENGL_INCLUDE_DIR})
include_directories(${LAPACK_INCLUDE_DIR})
include_directories(${BLAS_INCLUDE_DIR})
include_directories(${PNG_INCLUDE_DIRS})

# ビルド対象となるソース ファイルのリストを取得します。
file(GLOB SOURCES ./src/*.cpp)

# 実行ファイルの設定を作成します。
add_executable(mdspass3 ${SOURCES})

# 実行ファイルにリンクするライブラリを指定します。
target_link_libraries(mdspass3 ${FLTK_LIBRARIES})
target_link_libraries(mdspass3 ${OPENGL_LIBRARIES})
target_link_libraries(mdspass3 ${LAPACK_LIBRARIES})
target_link_libraries(mdspass3 ${BLAS_LIBRARIES})
target_link_libraries(mdspass3 ${PNG_LIBRARIES})

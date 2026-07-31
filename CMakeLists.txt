cmake_minimum_required(VERSION 3.18.3)

project(openBPM
  DESCRIPTION "A free and open source BPM solver"
  LANGUAGES C CXX)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)

# ── Windows (MSVC) 対応 ──────────────────────────────────────────────
# libm は存在しない (MATH_LIB を空にする)。ソースは UTF-8 (日本語コメント)。
# M_PI 等は _USE_MATH_DEFINES が必要。
if(MSVC)
    add_compile_options(/utf-8)
    add_compile_definitions(_USE_MATH_DEFINES _CRT_SECURE_NO_WARNINGS NOMINMAX)
    set(MATH_LIB "")
else()
    set(MATH_LIB m)
endif()

find_package(Threads REQUIRED)

find_package(Eigen3 REQUIRED)  # バージョン指定なし: Homebrew の Eigen 5.x も受理する

# CUDA を使用するかどうかを設定するオプションを追加
option(WITH_CUDA "Enable CUDA support" OFF)
if(WITH_CUDA)
    # CUDA は WITH_CUDA=ON のときだけ有効化する (nvcc の無い環境・CI でも
    # CPU 版をビルドできるようにする)
    enable_language(CUDA)
    #set(CMAKE_CUDA_HOST_COMPILER /usr/bin/g++)
    message(STATUS ${CMAKE_CUDA_COMPILER})
    message(STATUS ${CMAKE_CUDA_HOST_COMPILER})
    # CUDA を探す（CMake 3.18 より前）
    #find_package(CUDA REQUIRED)
    # CUDAToolkit を探す
    find_package(CUDAToolkit REQUIRED)
    set(CMAKE_CUDA_ARCHITECTURES 60)
    #message(STATUS "_MPI マクロが定義されました")
    # CUDA のコンパイラフラグを設定
    #set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -Xcompiler -pthread")
else()
    #message(STATUS "CUDA マクロは定義されていません")
endif()


# MPI を使用するかどうかを設定するオプションを追加
option(WITH_MPI "Enable MPI support" OFF)
# オプションに基づいて _MPI を定義
if(WITH_MPI)
    #before : -D_MPI
    find_package(MPI REQUIRED C)
    #+cuda
    #find_package(MPI REQUIRED)
    add_compile_definitions(_MPI)
    # MPI 版 (mpi/solve.c) は並列 HDF5 (H5Pset_fapl_mpio 等) を使うため、
    # 後段の find_package(HDF5) で並列ビルドの HDF5 を優先させる
    set(HDF5_PREFER_PARALLEL TRUE)
    message(STATUS "_MPI マクロが定義されました")
else()
    message(STATUS "_MPI マクロは定義されていません")
endif()

#対応 geometry数
#add_compile_definitions(_ID64)
add_compile_definitions(_ID32)
#add_compile_definitions(_ID16)

#set(CMAKE_SOURCE_DIR ..)
#include_directories(${CMAKE_SOURCE_DIR}/include)
#output folder
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})

#find_package(glfw3 REQUIRED)
#find_package(OpenGL REQUIRED)

# Find CUDA package
# use CUDA?
# Find OpenMP package
find_package(OpenMP REQUIRED)
# Find MPI package

if(WITH_MPI)
    file(GLOB MPI_SOURCES "mpi/*.c")
else()
    file(GLOB SOURCES "sol/*.c")
endif()

file(GLOB BPM_SOURCES "bpm/*.cpp")

# HDF5パッケージの検索

#conda
#set(HDF5_ROOT "$ENV{CONDA_PREFIX}")
find_package(HDF5 REQUIRED COMPONENTS C)
#find_package(HDF5 COMPONENTS)
#find_package(HDF5 COMPONENTS C)
#find_package(hdf5 REQUIRED CONFIG PATHS )
if (HDF5_FOUND)
    message(STATUS "HDF5 が見つかりました。HDF5 サポートを有効にします。")
    # インクルードディレクトリの追加
    #target_include_directories(TimeSeriesHDF5 PRIVATE ${HDF5_INCLUDE_DIRS})

    # ライブラリのリンク
    #target_link_libraries(TimeSeriesHDF5 PRIVATE ${HDF5_LIBRARIES})
    #https://forum.hdfgroup.org/t/cmake-parallel-hdf5-config-scripts-references-hdf-hdf-shared/4605/3
    #set(OP2_HDF5_DEFINITIONS ${HDF5_DEFINITIONS})
    #set(OP2_HDF5_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS})
    #set(OP2_HDF5_LIBRARIES ${HDF5_LIBRARIES})
    # If HDF5 is build with parallel support it requires MPI, which needs
    # to be enabled
    if (HDF5_IS_PARALLEL)
      if (WITH_MPI)
      #  set(OP2_HDF5_INCLUDE_DIRS ${OP2_HDF5_INCLUDE_DIRS} ${MPI_INCLUDE_PATH})
      #  set(OP2_HDF5_LIBRARIES ${OP2_HDF5_LIBRARIES} ${MPI_LIBRARIES})
      #  # MPI is always built with HDF5 if available
      #  set(OP2_MPI_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS} ${OP2_MPI_INCLUDE_DIRS})
      #  set(OP2_MPI_LIBRARIES ${HDF5_LIBRARIES} ${OP2_MPI_LIBRARIES})
        # 並列 HDF5 のヘッダ (H5public.h) は mpi.h を include するため、
        # MPI を使わないターゲット (obpm / obpm_post / テスト) が HDF5::HDF5 を
        # リンクしただけでは mpi.h が見つからずビルドが失敗する。
        # インポートターゲットの INTERFACE に MPI を足して全利用者へ伝播させる
        # (並列 libhdf5 自体が libmpi に依存するのでリンクも必要)。
        if (TARGET HDF5::HDF5 AND TARGET MPI::MPI_C)
          set_property(TARGET HDF5::HDF5 APPEND
                       PROPERTY INTERFACE_LINK_LIBRARIES MPI::MPI_C)
          message(STATUS "並列 HDF5 を検出: HDF5::HDF5 に MPI::MPI_C を伝播させます")
        endif()
      else()
        message(STATUS "HDF5 is built with parallel support requiring MPI, but MPI is disabled")
        message(STATUS "Disabling HDF5 support")
        #set(OP2_WITH_HDF5 OFF)
      endif()
    else()
      message(STATUS "HDF5 is not built with parallel support.")
    endif()
endif()

# CUDA が見つかった場合のみターゲットをビルド
if (WITH_CUDA)
    #if (CUDA_FOUND)
    if (CUDAToolkit_FOUND)
        #CC : nvcc
        message(STATUS "CUDA が見つかりました。CUDA サポートを有効にします。")
        #message(STATUS ${OBJECT_FILES})
        # CUDA ターゲットは以下の明示リスト (SOURCES2) で必要なソースに絞ってある。
        # CPU 版 (obpm) の file(GLOB SOURCES "sol/*.c") は全ソースをリンクする現状の
        # 構成で動作しており、これ以上の絞り込みはビルド時間の削減にしかならないため
        # 意図的に glob のまま維持している。

        # オブジェクトファイルを生成するカスタムコマンド
        #add_custom_command(
        #    OUTPUT ${OUTPUT_OBJECT}
        #    COMMAND ${CMAKE_C_COMPILER} -c ${SOURCE_FILE} -o ${OUTPUT_OBJECT} -Wall -O2
        #    DEPENDS ${SOURCE_FILE}
        #    COMMENT "Compiling ${SOURCE_FILE} to ${OUTPUT_OBJECT}"
        #    VERBATIM
        #)

        set(SOURCES2
            ${CMAKE_SOURCE_DIR}/sol/input_data.c
            ${CMAKE_SOURCE_DIR}/sol/input2.c
            ${CMAKE_SOURCE_DIR}/sol/monitor.c
            ${CMAKE_SOURCE_DIR}/sol/memallocfree1.c
            ${CMAKE_SOURCE_DIR}/sol/memallocfree3.c
            ${CMAKE_SOURCE_DIR}/sol/setup.c
            ${CMAKE_SOURCE_DIR}/sol/setupSize.c
            ${CMAKE_SOURCE_DIR}/sol/setupDispersion.c
            ${CMAKE_SOURCE_DIR}/sol/setupDft.c
            ${CMAKE_SOURCE_DIR}/sol/setupMurHx.c
            ${CMAKE_SOURCE_DIR}/sol/setupMurHy.c
            ${CMAKE_SOURCE_DIR}/sol/setupMurHz.c
            ${CMAKE_SOURCE_DIR}/sol/setupPmlEx.c
            ${CMAKE_SOURCE_DIR}/sol/setupPmlEy.c
            ${CMAKE_SOURCE_DIR}/sol/setupPmlEz.c
            ${CMAKE_SOURCE_DIR}/sol/setupPmlHx.c
            ${CMAKE_SOURCE_DIR}/sol/setupPmlHy.c
            ${CMAKE_SOURCE_DIR}/sol/setupPmlHz.c
            ${CMAKE_SOURCE_DIR}/sol/setupPml.c
            ${CMAKE_SOURCE_DIR}/sol/vfeed.c
            ${CMAKE_SOURCE_DIR}/sol/nearfield_c.c
            ${CMAKE_SOURCE_DIR}/sol/outputChars.c
            ${CMAKE_SOURCE_DIR}/sol/outputZin.c
            ${CMAKE_SOURCE_DIR}/sol/outputSpara.c
            ${CMAKE_SOURCE_DIR}/sol/outputCoupling.c
            ${CMAKE_SOURCE_DIR}/sol/outputCross.c
            ${CMAKE_SOURCE_DIR}/sol/setup_vector.c
            ${CMAKE_SOURCE_DIR}/sol/setupId.c
            ${CMAKE_SOURCE_DIR}/sol/ingeometry.c
            ${CMAKE_SOURCE_DIR}/sol/geomlines.c
            ${CMAKE_SOURCE_DIR}/sol/farfield.c
            ${CMAKE_SOURCE_DIR}/sol/writeout.c
            ${CMAKE_SOURCE_DIR}/sol/utils.c
            ${CMAKE_SOURCE_DIR}/sol/cputime.c
            ${CMAKE_SOURCE_DIR}/sol/plot3dGeom.c
            ${CMAKE_SOURCE_DIR}/sol/ev3d.c
        )

        # CUDA の設定を指定
        set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -gencode arch=compute_60,code=sm_60")

        # Add CUDA executable
        if (WITH_MPI)
            # Add CUDA source files(reject solve.cu)
            #file(GLOB CUDA_SOURCES "cuda/*.cu")
            set(CUDA_SOURCES2
                ${CMAKE_SOURCE_DIR}/cuda/average.cu
                ${CMAKE_SOURCE_DIR}/cuda/check_gpu.cu
                ${CMAKE_SOURCE_DIR}/cuda/cuda_memory.cu
                ${CMAKE_SOURCE_DIR}/cuda/dftNear3d.cu
                ${CMAKE_SOURCE_DIR}/cuda/dispersionEx.cu
                ${CMAKE_SOURCE_DIR}/cuda/dispersionEy.cu
                ${CMAKE_SOURCE_DIR}/cuda/dispersionEz.cu
                ${CMAKE_SOURCE_DIR}/cuda/efeed.cu
                ${CMAKE_SOURCE_DIR}/cuda/eload.cu
                ${CMAKE_SOURCE_DIR}/cuda/info_gpu.cu
                ${CMAKE_SOURCE_DIR}/cuda/initfield.cu
                ${CMAKE_SOURCE_DIR}/cuda/memallocfree2_gpu.cu
                ${CMAKE_SOURCE_DIR}/cuda/memallocfree3_gpu.cu
                ${CMAKE_SOURCE_DIR}/cuda/murH.cu
                ${CMAKE_SOURCE_DIR}/cuda/pbcx.cu
                ${CMAKE_SOURCE_DIR}/cuda/pbcy.cu
                ${CMAKE_SOURCE_DIR}/cuda/pbcz.cu
                ${CMAKE_SOURCE_DIR}/cuda/pmlEx.cu
                ${CMAKE_SOURCE_DIR}/cuda/pmlEy.cu
                ${CMAKE_SOURCE_DIR}/cuda/pmlEz.cu
                ${CMAKE_SOURCE_DIR}/cuda/pmlHx.cu
                ${CMAKE_SOURCE_DIR}/cuda/pmlHy.cu
                ${CMAKE_SOURCE_DIR}/cuda/pmlHz.cu
                ${CMAKE_SOURCE_DIR}/cuda/setup_gpu.cu
                ${CMAKE_SOURCE_DIR}/cuda/updateEx.cu
                ${CMAKE_SOURCE_DIR}/cuda/updateEy.cu
                ${CMAKE_SOURCE_DIR}/cuda/updateEz.cu
                ${CMAKE_SOURCE_DIR}/cuda/updateHx.cu
                ${CMAKE_SOURCE_DIR}/cuda/updateHy.cu
                ${CMAKE_SOURCE_DIR}/cuda/updateHz.cu
                ${CMAKE_SOURCE_DIR}/cuda/vpoint.cu
            )

            set(MPI_SOURCES2
                ${CMAKE_SOURCE_DIR}/mpi/comm.c
                ${CMAKE_SOURCE_DIR}/mpi/setup_mpi.c
            )

            file(GLOB CUDA_MPI_SOURCES "cuda_mpi/*.cu")
            add_library(obpm_cuda_libmpi ${SOURCES2} ${MPI_SOURCES2} ${CUDA_MPI_SOURCES} ${CUDA_SOURCES2})
            target_include_directories(obpm_cuda_libmpi PUBLIC "${PROJECT_SOURCE_DIR}/include")
            target_link_libraries(obpm_cuda_libmpi PRIVATE ${CUDA_LIBRARIES})
            target_link_libraries(obpm_cuda_libmpi PRIVATE OpenMP::OpenMP_C pthread m)
            target_link_libraries(obpm_cuda_libmpi PRIVATE MPI::MPI_C HDF5::HDF5)
            #error
            #target_link_libraries(obpm_cuda_libmpi PRIVATE MPI::MPI_C HDF5::HDF5-shared)

            #add_executable(obpm_cuda_mpi ${SOURCES2} ${MPI_SOURCES2} ${CUDA_MPI_SOURCES} ${CUDA_SOURCES2} src/cuda_mpi_Main.cu)
            add_executable(obpm_cuda_mpi src/cuda_mpi_Main.cu)

            #
            target_include_directories(obpm_cuda_mpi PUBLIC "${PROJECT_SOURCE_DIR}/include")

            # Link CUDA libraries
            #obpm_cuda_libmpi
            target_link_libraries(obpm_cuda_mpi PRIVATE obpm_cuda_libmpi)
            target_link_libraries(obpm_cuda_mpi PRIVATE ${CUDA_LIBRARIES})
            target_link_libraries(obpm_cuda_mpi PRIVATE OpenMP::OpenMP_C pthread m)
            target_link_libraries(obpm_cuda_mpi PRIVATE MPI::MPI_C HDF5::HDF5)
        else ()
            # Add CUDA source files
            file(GLOB CUDA_SOURCES "cuda/*.cu")
            set(CUDA_SOURCES
                ${CMAKE_SOURCE_DIR}/cuda/average.cu
                ${CMAKE_SOURCE_DIR}/cuda/check_gpu.cu
                ${CMAKE_SOURCE_DIR}/cuda/cuda_memory.cu
                ${CMAKE_SOURCE_DIR}/cuda/dftNear3d.cu
                ${CMAKE_SOURCE_DIR}/cuda/dispersionEx.cu
                ${CMAKE_SOURCE_DIR}/cuda/dispersionEy.cu
                ${CMAKE_SOURCE_DIR}/cuda/dispersionEz.cu
                ${CMAKE_SOURCE_DIR}/cuda/efeed.cu
                ${CMAKE_SOURCE_DIR}/cuda/eload.cu
                ${CMAKE_SOURCE_DIR}/cuda/info_gpu.cu
                ${CMAKE_SOURCE_DIR}/cuda/initfield.cu
                ${CMAKE_SOURCE_DIR}/cuda/memallocfree2_gpu.cu
                ${CMAKE_SOURCE_DIR}/cuda/memallocfree3_gpu.cu
                ${CMAKE_SOURCE_DIR}/cuda/murH.cu
                ${CMAKE_SOURCE_DIR}/cuda/pbcx.cu
                ${CMAKE_SOURCE_DIR}/cuda/pbcy.cu
                ${CMAKE_SOURCE_DIR}/cuda/pbcz.cu
                ${CMAKE_SOURCE_DIR}/cuda/pmlEx.cu
                ${CMAKE_SOURCE_DIR}/cuda/pmlEy.cu
                ${CMAKE_SOURCE_DIR}/cuda/pmlEz.cu
                ${CMAKE_SOURCE_DIR}/cuda/pmlHx.cu
                ${CMAKE_SOURCE_DIR}/cuda/pmlHy.cu
                ${CMAKE_SOURCE_DIR}/cuda/pmlHz.cu
                ${CMAKE_SOURCE_DIR}/cuda/setup_gpu.cu
                ${CMAKE_SOURCE_DIR}/cuda/solve.cu
                ${CMAKE_SOURCE_DIR}/cuda/solve_bpm.cu
                ${CMAKE_SOURCE_DIR}/cuda/updateEx.cu
                ${CMAKE_SOURCE_DIR}/cuda/updateEy.cu
                ${CMAKE_SOURCE_DIR}/cuda/updateEz.cu
                ${CMAKE_SOURCE_DIR}/cuda/updateHx.cu
                ${CMAKE_SOURCE_DIR}/cuda/updateHy.cu
                ${CMAKE_SOURCE_DIR}/cuda/updateHz.cu
                ${CMAKE_SOURCE_DIR}/cuda/vpoint.cu
            )
            #reduction_sum.cu

            #add_executable(obpm_cuda ${CUDA_SOURCES} src/cuda_Main.cu)
            #add_executable(obpm_cuda ${CUDA_SOURCES} ${OBJECT_FILES} src/cuda_Main.cu)
            add_executable(obpm_cuda ${SOURCES2} ${CUDA_SOURCES} bpm/FDBPMpropagator.cu bpm/wabpm.cu bpm/wabpm.cpp bpm/modes.cpp src/cuda_Main.cu)
            #add_executable(obpm_cuda ${SOURCES2} ${CUDA_SOURCES} ${OBJECT_FILES} src/cuda_Main.cu)

            target_include_directories(obpm_cuda PUBLIC "${PROJECT_SOURCE_DIR}/include")

            # Link CUDA libraries
            target_link_libraries(obpm_cuda PRIVATE ${CUDA_LIBRARIES})

            # オプションとして C++ 実行ファイルも追加（CUDA 無効時のデフォルト動作用）
            #set(DEFAULT_TARGET obpm_cuda)

            # Optional: Specify CUDA architectures
            #set(CUDA_ARCHITECTURES "60;70;75;80")
            #add_compile_options(--expt-relaxed-constexpr)

            # Link OpenMP libraries to the target
            if(OpenMP_C_FOUND)
                message(STATUS "OpenMP FOUND.")
                #target_link_libraries(obpm_cuda PRIVATE OpenMP::OpenMP_C)
                #target_link_libraries(obpm_cuda PRIVATE OpenMP::OpenMP_C m)
                #target_link_libraries(obpm_cuda PRIVATE OpenMP::OpenMP_C Threads::Threads)
                #target_link_libraries(obpm_cuda PRIVATE OpenMP::OpenMP_C Threads::Threads m)
                target_link_libraries(obpm_cuda PRIVATE OpenMP::OpenMP_C pthread)
                #target_link_libraries(obpm PRIVATE MPI::MPI_C OpenMP::OpenMP_C m)
                target_link_libraries(obpm_cuda PRIVATE m HDF5::HDF5)
            else()
                message(STATUS "OpenMP not FOUND.")
                #target_link_libraries(obpm_cuda PUBLIC m)
                #target_link_libraries(obpm_cuda PRIVATE Threads::Threads m)
                #target_link_libraries(obpm_cuda PRIVATE Threads::Threads)
                target_link_libraries(obpm_cuda PRIVATE pthread)
                target_link_libraries(obpm_cuda PRIVATE m HDF5::HDF5)
            endif()

        endif()

        # CUDA のプロパティ設定 (オプション)(x)
        #set_target_properties(obpm_cuda PROPERTIES
        #    CUDA_SEPARABLE_COMPILATION ON          # 分離コンパイルを有効にする(x?)
        #    #CUDA_ARCHITECTURES "60;70"            # 対応する CUDA アーキテクチャを指定
        #    CUDA_ARCHITECTURES "60"                # 対応する CUDA アーキテクチャを指定
        #    CUDA_STANDARD 14                       # CUDA の標準 (C++14) を指定
        #)
    else()
        message(WARNING "CUDA が見つかりませんでした。CUDA サポートは無効です。")
    endif()
else()
    if(WITH_MPI)
        #solve.c reject
        set(SOURCES3
            ${CMAKE_SOURCE_DIR}/sol/average.c
            ${CMAKE_SOURCE_DIR}/sol/cputime.c
            ${CMAKE_SOURCE_DIR}/sol/dftNear3d.c
            ${CMAKE_SOURCE_DIR}/sol/dispersionEx.c
            ${CMAKE_SOURCE_DIR}/sol/dispersionEy.c
            ${CMAKE_SOURCE_DIR}/sol/dispersionEz.c
            ${CMAKE_SOURCE_DIR}/sol/efeed.c
            ${CMAKE_SOURCE_DIR}/sol/eload.c
            ${CMAKE_SOURCE_DIR}/sol/ev3d.c
            ${CMAKE_SOURCE_DIR}/sol/farfield.c
            ${CMAKE_SOURCE_DIR}/sol/geomlines.c
            ${CMAKE_SOURCE_DIR}/sol/ingeometry.c
            ${CMAKE_SOURCE_DIR}/sol/initfield.c
            ${CMAKE_SOURCE_DIR}/sol/input_data.c
            ${CMAKE_SOURCE_DIR}/sol/input2.c
            ${CMAKE_SOURCE_DIR}/sol/memallocfree1.c
            ${CMAKE_SOURCE_DIR}/sol/memallocfree2.c
            ${CMAKE_SOURCE_DIR}/sol/memallocfree3.c
            ${CMAKE_SOURCE_DIR}/sol/monitor.c
            ${CMAKE_SOURCE_DIR}/sol/murH.c
            ${CMAKE_SOURCE_DIR}/sol/nearfield_c.c
            ${CMAKE_SOURCE_DIR}/sol/outputChars.c
            ${CMAKE_SOURCE_DIR}/sol/outputCoupling.c
            ${CMAKE_SOURCE_DIR}/sol/outputCross.c
            ${CMAKE_SOURCE_DIR}/sol/outputSpara.c
            ${CMAKE_SOURCE_DIR}/sol/outputZin.c
            ${CMAKE_SOURCE_DIR}/sol/pbcx.c
            ${CMAKE_SOURCE_DIR}/sol/pbcy.c
            ${CMAKE_SOURCE_DIR}/sol/pbcz.c
            ${CMAKE_SOURCE_DIR}/sol/plot3dGeom.c
            ${CMAKE_SOURCE_DIR}/sol/pmlEx.c
            ${CMAKE_SOURCE_DIR}/sol/pmlEy.c
            ${CMAKE_SOURCE_DIR}/sol/pmlEz.c
            ${CMAKE_SOURCE_DIR}/sol/pmlHx.c
            ${CMAKE_SOURCE_DIR}/sol/pmlHy.c
            ${CMAKE_SOURCE_DIR}/sol/pmlHz.c
            ${CMAKE_SOURCE_DIR}/sol/setup.c
            ${CMAKE_SOURCE_DIR}/sol/setup_vector.c
            ${CMAKE_SOURCE_DIR}/sol/setupDft.c
            ${CMAKE_SOURCE_DIR}/sol/setupDispersion.c
            ${CMAKE_SOURCE_DIR}/sol/setupId.c
            ${CMAKE_SOURCE_DIR}/sol/setupMurHx.c
            ${CMAKE_SOURCE_DIR}/sol/setupMurHy.c
            ${CMAKE_SOURCE_DIR}/sol/setupMurHz.c
            ${CMAKE_SOURCE_DIR}/sol/setupPml.c
            ${CMAKE_SOURCE_DIR}/sol/setupPmlEx.c
            ${CMAKE_SOURCE_DIR}/sol/setupPmlEy.c
            ${CMAKE_SOURCE_DIR}/sol/setupPmlEz.c
            ${CMAKE_SOURCE_DIR}/sol/setupPmlHx.c
            ${CMAKE_SOURCE_DIR}/sol/setupPmlHy.c
            ${CMAKE_SOURCE_DIR}/sol/setupPmlHz.c
            ${CMAKE_SOURCE_DIR}/sol/setupSize.c
            ${CMAKE_SOURCE_DIR}/sol/updateEx.c
            ${CMAKE_SOURCE_DIR}/sol/updateEy.c
            ${CMAKE_SOURCE_DIR}/sol/updateEz.c
            ${CMAKE_SOURCE_DIR}/sol/updateHx.c
            ${CMAKE_SOURCE_DIR}/sol/updateHy.c
            ${CMAKE_SOURCE_DIR}/sol/updateHz.c
            ${CMAKE_SOURCE_DIR}/sol/utils.c
            ${CMAKE_SOURCE_DIR}/sol/vfeed.c
            ${CMAKE_SOURCE_DIR}/sol/vpoint.c
            ${CMAKE_SOURCE_DIR}/sol/writeout.c
        )

        #file(GLOB MPI_SOURCES "mpi/*.c")
        #target_compile_definitions(obpm_mpi PRIVATE _MPI)
        add_executable(obpm_mpi ${SOURCES3} ${MPI_SOURCES} src/mpi_Main.c)

        # Include directories for header files
        target_include_directories(obpm_mpi PUBLIC "${PROJECT_SOURCE_DIR}/include")

        #target_link_libraries(obpm_mpi PRIVATE MPI::MPI_C OpenMP::OpenMP_C m)
        target_link_libraries(obpm_mpi PRIVATE OpenMP::OpenMP_C Threads::Threads m)
        target_link_libraries(obpm_mpi PRIVATE MPI::MPI_C HDF5::HDF5)
    else ()
        #add_executable(obpm ${SOURCES} sol/solve_bpm.cpp src/sol_Main.c)
        # MSVC の C コンパイラは C99 _Complex を扱えないため、
        # floatcomplex (std::complex) を使う C ファイルは C++ としてコンパイルする
        if(MSVC)
            set_source_files_properties(bpm/FDBPMpropagator.c PROPERTIES LANGUAGE CXX)
        endif()
        add_executable(obpm ${SOURCES} bpm/initfield.cpp bpm/FDBPMpropagator.c bpm/wabpm.cpp bpm/modes.cpp sol/solve_bpm.cpp src/sol_Main.cpp)
        #target_link_libraries(obpm glfw)
        #target_link_libraries(obpm OpenGL::GL)

        # Include directories for header files
        target_include_directories(obpm PUBLIC "${PROJECT_SOURCE_DIR}/include")
        #target_link_libraries(obpm ${CMAKE_DL_LIBS})
        #target_link_libraries(obpm ${CMAKE_DL_LIBS})
        # Link OpenMP libraries to the target
        if(OpenMP_C_FOUND)
            #target_link_libraries(obpm PRIVATE OpenMP::OpenMP_C m)
            target_link_libraries(obpm PRIVATE OpenMP::OpenMP_C ${MATH_LIB} HDF5::HDF5 Eigen3::Eigen)
            #target_link_libraries(obpm PRIVATE OpenMP::OpenMP_CXX m)
        else()
            #target_link_libraries(obpm PRIVATE m)
            target_link_libraries(obpm PRIVATE ${MATH_LIB} HDF5::HDF5)
        endif()

        # ターゲットのオブジェクトディレクトリを設定
        set_target_properties(obpm PROPERTIES
            OBJECT_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/obj"
        )
    endif()
endif()

option(WITH_POST "Enable Prepost Tool" ON)
if (WITH_POST)
    file(GLOB POST_SOURCES "post/*.c")
    set(SOURCES3
        ${CMAKE_SOURCE_DIR}/sol/nearfield_c.c
        ${CMAKE_SOURCE_DIR}/sol/farfield.c
        ${CMAKE_SOURCE_DIR}/sol/outputCoupling.c
        ${CMAKE_SOURCE_DIR}/sol/utils.c
        ${CMAKE_SOURCE_DIR}/sol/ev3d.c
    )

    add_executable(obpm_post ${POST_SOURCES} ${SOURCES3} src/post_Main.c)

    # Include directories for header files
    target_include_directories(obpm_post PUBLIC "${PROJECT_SOURCE_DIR}/include")

    #target_link_libraries(obpm_post PRIVATE OpenMP::OpenMP_C m)
    target_link_libraries(obpm_post PRIVATE OpenMP::OpenMP_C ${MATH_LIB} HDF5::HDF5)

endif()

# ---------------------------------------------------------------
# テスト (解析解との比較検証)
#   cmake -S . -B build -DWITH_TESTS=ON
#   cmake --build build && ctest --test-dir build --output-on-failure
#   - test_wabpm  : BPM 伝搬カーネル (回折/エネルギー保存/吸収/曲げ/半ベクトル)
#   - test_modes  : モードソルバ (LP01/LP11 の実効屈折率・直交性)
#   - test_allset : P_Struct API (findModes/modeSuperposition/offset/tilt)
# ---------------------------------------------------------------
option(WITH_TESTS "Build unit tests" OFF)
if (WITH_TESTS)
    enable_testing()

    add_executable(test_wabpm tests/test_wabpm.cpp bpm/wabpm.cpp)
    add_executable(test_modes tests/test_modes.cpp bpm/modes.cpp bpm/wabpm.cpp)
    add_executable(test_allset tests/test_allset.cpp bpm/modes.cpp bpm/wabpm.cpp)
    # スカラー近軸カーネルは C (float _Complex) のため、テストも C で書く。
    # MSVC の C コンパイラは C99 _Complex を扱えないため、カーネル本体と同様に
    # C++ (std::complex<float>) としてコンパイルする (bpm_prototype.h が切替)
    if(MSVC)
        set_source_files_properties(tests/test_fdbpm.c PROPERTIES LANGUAGE CXX)
    endif()
    add_executable(test_fdbpm tests/test_fdbpm.c bpm/FDBPMpropagator.c)
    target_link_libraries(test_fdbpm PRIVATE ${MATH_LIB})
    target_link_libraries(test_allset PRIVATE Eigen3::Eigen)
    foreach(t test_wabpm test_modes test_allset test_fdbpm)
        target_include_directories(${t} PRIVATE "${PROJECT_SOURCE_DIR}/include")
        if (OpenMP_CXX_FOUND)
            target_link_libraries(${t} PRIVATE OpenMP::OpenMP_CXX)
        endif()
    endforeach()

    add_test(NAME wabpm_unit COMMAND test_wabpm)
    add_test(NAME modes_unit COMMAND test_modes)
    add_test(NAME allset_unit COMMAND test_allset)
    add_test(NAME fdbpm_unit COMMAND test_fdbpm)

    # ONN 光活性化関数 (TPA + powersweep) の結合テスト:
    # ソルバーを実際に走らせ、activation_curve.csv を平面波近似の解析解
    # T = 1/(1 + beta*(P_in/A_eff)*L) と比較する (Python は標準ライブラリのみ)。
    # obpm 本体を使うため WITH_MPI/WITH_CUDA 無効の CPU ビルドでのみ登録する。
    find_package(Python3 COMPONENTS Interpreter QUIET)
    if (Python3_FOUND AND TARGET obpm)
        set(ONN_TEST_DIR "${CMAKE_BINARY_DIR}/onn_test")
        file(MAKE_DIRECTORY "${ONN_TEST_DIR}")
        add_test(NAME onn_activation_solver
                 COMMAND $<TARGET_FILE:obpm> -n 2 -no-fdtd-out
                         "${CMAKE_SOURCE_DIR}/data/sample/onn_activation.ofd"
                 WORKING_DIRECTORY "${ONN_TEST_DIR}")
        add_test(NAME onn_activation_check
                 COMMAND ${Python3_EXECUTABLE}
                         "${CMAKE_SOURCE_DIR}/tests/check_activation.py"
                         --dir "${ONN_TEST_DIR}"
                 WORKING_DIRECTORY "${ONN_TEST_DIR}")
        # ソルバー実行 -> 検証 の順序を保証する
        set_tests_properties(onn_activation_solver PROPERTIES FIXTURES_SETUP onn_run)
        set_tests_properties(onn_activation_check  PROPERTIES FIXTURES_REQUIRED onn_run)

        # 多モード重ね合わせ励振 (launch = mode 0 1) : モードビート長を
        # 理論値 lambda/|neff0-neff1| と比較する
        set(BEAT_TEST_DIR "${CMAKE_BINARY_DIR}/beat_test")
        file(MAKE_DIRECTORY "${BEAT_TEST_DIR}")
        add_test(NAME mode_beat_solver
                 COMMAND $<TARGET_FILE:obpm> -n 2 -no-fdtd-out
                         "${CMAKE_SOURCE_DIR}/data/fiber_mode_superposition.ofd"
                 WORKING_DIRECTORY "${BEAT_TEST_DIR}")
        add_test(NAME mode_beat_check
                 COMMAND ${Python3_EXECUTABLE}
                         "${CMAKE_SOURCE_DIR}/tests/check_modes_beat.py"
                         --h5 "${BEAT_TEST_DIR}/time_series_data.h5"
                         --log "${BEAT_TEST_DIR}/obpm.log"
                 WORKING_DIRECTORY "${BEAT_TEST_DIR}")
        set_tests_properties(mode_beat_solver PROPERTIES FIXTURES_SETUP beat_run)
        set_tests_properties(mode_beat_check  PROPERTIES FIXTURES_REQUIRED beat_run)

        # 波長掃引 (wlsweep) : 掃引結果が単独実行と厳密に一致することを確認する
        # (スクリプト内でソルバーを 2 回実行する)
        add_test(NAME wlsweep_check
                 COMMAND ${Python3_EXECUTABLE}
                         "${CMAKE_SOURCE_DIR}/tests/check_wlsweep.py"
                         --bin $<TARGET_FILE:obpm>
                         --ofd "${CMAKE_SOURCE_DIR}/data/spectrum_sweep.ofd"
                         --workdir "${CMAKE_BINARY_DIR}/wlsweep_test")
    endif()
endif()

#custom viewer
#option(WITH_HMI "Enable HMI Tool" ON)
#if (WITH_HMI)
#    file(GLOB HMI_SOURCES "hmi/*.c")
#    set(SOURCES4
#        ${CMAKE_SOURCE_DIR}/sol/nearfield_c.c
#        ${CMAKE_SOURCE_DIR}/sol/farfield.c
#        ${CMAKE_SOURCE_DIR}/sol/outputCoupling.c
#        ${CMAKE_SOURCE_DIR}/sol/utils.c
#        ${CMAKE_SOURCE_DIR}/sol/ev3d.c
#    )
#
#    add_executable(obpm_post ${HMI_SOURCES} ${SOURCES4} src/post_Main.c)
#
#    # Include directories for header files
#    target_include_directories(obpm_post PUBLIC "${PROJECT_SOURCE_DIR}/include")
#
#    #target_link_libraries(obpm_post PRIVATE OpenMP::OpenMP_C m)
#    target_link_libraries(obpm_post PRIVATE OpenMP::OpenMP_C m HDF5::HDF5)
#
#endif()

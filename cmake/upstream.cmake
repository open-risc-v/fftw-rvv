# 上游接入：fftw（autotools）→ <build>/pkg/prefix。

set(PKG_ROOT "${CMAKE_BINARY_DIR}/pkg")
if(NOT DEFINED PKG_PREFIX)
  set(PKG_PREFIX "${PKG_ROOT}/prefix")
  add_custom_target(pkg ALL
    COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/cmake/pkg-build.sh" "${SRC_DIR}" "${PKG_ROOT}"
            "$<IF:$<BOOL:${RVV_KERNELS}>,${RVV_MARCH},rv64gc>"
            "$<IF:$<BOOL:${RVV_KERNELS}>,1,0>"
    COMMENT "构建上游 fftw（${SRC_DIR}）"
    VERBATIM)
  set(HAVE_PKG_TARGET TRUE)
endif()

add_executable(fftw_harness test/pkg_fftw.c)
if(HAVE_PKG_TARGET)
  add_dependencies(fftw_harness pkg)
endif()
target_include_directories(fftw_harness PRIVATE "${PKG_PREFIX}/include")
target_link_directories(fftw_harness PRIVATE "${PKG_PREFIX}/lib")
target_link_libraries(fftw_harness PRIVATE fftw3 m)
target_compile_options(fftw_harness PRIVATE ${RVV_FLAGS})

add_test(NAME fftw_harness COMMAND fftw_harness)

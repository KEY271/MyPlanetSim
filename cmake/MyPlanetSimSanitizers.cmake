option(MPS_ENABLE_SANITIZERS "Enable address and undefined behavior sanitizers" OFF)

add_library(myplanetsim_sanitizers INTERFACE)

if(MPS_ENABLE_SANITIZERS)
  if(MSVC)
    message(FATAL_ERROR "MPS_ENABLE_SANITIZERS is not supported with MSVC")
  endif()

  if(NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    message(
      FATAL_ERROR
      "MPS_ENABLE_SANITIZERS requires a Clang or GCC compatible compiler"
    )
  endif()

  target_compile_options(
    myplanetsim_sanitizers
    INTERFACE -fsanitize=address,undefined -fno-omit-frame-pointer
  )
  target_link_options(
    myplanetsim_sanitizers
    INTERFACE -fsanitize=address,undefined -fno-omit-frame-pointer
  )
endif()

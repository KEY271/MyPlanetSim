option(MPS_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)

add_library(myplanetsim_warnings INTERFACE)

if(MSVC)
  target_compile_options(
    myplanetsim_warnings
    INTERFACE /W4 /permissive-
  )
  if(MPS_WARNINGS_AS_ERRORS)
    target_compile_options(myplanetsim_warnings INTERFACE /WX)
  endif()
else()
  target_compile_options(
    myplanetsim_warnings
    INTERFACE -Wall -Wextra -Wpedantic
  )
  if(MPS_WARNINGS_AS_ERRORS)
    target_compile_options(myplanetsim_warnings INTERFACE -Werror)
  endif()
endif()

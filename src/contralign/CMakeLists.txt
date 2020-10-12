add_definitions(-DRNA=1)
add_library(contralign STATIC
    wrapper.cpp MultiSequence.cpp Sequence.cpp Utilities.cpp
)
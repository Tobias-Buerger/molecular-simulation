# Defining the executable target
add_executable(01 main.cpp)

# Linking against our MD code (propagates all dependencies)
target_link_libraries(01 PUBLIC my_md_lib)

# Copying "input" file simulation_test_input.xyz to the same directory as executable
add_input_file(01 simulation_test_input.txt)

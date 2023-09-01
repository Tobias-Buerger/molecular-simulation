# Defining the executable target
add_executable(07 main.cpp)

# Linking against our MD code (propagates all dependencies)
target_link_libraries(07 PUBLIC my_md_lib argparse)

add_input_file(07 cluster_923.xyz)
add_input_file(07 cluster_1415.xyz)
add_input_file(07 cluster_2057.xyz)
add_input_file(07 cluster_2869.xyz)
add_input_file(07 cluster_3871.xyz)

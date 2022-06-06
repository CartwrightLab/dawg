libryml = static_library(
	'ryml', ['rapidyaml/rapidyaml.hpp'],
	include_directories : include_directories('rapidyaml')
)

ryml_dep = declare_dependency(
	include_directories : include_directories('rapidyaml'),
	link_with : libryml
)

librandom = static_library(
    'random', ['random/random.cpp'],
    include_directories : include_directories('random')
)

random_dep = declare_dependency(
    include_directories : include_directories('random'),
    link_with : librandom
)

cli_dep = declare_dependency(include_directories : include_directories('CLI'))

doctest_dep = declare_dependency(include_directories: include_directories('doctest'))

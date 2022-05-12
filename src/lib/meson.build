
libdawg_sources = files([
  'ma.cpp',
  'matic.cpp',
  'models.cpp',
  'output.cpp',
  'parse.cpp'
])

libdawg_deps = [boost_dep, cli_dep, doctest_dep, random_dep, ryml_dep]

libdawg = static_library('libdawg', [libdawg_sources, version_file],
  include_directories : inc,
  dependencies : libdawg_deps,
  cpp_args : ['-DDOCTEST_CONFIG_DISABLE']
)

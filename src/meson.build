# Process version information
semver = meson.project_version().split('.')
version_src = configure_file(
  input : 'version.h.in.in',
  output : 'version.h.in',
  configuration : {
    'version_major' : semver[0],
    'version_minor' : semver[1],
    'version_patch' : semver[2],
    'version_int' : (semver[0].to_int()*1000+ semver[1].to_int())*1000+ \
                     semver[2].to_int(),
    'version_vcs' : '@VCS_TAG@',
    'package_string' : 'dawg ' + '@VCS_TAG@'
  })

version_file = vcs_tag(input : version_src, output : 'version.h')

inc = include_directories(['.', 'include'])

boost_dep = dependency('boost',
  version : '>=1.47.0',
  required : true
)

subdir('lib')

executable('dawg', ['dawg.cpp', version_file],
  include_directories : inc,
  link_with : [libdawg],
  cpp_args : ['-DDOCTEST_CONFIG_DISABLE'],
  dependencies : [libdawg_deps],
  install : true
)


# Dawg's Project Configuration

# Init Project
project('Dawg', 'cpp',
  version : '2.0.9000',
  license : 'GPL-2.0-or-later',
  meson_version : '>=0.58.0',
  default_options : [ 'buildtype=debugoptimized', 'cpp_std=c++17' ]
)

subdir('contrib')
subdir('src')
subdir('doc')
# subdir('tests')

configure_file(input : 'dawg.sublime-project.in',
  output : 'dawg.sublime-project',
  configuration : {
    'source_dir' : meson.project_source_root(),
    'project_name' : meson.project_name(),
    'build_dir' : meson.project_build_root()
  }
)

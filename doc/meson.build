fs = import('fs')

doxygen = find_program('doxygen', required : false)

if doxygen.found()
	semver = meson.project_version().split('.')

	doxy_conf = configuration_data()
	doxy_conf.set('DOXYGEN_INPUT_DIR', meson.project_source_root())
	doxy_conf.set('DOXYGEN_OUTPUT_DIR', meson.current_build_dir())
	doxy_conf.set('DAWG_DESCRIPTION', 'Sequence Evolution Simulator')
	doxy_conf.set('SOURCE_DIR', meson.project_source_root())
	doxy_conf.set('DAWG_VERSION_MAJOR', semver[0])
	doxy_conf.set('DAWG_VERSION_MINOR', semver[1])
	doxy_conf.set('DAWG_VERSION_PATCH', semver[2])

	doxyfile = configure_file(
		input : 'doxyfile.in',
		output : 'doxyfile',
		configuration : doxy_conf,
		install : false
	)

	html_header = configure_file(
		input : 'header.html.in',
		output : 'header.html',
		configuration : doxy_conf,
		install : false
	)

	custom_target('doc',
		command : [doxygen, doxyfile, html_header],
		output : ['html'],
		build_by_default : false
	)
else
	message('Documentation disabled, doxygen not found')
endif

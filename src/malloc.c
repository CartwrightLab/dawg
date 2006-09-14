#ifdef HAVE_CONFIG_H
#	include <config.h>
#endif

#undef malloc
#undef realloc

#ifdef HAVE_SYS_TYPES_H
#	include <sys/types.h>
#endif
#ifdef HAVE_STDDEF_H
#	include <stddef.h>
#endif
#ifdef HAVE_STDLIB_H
#	include <stdlib.h>
#endif

void *rpl_malloc(size_t n)
{
	if (n == 0) n = 1;
	return malloc(n);
} 

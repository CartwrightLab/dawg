#include <stdio.h>
#include <stdarg.h>

bool DawgError(const char* csErr, ...)
{
	fprintf(stderr, "Error: ");
	va_list args;
	va_start(args, csErr);
	vfprintf(stderr, csErr, args);
	va_end(args);
	fprintf(stderr, "\n");
	return false;
}

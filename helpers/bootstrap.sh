#! /bin/sh

aclocal \
&& autoheader \
&& automake --gnu --add-missing --warning=none \
&& autoconf -f 

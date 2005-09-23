#! /bin/sh

aclocal 2> /dev/null \
&& autoheader --warnings=none \
&& automake --gnu --add-missing --warnings=none \
&& autoconf -f --warnings=none

#! /bin/sh

: ${DAWG=dawg}
: ${DIFF=diff}

${DAWG} ${top_srcdir}/tests/test0.dawg | ${DIFF} - \
	${top_srcdir}/tests/test0.fasta
result=$?

exit $result

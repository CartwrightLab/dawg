#ifndef BISON_PARSER_HPP
# define BISON_PARSER_HPP

#ifndef YYSTYPE
typedef union {
	double d;	/* number values */
	char  cs[1024];  /* string values */
	Node  *node;
	char ch;
} yystype;
# define YYSTYPE yystype
# define YYSTYPE_IS_TRIVIAL 1
#endif
# define	LENGTH	257
# define	LABEL	258
# define	END	259
# define	UNKNOWN	260
# define	LPARTH	261
# define	RPARTH	262
# define	COLON	263
# define	SEMICOLON	264


extern YYSTYPE yylval;

#endif /* not BISON_PARSER_HPP */

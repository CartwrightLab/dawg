#ifndef BISON_PARSER_HPP
# define BISON_PARSER_HPP

#ifndef YYSTYPE
typedef union {
	double d;	/* number values */
	char*  cs;  /* string values */
	char   ch;  /* characters */
	bool   b;   /* booleans */
	DawgVar::Vec *pvec; /*vector*/
	DawgVar *pvar; /*DawgVar*/
	Node	*pnode; /*Tree*/
} yystype;
# define YYSTYPE yystype
# define YYSTYPE_IS_TRIVIAL 1
#endif
# define	NUM	257
# define	STRING	258
# define	DOT	259
# define	LABEL	260
# define	EQ	261
# define	END	262
# define	ID	263
# define	BOOL	264
# define	COLON	265
# define	LBRACE	266
# define	RBRACE	267
# define	LPARTH	268
# define	RPARTH	269
# define	UNKNOWN	270


extern YYSTYPE yylval;

#endif /* not BISON_PARSER_HPP */

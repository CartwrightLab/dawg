#ifndef BISON_PARSER_HPP
# define BISON_PARSER_HPP

#ifndef YYSTYPE
typedef union {
	double d;	/* number values */
	char  cs[1024];  /* string values */
	char   ch;  /* characters */
	bool   b;   /* booleans */
	DawgVar::Vec *pvec; /*vector*/
	DawgVar *pvar; /*DawgVar*/
	Node	*pnode; /*Tree*/
	std::string* pstr;
} yystype;
# define YYSTYPE yystype
# define YYSTYPE_IS_TRIVIAL 1
#endif
# define	NUM	257
# define	LENGTH	258
# define	STRING	259
# define	LABEL	260
# define	ID	261
# define	BOOL	262
# define	DOT	263
# define	EQ	264
# define	LBRACE	265
# define	RBRACE	266
# define	LPARTH	267
# define	RPARTH	268
# define	UNKNOWN	269
# define	END	270


extern YYSTYPE yylval;

#endif /* not BISON_PARSER_HPP */

/* A Bison parser, made by GNU Bison 1.875.  */

/* Skeleton parser for Yacc-like parsing with Bison,
   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002 Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330,
   Boston, MA 02111-1307, USA.  */

/* As a special exception, when this file is copied by Bison into a
   Bison output file, you may use that output file without restriction.
   This special exception was added by the Free Software Foundation
   in version 1.24 of Bison.  */

/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     NUM = 258,
     STRING = 259,
     DOT = 260,
     LABEL = 261,
     EQ = 262,
     END = 263,
     ID = 264,
     BOOL = 265,
     COLON = 266,
     LBRACE = 267,
     RBRACE = 268,
     LPARTH = 269,
     RPARTH = 270,
     UNKNOWN = 271
   };
#endif
#define NUM 258
#define STRING 259
#define DOT 260
#define LABEL 261
#define EQ 262
#define END 263
#define ID 264
#define BOOL 265
#define COLON 266
#define LBRACE 267
#define RBRACE 268
#define LPARTH 269
#define RPARTH 270
#define UNKNOWN 271




#if ! defined (YYSTYPE) && ! defined (YYSTYPE_IS_DECLARED)
#line 20 "parser.b"
typedef union YYSTYPE {
	double d;	/* number values */
	char*  cs;  /* string values */
	char   ch;  /* characters */
	bool   b;   /* booleans */
	DawgVar::Vec *pvec; /*vector*/
	DawgVar *pvar; /*DawgVar*/
	Node	*pnode; /*Tree*/
} YYSTYPE;
/* Line 1248 of yacc.c.  */
#line 78 "parser.hpp"
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE yylval;




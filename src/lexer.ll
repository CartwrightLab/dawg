%{
// lexer.ll - Copyright (C) 2004-2005 Reed A. Cartwright (all rights reserved)

#include "dawg.h"
#include "var.h"
#include "parser.h"

using namespace std;

bool g_bTerminate = false;
bool g_bParseOkay = true;
size_t g_uLine = 1; 
VarDB *g_pDB = NULL;

#ifdef _MSC_VER
#	pragma warning(disable: 4127 4244 4267 )
#	if _MSC_VER >= 1400
#		define isatty _isatty
#		define fileno _fileno	
#	endif
#endif

int yyparse();

bool RunParser(FILE *fin, VarDB *db)
{
	g_pDB = db;
	yyin = fin;
	g_bParseOkay = true;
	g_bTerminate = false;
	g_uLine = 1;
	yyparse();
	g_pDB = NULL;
	yyin = NULL;
	
	return g_bParseOkay;
}

void yyerror (char *s)
{
	g_bParseOkay = false;
	g_pDB->ParseError(s, g_uLine, yytext);

}

%}

%option nounput
%option noyywrap

DIGIT  [0-9]
SPACE [ \t\r\v\f]
BIDWORD ^{SPACE}*[A-Za-z][A-Za-z_0-9.]*
IDWORD [A-Za-z.][A-Za-z_0-9.]*
LABELCH [^ \t\n\r\v\f\(\)\[\]:;,\'\"]
NUMBER [-+]?{DIGIT}+("."{DIGIT}+)?([eE][+-]?{DIGIT}+)?
STR  \"[^\"\n]*\"|\'[^\'\n]*\'

%x tree
%x quote

%%

[=\{\}] {
	yylval.ch = yytext[0];
	return yytext[0];
}

^\s*"[" {
	yylval.ch = '[';
	return yytext[0];
}

"]"\s*\r?$ {
	yylval.ch = '[';
	return yytext[0];
}

[?+]"=" {
	yylval.ch = yytext[0];
	return yytext[0];
}

[Ff][Aa][Ll][Ss][Ee]|[Nn][Oo]|[Nn]|[Ff] {
	yylval.b = false;
	return BOOL;
}

[Tt][Rr][Uu][Ee]|[Yy][Ee][Ss]|[Tt]|[Yy] {
	yylval.b = true;
	return BOOL;
}

{BIDWORD} {
	yylval.cs = strdup(yytext+strspn(yytext, " \t\r\v\f"));
	return BID;
}

{IDWORD} {
	yylval.cs = strdup(yytext);
	return ID;
}


{NUMBER} {
	yylval.d = atof(yytext);
	return NUM;
}

\"\n(.*\n)*\"\n {
	yytext[strlen(yytext)-3] = '\0';
	yylval.cs = strdup(yytext+2);
	return STRING;
}


{STR} {
	yytext[strlen(yytext)-1] = '\0';
	yylval.cs = strdup(yytext+1);
	return STRING;
}

<INITIAL>"\"\"\""\n? {
	yylval.ch = yytext[0];
	g_uLine++;
	BEGIN(quote);
	return DQUOTE;
}

<quote>\n?"\"\"\"" {
	yylval.ch = yytext[0];
	g_uLine++;
	BEGIN(INITIAL);
	return DQUOTE;
}

<INITIAL>"\'\'\'" {
	yylval.ch = yytext[0];
	BEGIN(quote);
	return SQUOTE;
}

<quote>"\'\'\'" {
	yylval.ch = yytext[0];
	BEGIN(INITIAL);
	return SQUOTE;
}

<quote>\n {
	yylval.ch = yytext[0];
	g_uLine++;
	return CHAR;
}

<quote>\r {
	//ignore
}

<quote>. {
	yylval.ch = yytext[0];
	return CHAR;
}

"#".* | 
"//".* {
	// Comments
}

"(" {
	yylval.ch = yytext[0];
	BEGIN(tree);
	return yytext[0];
}

<tree>[\(\)] {
	yylval.ch = yytext[0];
	return yytext[0];
}

<tree>";" { BEGIN(INITIAL); }

<tree>":"{NUMBER} {
	yylval.d = atof(yytext+1);
	return LENGTH;
}

<tree>{STR} {
	yytext[strlen(yytext)-1] = '\0';
	yylval.cs = strdup(yytext+1);
	return LABEL;
}

<tree>{LABELCH}+ {
	yylval.cs = strdup(yytext);
	return LABEL;	
}
<tree>"[".+"]" { }

<*><<EOF>> {
	//if(g_bTerminate)
		yyterminate();
	//g_bTerminate = true;
	//return END;
}

";" { }

<*>[, \t\r\v\f]+ { }

<*>\n {
	g_uLine++;
}

<*>. {
	//yylval.ch = yytext[0];
	return UNKNOWN;	
}

%%

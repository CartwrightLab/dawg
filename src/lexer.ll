%{
// lexer.ll - Copyright (C) 2004-2005 Reed A. Cartwright (all rights reserved)

#include "dawg.h"
#include "var.h"
#include "parser.h"

using namespace std;
struct State
{
	int    nLine;	
	string ssFile;
} g_state;

bool g_bParseOkay = true;
bool g_bTerminate = false;

void yyerror (char *s)
{
	g_bParseOkay = false;
	cerr << "ALERT: " << s << " in " << g_state.ssFile << " at line " << g_state.nLine;
	cerr << ": \"" << yytext << "\"." << endl;
}
int yyparse (void);
bool Parse(const char* cs)
{
	FILE* stream = (cs==NULL || !strcmp(cs, "-")) ? stdin : fopen(cs, "r");
	if(stream == NULL)
		return false;
	g_state.nLine = 1;
	g_state.ssFile = (cs==NULL || !strcmp(cs, "-")) ? "stdin" : cs;
	yyin = stream;
	g_bParseOkay = true;
	yyparse();
	if(cs!=NULL)
		fclose(stream);
	return g_bParseOkay;
}

%}

%option nounput
%option noyywrap

DIGIT  [0-9]
IDWORD [A-Za-z][A-Za-z_0-9.]*
LABELCH [^ \t\n\r\v\f\(\)\[\]:;,\'\"]
NUMBER [-+]?{DIGIT}+("."{DIGIT}+)?([eE][+-]?{DIGIT}+)?
SPACE [ \t\r\v\f]
STR  \"[^\"\n]*\"|\'[^\'\n]*\'

%x tree
%x tostr

%%

[=\{\}\[\]] {
	yylval.ch = yytext[0];
	return yytext[0];
}

[?+]"=" {
	yylval.ch = yytext[0];
	return yytext[0];
}

[Ff][Aa][Ll][Ss][Ee] {
	yylval.b = false;
	return BOOL;
}

[Tt][Rr][Uu][Ee] {
	yylval.b = true;
	return BOOL;
}

{IDWORD} {
	yylval.pss = new string(yytext);
	return ID;
}

{NUMBER} {
	yylval.d = atof(yytext);
	return NUM;
}

\"\n(.*\n)*\"\n {
	yytext[strlen(yytext)-3] = '\0';
	yylval.pss = new string(yytext+2);
	return STRING;
}


{STR} {
	yytext[strlen(yytext)-1] = '\0';
	yylval.pss = new string(yytext+1);
	return STRING;
}

"<<"{IDWORD}{SPACE}*\n {
	yytext += 2;

	int s;
	for(s=0;!isspace(yytext[s]);++s) { }
	yytext[s] = '\0';
	
	yylval.pss = new string;
	string ssTemp;
	string ssEnd(yytext);
	while(1)
	{
		int c = yyinput();
		if(c == '\r')
			continue;
		if(c == '\n' || c == EOF)
		{
			g_state.nLine++;
			if(ssTemp == ssEnd)
				break;
			yylval.pss->append(ssTemp);
			yylval.pss->append("\n");
			ssTemp.clear();
		}
		else
		{
			ssTemp += c;
		}
		
		if(c == EOF)
			return UNKNOWN;
	}
	return STRING;
}

"#"[^\n]* | 
"//"[^\n]* {
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
	yylval.pss = new string(yytext+1);
	return LABEL;
}

<tree>{LABELCH}+ {
	yylval.pss = new string(yytext);
	return LABEL;	
}
<tree>"[".+"]" { }

<*><<EOF>> {
	if(g_bTerminate)
		yyterminate();
	g_bTerminate = true;
	return END;
}

\n {
	g_state.nLine++;
	yylval.ch = yytext[0];
	return yytext[0];
}

";" { }

<*>[, \t\r\v\f]+ { }

<*>\n {
	g_state.nLine++;
}

<*>. {
	//yylval.ch = yytext[0];
	return UNKNOWN;	
}

%%

%{
// cmd: "flex -t -CF lexer.fl > lexer.cpp"
//
// Error Warning: Does not check for circular includes.

#include "dawg.h"
#include "var.h"
#include "parser.h"

#pragma warning(disable: 4127 4244)

using namespace std;
struct State
{
	int    nLine;	
	string ssFile;
} g_state;

void yyerror (char *s)
{
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
	yyparse();
	if(cs!=NULL)
		fclose(stream);
	return true;
}

%}

%option nounput
%option noyywrap

DIGIT  [0-9]
IDWORD [A-Za-z][A-Za-z_0-9]*
STR    \"[^\"\n]*\"
LABELCH [^ \t\n\r\v\f\(\)\[\]:;,\'\"]
NUMBER [-+]?{DIGIT}+("."{DIGIT}+)?([eE][+-]?{DIGIT}+)?

%x tree

%%

[Ff]"alse" {
	yylval.b = false;
	return BOOL;
}

[Tt]"rue" {
	yylval.b = true;
	return BOOL;
}

{IDWORD} {
	strncpy(yylval.cs, yytext, 1023);
	yylval.cs[1023] = '\0';
	return ID;
}

{NUMBER} {
	yylval.d = atof(yytext);
	return NUM;
}

{STR} {
	size_t t = strlen(yytext);
	yytext[t-1] = '\0';
	strncpy(yylval.cs, yytext+1, 1023);
	yylval.cs[1023] = '\0';
	return STRING;
}

[=.\{\}] {
	yylval.ch = yytext[0];
	return yytext[0];
}

"#"[^\n]* {
}

"/""/"[^\n]* {
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
	size_t t = strlen(yytext);
	yytext[t-1] = '\0';
	strncpy(yylval.cs, yytext+1, 1023);
	yylval.cs[1023] = '\0';
	return LABEL;
}

<tree>{LABELCH}+ {
	strncpy(yylval.cs, yytext, 1023);
	yylval.cs[1023] = '\0';
	return LABEL;	
}
<tree>"[".+"]" { }

<*><<EOF>> {
	yyterminate();
	return END;
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

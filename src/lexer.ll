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

bool g_bParseOkay = true;

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
IDWORD [A-Za-z][A-Za-z_0-9]*
STR    \"[^\"\n]*\"
LABELCH [^ \t\n\r\v\f\(\)\[\]:;,\'\"]
NUMBER [-+]?{DIGIT}+("."{DIGIT}+)?([eE][+-]?{DIGIT}+)?
SPACE [ \t\r\v\f]

%x tree
%x tostr

%%

[=\{\}] {
	yylval.ch = yytext[0];
	return yytext[0];
}

[Ff]"alse" {
	yylval.b = false;
	return BOOL;
}

[Tt]"rue" {
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

{STR} {
	yytext[strlen(yytext)-1] = '\0';
	yylval.pss = new string(yytext+1);
	return STRING;
}

"<<"{IDWORD}{SPACE}*"\n" {
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

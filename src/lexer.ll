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

"<<"{IDWORD}{SPACE}+"\n" {
	yytext += 2;
	for(int i=0;!isspace(*(yytext+i));++i) { }
	yytext[i] = '\n';
	yytext[i+1] = '\0';
	
	yylval.pss = new string;
	string ssTemp;
	while(1)
	{
		int c = yyinput();
		ssTemp += c;
		if(c == '\n')
		{
			if(ssTemp == yytext)
				break;
			yylval.pss->append(ssTemp);
			ssTemp.clear();
		}
		else if(c == EOF)
		{
			yylval.pss->append(ssTemp);
			ssTemp.clear();
			break;			
		}
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

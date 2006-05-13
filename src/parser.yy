%{
// parser.yy - Copyright (C) 2004-2005 Reed A. Cartwright (all rights reserved)

#include "dawg.h"
#include "var.h"

#define YYERROR_VERBOSE 1

extern char yytext[];
extern FILE *yyin;
int yylex(void);
void yyerror (char *s);

using namespace std;

string g_ssSection("");
extern VarDB *g_pDB;

#ifdef _MSC_VER
#	pragma warning(disable: 4065 4244 4127 4102 4706)
#endif

string varName(const char *cs)
{
	string ss = g_ssSection;
	if(!ss.empty())
		ss += ".";
	ss += cs;
	return ss;
}

%}

%union {
	double d;	/* number values */
	char   *cs;  /* string values */
	char   ch;  /* characters */
	bool   b;   /* booleans */
	NewickNode	*pnode; /*Tree*/
	Variable *pvar;
}

%token <d>  NUM
%token <d>  LENGTH
%token <cs> STRING
%token <cs> LABEL
%token <cs> BID
%token <cs> ID
%token <b>  BOOL
%token <ch> CHAR
%token <ch> EQ		'='
%token <ch> QEQ		'?'
%token <ch> AEQ		'+'
%token <ch> LBRACE	'{'
%token <ch> RBRACE	'}'
%token <ch> LPARTH	'('
%token <ch> RPARTH	')'
%token <ch> TO		'<'
%token <ch> LBRACKET '['
%token <ch> RBRACKET ']'
%token <ch> ENDL	 '\n'
%token <ch> SQUOTE	
%token <ch> DQUOTE	
%token <ch> UNKNOWN
%token      END

%type <pvar> dvar
%type <pvar> vseq
%type <pnode> tree
%type <pnode> nodeseq
%type <pnode> node
%type <cs>	 chseq
%type <cs>  qstring 


%expect 1

%%

input:
	/* empty */
| input statement
;

statement:
'['']' { g_ssSection = ""; }
| '[' ID ']' {
	if($2[0] == '.') {
		g_ssSection.append($2);
	} else
		g_ssSection = $2;
	free($2);
}
| BID '=' dvar { g_pDB->SetVar(varName($1), $3, 0); free($1); }
| BID '?' dvar { g_pDB->SetVar(varName($1), $3, 1); free($1); }
| BID '+' dvar { g_pDB->SetVar(varName($1), $3, 2); free($1); }
;

dvar:
'{' vseq '}' { $$ = $2; }
| NUM     { $$ = static_cast<Variable*>(new NumberVar($1)); }
| BOOL    { $$ = static_cast<Variable*>(new BooleanVar($1)); }
| STRING  { $$ = static_cast<Variable*>(new StringVar($1)); free($1); }
| qstring { $$ = static_cast<Variable*>(new StringVar($1)); free($1); }
| tree    { $$ = static_cast<Variable*>(new TreeVar($1)); }
;

qstring:
DQUOTE chseq DQUOTE { $$ = $2; }
| SQUOTE chseq SQUOTE { $$ = $2; }
;

chseq:
CHAR { $$ = new char[1024]; $$[0] = $1; $$[1] = '\0'; }
| chseq CHAR {
	size_t len = strlen($1);
	size_t i;
	for(i = 1024; i < len+1; i*= 2) { }	
	if(len+1 == i)
	{
		$$ = new char[2*i];
		strcpy($$, $1);
		free($1);
		break;
	}
	else
	{
		$$ = $1;
	}
	$$[len] = $2;
	$$[len+1] = '\0';
}
;

vseq:
  dvar { $$ = static_cast<Variable *>(new VectorVar($1)); }
| vseq dvar { $$ = $1; $$->Append($2); }
;

tree:
node
;

node:
  '(' nodeseq ')' LABEL LENGTH { $$ = new NewickNode($2, $4, $5 ); free($4); }
| '(' nodeseq ')' LABEL { $$ = new NewickNode($2, $4, 0.0); free($4); }
| '(' nodeseq ')' LENGTH { $$ = new NewickNode($2, NULL, $4); }
| '(' nodeseq ')' {	$$ = new NewickNode($2, NULL, 0.0); }
| LABEL LENGTH { $$ = new NewickNode(NULL, $1, $2); free($1); }
| LABEL { $$ = new NewickNode(NULL, $1, 0.0); free($1); }
;

nodeseq:
nodeseq node { $$ = $2; $$->m_pSib.reset($1); }
| node
;

%%



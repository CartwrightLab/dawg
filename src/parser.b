%{
#pragma warning(disable: 4244 4065 4255 4127 4102)
#include "dawg.h"
#include "var.h"


#define YYERROR_VERBOSE 1

extern char yytext[];
extern FILE *yyin;
int yylex(void);
void yyerror (char *s);

using namespace std;
%}

%union {
	double d;	/* number values */
	char  cs[1024];  /* string values */
	char   ch;  /* characters */
	bool   b;   /* booleans */
	DawgVar::Vec *pvec; /*vector*/
	DawgVar *pvar; /*DawgVar*/
	NewickNode	*pnode; /*Tree*/
	std::string* pstr;
}

%token <d>  NUM
%token <d>  LENGTH
%token <cs> STRING
%token <cs> LABEL
%token <cs> ID
%token <b>  BOOL
%token <ch> DOT		'.'
%token <ch> EQ		'='
%token <ch> LBRACE	'{'
%token <ch> RBRACE	'}'
%token <ch> LPARTH	'('
%token <ch> RPARTH	')'
%token <ch> UNKNOWN
%token      END

%type <pstr> str
%type <pvar> dvar
%type <pvec> vvector
%type <pvec> vseq
%type <pnode> tree
%type <pnode> nodeseq
%type <pnode> node


%expect 1

%%

input:
	/* empty */
| input statement
;

statement:
END { }
| ID '=' dvar { DawgVar::SetVar(string($1), $3); }
;

str:
str '.' STRING { $$ = $1; $$->append("\n");	$$->append($3); }
| STRING { $$ = new std::string($1); }
;

dvar:
  vvector { $$ = new DawgVar($1); }
| NUM { $$= new DawgVar($1); }
| BOOL { $$= new DawgVar($1); }
| str {	$$ = new DawgVar(*$1); delete $1; }
| tree { $$ = new DawgVar($1); }
;

vvector: '{' vseq '}' { $$ = $2; }
;

vseq:
  dvar { $$ = new DawgVar::Vec; $$->push_back($1); }
| vseq dvar { $$ = $1; $$->push_back($2); }
;

tree:
node
;

node:
  '(' nodeseq ')' LABEL LENGTH { $$ = new NewickNode($2, $4, $5 ); }
| '(' nodeseq ')' LABEL { $$ = new NewickNode($2, $4, 0.0); }
| '(' nodeseq ')' LENGTH { $$ = new NewickNode($2, NULL, $4); }
| '(' nodeseq ')' {	$$ = new NewickNode($2, NULL, 0.0); }
| LABEL LENGTH { $$ = new NewickNode(NULL, $1, $2); }
| LABEL { $$ = new NewickNode(NULL, $1, 0.0); }
;

nodeseq:
nodeseq node { $$ = $2; $$->m_pSib.reset($1); }
| node
;

%%



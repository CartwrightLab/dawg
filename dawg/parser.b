%{
// cmd: bison -d parser.b -o parser.c

#include "dawg.h"
#include <vector>
#include <stdio.h>

#pragma warning(disable: 4244 4702)

#define YYERROR_VERBOSE 1

using namespace std;

extern char yytext[];
extern FILE *yyin;
int yylex(void);
void yyerror (char *s);
%}

%union {
	double d;	/* number values */
	char*  cs;  /* string values */
	char   ch;  /* characters */
	bool   b;   /* booleans */
	DawgVar::Vec *pvec; /*vector*/
	DawgVar *pvar; /*DawgVar*/
	Node	*pnode; /*Tree*/
}

%token <d>  NUM
%token <cs> STRING
%token <ch> DOT
%token <cs> LABEL
%token <ch> EQ
%token      END
%token <cs> ID
%token <b>  BOOL
%token <ch> COLON
%token <ch> LBRACE
%token <ch> RBRACE
%token <ch> LPARTH
%token <ch> RPARTH
%token <ch> UNKNOWN

%type <cs>	  label
%type <cs>	  str
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

END
{

}

| ID EQ dvar
{
	DawgVar::SetVar(string($1), $3);
	delete[] $1;
}
;

str:
str DOT STRING {
	size_t t1 = strlen($1);
	size_t t2 = strlen($3);
	$$ = strcpy(new char[t1+t2+2], $1);
	strcat($$, "\n");
	strcat($$, $3);
	delete[] $1;
	delete[] $3;
}
| STRING { $$ = $1; }
;

dvar:

vvector
{
	$$ = new DawgVar($1);
}
| NUM { $$= new DawgVar($1); }
| BOOL { $$= new DawgVar($1); }
| str
{
	$$= new DawgVar(string($1));
	delete[] $1;
}
| tree { $$ = new DawgVar($1); }
;

vvector: LBRACE vseq RBRACE { $$ = $2; };

vseq:
vseq dvar
{
	$$->push_back($2);
}

| dvar
{
	$$ = new DawgVar::Vec;
	$$->push_back($1);
};

label:
LABEL { $$ = $1; };

tree: node { $$ = $1; };

node:
LPARTH nodeseq RPARTH label COLON NUM
{
	$$ = new Node($4, $2);
	$$->BranchLength($6);
	delete[] $4;	
}

| LPARTH nodeseq RPARTH label
{
	$$ = new Node($4, $2);
	delete[] $4;	
}

| LPARTH nodeseq RPARTH COLON NUM
{
	$$ = new Node("", $2);
	$$->BranchLength($5);
}

| LPARTH nodeseq RPARTH
{
	$$ = new Node("", $2);
}

| label COLON NUM
{
	$$ = new Node($1);
	$$->BranchLength($3);
	delete[] $1;		
}

| label
{
	$$ = new Node($1);
	delete[] $1;		
};

nodeseq:
nodeseq node
{
	$1->AddSib($2);
}

| node
{
	$$ = $1;
};

%%



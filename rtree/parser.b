%{
#pragma warning(disable: 4244 4065 4255 4127 4102)

#include <stdio.h>
#include <string.h>
#include <iostream>
#include "rtree.h"

extern char yytext[];
extern FILE *yyin;
int yylex(void);
void yyerror (char *s);

using namespace std;

void PrintNode(Node* pNode)
{
	
	if(pNode->pSub)
	{
		cout << "(";
		PrintNode(pNode->pSub);
		cout << ")";
	}
	cout << pNode->csLabel;
	cout << ":";
	cout << pNode->dLen;
	if(pNode->pSib)
	{
		cout << ",";
		PrintNode(pNode->pSib);
	}

}

%}

%union {
	double d;	/* number values */
	char  cs[1024];  /* string values */
	Node  *node;
	char ch;
}

%token <d>  LENGTH
%token <cs> LABEL
%token      END
%token      UNKNOWN
%token <ch> LPARTH '('
%token <ch> RPARTH ')'
%token <ch> COLON ':'
%token <ch> SEMICOLON ';'

%type <node> node
%type <node> nodeseq

%start input

%expect 1

%%

input:
input tree
| tree
;

tree:
node ';' { Node::s_stack.push_back($1); }
| END
;

node:
'(' nodeseq ')' LABEL LENGTH {	$$ = new Node($4, $5, $2); }
| '(' nodeseq ')' LABEL { $$ = new Node($4, 0.0, $2); }
| '(' nodeseq ')' LENGTH { $$ = new Node(NULL, $4, $2); }
| '(' nodeseq ')' {$$ = new Node(NULL, 0.0, $2); }
| LABEL LENGTH { $$ = new Node($1, $2, NULL); }
| LABEL { $$ = new Node($1, 0.0, NULL); }
;

nodeseq:
nodeseq node
{
	$$ = $2;
	$$->pSib = $1;
}
| node
{
	$$ = $1;
}
;

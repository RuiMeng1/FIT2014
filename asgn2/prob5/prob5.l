 /* Author : Rui Meng (33095612) */

/* Definition section */
%{
/* C code to be copied verbatim */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#undef  I
void yyerror(char *);
#include "quant.h"
#include "y.tab.h"
%}

%%
    /* Rules section */

(I|H|X|Y|Z|CNOT|TOF) {
    yylval.str = strdup(yytext);

    if (strcmp(yytext, "I") == 0) return I;
    if (strcmp(yytext, "H") == 0) return H;
    if (strcmp(yytext, "X") == 0) return X;
    if (strcmp(yytext, "Y") == 0) return Y;
    if (strcmp(yytext, "Z") == 0) return Z;
    if (strcmp(yytext, "CNOT") == 0) return CNOT;
    if (strcmp(yytext, "TOF") == 0) return TOF;
}
    /* TODO: Do we need to seperate this into individual statements?? */
[*\(\)]   {
            return  *yytext;
        }
"(x)" {
    // printf("Token: KRONECKERPROD;  Lexeme: %s\n", yytext);
    return KRONECKERPROD;
}
\n      {
    //printf("Token and Lexeme: <newline>\n");
    return *yytext;
    }
[ \t]   {    }     /*  skip whitespace  */

.       {
           yyerror("invalid character");

        }

%%
    /* C Code section */

void yyerror(char *s) {
      fprintf(stderr, "line %d: yytext = %s.  Error msg: %s.\n", yylineno, yytext, s);
}


int yywrap(void) {
   return 1;
}


/*  comment out the function  main()  if using  lex  with  yacc   */
// int main(void) {
//     yylex();
//     return 0;
// }
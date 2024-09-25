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
/*
#include "y.tab.h"
*/
%}

%%
    /* Rules section */

(I|H|X|Y|Z|CNOT|TOF) {
    // printf("Token: %s;  Lexeme: %s\n", yytext, yytext);

    // return token
    return *yytext
}
[*()] {
    // printf("Token and Lexeme: %s\n", yytext);
    return *yytext
}
"(x)" {
    // printf("Token: KRONECKERPROD;  Lexeme: %s\n", yytext);
    return KRONECKERPROD
}
\n      {
    //printf("Token and Lexeme: <newline>\n");
    return *yytext
    }
[ \t]   {    }     /*  skip whitespace  */

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
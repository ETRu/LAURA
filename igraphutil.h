
#pragma once

#if defined (WIN32)
#include "igraph/igraph.h"
#else
#include "igraph.h"
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <locale.h>

/*
 * ------------------------- PRINTING AND UTILITIES -----------------------
 *
 */


/*This function PRINTS A IGRAPH STRUCT VECTOR*/
void strvector_print(const igraph_strvector_t *sv);

/*This function PRINTS A IGRAPH VECTOR as a COLUMN*/
void print_vector_column(igraph_vector_t *v, FILE *f);

/*This function PRINTS A IGRAPH VECTOR as a LINE*/
void print_vector_line(igraph_vector_t *v, FILE *f);

/*This function PRINTS A IGRAPH MATRIX ROW as a LINE*/
void print_matrix_row_line(igraph_matrix_t *m, int row, FILE *f);

/*This function PRINTS A IGRAPH MATRIX COLUMN as a LINE*/
void print_matrix_col_line(igraph_matrix_t *m, int col, FILE *f);

/*This function PRINTS A IGRAPH MATRIX*/
void print_matrix(const igraph_matrix_t *m, FILE *f);

/*This function PRINTS A _USAER EASILY READABLE_ IGRAPH MATRIX*/
void print_matrix_ur(const igraph_matrix_t *m, FILE *f);

/*This function PRINTS A INDEXED IGRAPH VECTOR*/
void print_vector_indexed(igraph_vector_t *v, FILE *f);



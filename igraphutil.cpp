
#include "igraphutil.h"




/*
 * ------------------------- PRINTING AND UTILITIES -----------------------
 *
 */


/*This function PRINTS A IGRAPH STRUCT VECTOR*/
void strvector_print(const igraph_strvector_t *sv) {
    long int i, s=igraph_strvector_size(sv);
    for (i=0; i<s; i++) {
        printf("---%s---\n", STR(*sv, i));
    }
    
}

/*This function PRINTS A IGRAPH VECTOR as a COLUMN*/
void print_vector_column(igraph_vector_t *v, FILE *f) {
    long int i;
    for (i=0; i<igraph_vector_size(v); i++) {
        fprintf(f, "%f \n", (double)VECTOR(*v)[i]);
    }
    fprintf(f, "\n");
}


/*This function PRINTS A IGRAPH VECTOR as a LINE*/
void print_vector_line(igraph_vector_t *v, FILE *f) {
    long int i;
    for (i=0; i<igraph_vector_size(v); i++) {
        fprintf(f, "%f ", (double)VECTOR(*v)[i]);
    }
    fprintf(f, "\n");
}

/*This function PRINTS A IGRAPH MATRIX ROW as a LINE*/
void print_matrix_row_line(igraph_matrix_t *m, int row, FILE *f) {
    long int i;
    for (i=0; i<igraph_matrix_ncol(m); i++) {
        fprintf(f, "%f ", (double)MATRIX(*m,row,i));
    }
    fprintf(f, "\n");
}

/*This function PRINTS A IGRAPH MATRIX COLUMN as a LINE*/
void print_matrix_col_line(igraph_matrix_t *m, int col, FILE *f) {
    long int i;
    for (i=0; i<igraph_matrix_nrow(m); i++) {
        fprintf(f, "%f ", (double)MATRIX(*m,i,col));
    }
    fprintf(f, "\n");
}




/*This function PRINTS A IGRAPH MATRIX*/
void print_matrix(const igraph_matrix_t *m, FILE *f) {
    long int i, j, nrow=igraph_matrix_nrow(m), ncol=igraph_matrix_ncol(m);
    for (i=0; i<nrow; i++) {
        for (j=0; j<ncol; j++) {
            fprintf(f,"%lf", (double)MATRIX(*m, i, j));
            if (j!=ncol-1) { fprintf(f," "); }
        }
        fprintf(f,"\n");
    }
}




/*This function PRINTS A _USAER EASILY READABLE_ IGRAPH MATRIX*/
void print_matrix_ur(const igraph_matrix_t *m, FILE *f) {
    long int i, j, nrow=igraph_matrix_nrow(m), ncol=igraph_matrix_ncol(m);
    for (i=0; i<nrow; i++) {
        for (j=0; j<ncol; j++) {
//            fprintf(f,"%.2g", (double)MATRIX(*m, i, j));
            fprintf(f,"%.3f", (double)MATRIX(*m, i, j));
            if (j!=ncol-1) { fprintf(f," "); }
        }
        fprintf(f,"\n");
    }
}



/*This function PRINTS A INDEXED IGRAPH VECTOR*/
void print_vector_indexed(igraph_vector_t *v, FILE *f) {
    long int i;
    for (i=0; i<igraph_vector_size(v); i++) {
        fprintf(f, "%li %f \n", i, (double)VECTOR(*v)[i]);
    }
    fprintf(f, "\n");
}



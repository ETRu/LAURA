
#include "networkop.h"


#ifndef PI
#define PI 3.141593
#endif

extern int nodesnumber;
extern int sourcenode;
extern igraph_matrix_t admatrix;
extern igraph_t graph;
extern igraph_t sgraph;
extern igraph_matrix_t layout;
extern igraph_matrix_t density, densityold, state, statenew, loadedstate;
extern igraph_matrix_t flux;
extern igraph_matrix_t bridgeslinks;
extern igraph_vector_t gain;
extern igraph_matrix_t loss;
extern igraph_vector_t statstate;
extern igraph_matrix_t estates;

extern int totrun;

extern int framexmin, framexmax, frameymin, frameymax;

extern char errorstring[100];

//----- CONTROL VARIABLES - NEEDED FROM MAIN.CPP
extern int rewrite;

extern int havepath;
extern int graphisloaded;
extern int error;

extern int usesteady;
extern int haveloadedstate;
extern int isdissipating;


void adjustlayout(){
    
    igraph_vector_t coords;
    igraph_vector_init(&coords,0);
    
    
    //arrange x
    igraph_matrix_get_col(&layout,&coords,0);
    //set positive
    //if(igraph_vector_min(&coords)<0)
    {igraph_vector_add_constant(&coords,(-1.)*igraph_vector_min(&coords));}
    //set 0
    
    igraph_vector_scale(&coords,(framexmax-framexmin-(2*SHIFT))/(igraph_vector_max(&coords)));
    igraph_matrix_set_col(&layout,&coords,0);
    
    //arrange y
    igraph_matrix_get_col(&layout,&coords,1);
    //set positive
    //if(igraph_vector_min(&coords)<0)
    {igraph_vector_add_constant(&coords,(-1.)*igraph_vector_min(&coords));}
    igraph_vector_scale(&coords,(frameymax-frameymin-(2*SHIFT))/(igraph_vector_max(&coords)));
    igraph_matrix_set_col(&layout,&coords,1);
    
    //shift
    igraph_matrix_add_constant(&layout, SHIFT);
    
    igraph_vector_destroy(&coords);
    
}



//--------------------------------------------
void generatelattice(int mylatticedim, int mylatticeside, int myistoro, int myrand, int diss, double drate, int dnode) {
    
    igraph_vector_t row;
    init_genrand(0);
    
    igraph_vector_t dimvector;
    igraph_bool_t directed=1, mutual=1;
    
    igraph_vector_init(&dimvector,mylatticedim);
    igraph_vector_fill(&dimvector,mylatticeside);
    igraph_lattice(&graph, &dimvector, 1, directed, mutual, myistoro);
    
    
    igraph_matrix_init(&admatrix, 0, 0);
    igraph_get_adjacency(&graph, &admatrix,IGRAPH_GET_ADJACENCY_BOTH, 0);
    
    igraph_copy(&sgraph,&graph);
    igraph_simplify(&sgraph, 1,1,0);
    igraph_to_undirected(&sgraph, IGRAPH_TO_UNDIRECTED_COLLAPSE,/*edge_comb=*/ 0);
    
    nodesnumber=igraph_matrix_nrow(&admatrix);
    
    
    logDebug("\nnodesnumber=%i\n\n", nodesnumber);
    
    igraph_vector_destroy(&dimvector);
    
    //----------- prob matrix-----------------
    //ROW NORMALIZATION:
    igraph_vector_init(&row,nodesnumber);
    for (int i=0; i<nodesnumber; ++i) {
        igraph_matrix_get_row(&admatrix,&row,i);
        
        //set random connectivities
        if(myrand==1){
            for(int j=0; j<nodesnumber;++j){
                if(VECTOR(row)[j]!=0){
                    VECTOR(row)[j]=genrand_real1();
                }
            }
        }
        
        igraph_vector_scale(&row, 1./(igraph_vector_sum(&row)));
        igraph_matrix_set_row(&admatrix, &row,i);
        
    }
    igraph_vector_destroy(&row);
    
    
    isdissipating=0;
    sourcenode=0;
    
    if(diss==1 && drate>0){
        
        isdissipating=1;
        sourcenode=dnode;
        
        //add the dissipation node to admatrix:
        igraph_matrix_add_rows(&admatrix, 1);
        for(int i=0; i<igraph_matrix_ncol(&admatrix);++i){MATRIX(admatrix,igraph_matrix_nrow(&admatrix)-1,i)=0;}
        igraph_matrix_add_cols(&admatrix, 1);
        for(int i=0; i<igraph_matrix_nrow(&admatrix);++i){MATRIX(admatrix,i,igraph_matrix_ncol(&admatrix)-1)=0;}
        
        printf("\n clean version \n");
        //print_matrix_ur(&admatrix,stdout);printf("\n\n\n\n");
        
        //connect the sink
        for(int i=0;i<nodesnumber; ++i){
            
           /* double totsum=0;
            for(int j=0;j<igraph_matrix_ncol(&admatrix);++j){
                totsum=totsum+MATRIX(admatrix,i,j);
            }
            
            //rescale the previous connections
            for(int j=0; j<igraph_matrix_ncol(&admatrix)-1;++j){
                double tempval;
                tempval=MATRIX(admatrix,i,j);
                tempval=tempval-(tempval*drate);
                MATRIX(admatrix,i,j)=tempval;
            }
            
            //add dissipation
            */
            
            MATRIX(admatrix,i,nodesnumber)=drate;
            
        }
        
        //connect the source (dissipation node---->cols of the source)
        //connect the sink (rows of the sink cluster--->dissipation node)
        if(dnode>nodesnumber){MATRIX(admatrix,nodesnumber,((nodesnumber/2)))=1;}
        else{MATRIX(admatrix,nodesnumber,dnode)=1;}
        
        
        printf("\n connected version \n");
       // print_matrix_ur(&admatrix,stdout);printf("\n\n\n\n");
        
        
        //generate graph from admatrix
        igraph_weighted_adjacency(&graph, &admatrix,IGRAPH_ADJ_UPPER, "w",0);
        
        
        nodesnumber=igraph_matrix_nrow(&admatrix);
        
        
        logDebug("\n\nnodesnumber=%i\n\n", nodesnumber);
        
        //print_matrix(&admatrix, stdout);
        
        igraph_copy(&sgraph,&graph);
        igraph_simplify(&sgraph, 1,1,0);
        igraph_to_undirected(&sgraph, IGRAPH_TO_UNDIRECTED_COLLAPSE,/*edge_comb=*/ 0);
        
        

        
        //----------- prob matrix-----------------
        //ROW NORMALIZATION:
        igraph_vector_init(&row,nodesnumber);
        for (int i=0; i<nodesnumber; ++i) {
            igraph_matrix_get_row(&admatrix,&row,i);
            
            igraph_vector_scale(&row, 1./(igraph_vector_sum(&row)));
            igraph_matrix_set_row(&admatrix, &row,i);
            
        }
        igraph_vector_destroy(&row);
        
        
        //print_matrix_ur(&admatrix,stdout);

        
        
        
    }
    
    
    //print_matrix_ur(&admatrix,stdout);
    
    //layout
    
    igraph_matrix_init(&layout, 0, 0);
    igraph_layout_grid(&graph, &layout, mylatticeside);
    
    adjustlayout();
    
    error=0;
    
    haveloadedstate=0;
}



//--------------------------------------------
void generatepotlat2d(int mylatticeside, int myistoro, int diss, double drate, int dnode){
    
    double xmax, xmin, ymax, ymin;
    double stepx, stepy;
    igraph_vector_t neis;
    float pot1, pot2;
    
    igraph_vector_t row;
    
    igraph_matrix_t latticemat;
    
    igraph_vector_t dimvector;
    igraph_bool_t directed=1, mutual=1;
    
    igraph_vector_init(&dimvector,2);
    igraph_vector_fill(&dimvector,mylatticeside);
    igraph_lattice(&graph, &dimvector, 1, directed, mutual, myistoro);
    
    
    igraph_matrix_init(&admatrix, 0, 0);
    igraph_get_adjacency(&graph, &admatrix,IGRAPH_GET_ADJACENCY_BOTH, 0);
    
    igraph_copy(&sgraph,&graph);
    igraph_simplify(&sgraph, 1,1,0);
    igraph_to_undirected(&sgraph, IGRAPH_TO_UNDIRECTED_COLLAPSE,/*edge_comb=*/ 0);
    
    nodesnumber=igraph_matrix_nrow(&admatrix);
    
    
    logDebug("\nnodesnumber=%i\n\n", nodesnumber);
    
    igraph_vector_destroy(&dimvector);
    
    
    
    //---------- SET POTENTIAL -----------
    
    
    //set admatrix to zero.
    igraph_matrix_null(&admatrix);
    
    //generate the matrix that represents the lattice.
    igraph_matrix_init(&latticemat,mylatticeside,mylatticeside);
    
    //set maxs and mins, steps
    xmin=-1; xmax=1;
    ymin=-1; ymax=1;
    
    stepx=(xmax-xmin)/(mylatticeside-1);
    stepy=(ymax-ymin)/(mylatticeside-1);
    
    
    // by now, we will define statically the function vi.
    //
    // in this case I choose V(x,y)=x^2 + y^2
    
    // printf("\n");
    for(int xcoord=0; xcoord<mylatticeside; ++xcoord){
        for(int ycoord=0; ycoord<mylatticeside; ++ycoord){
            //printf("xmin=%f stepx=%f xcoord=%i=> x=%f; ymin=%f stepy=%f ycoord=%i =>y=%f   ", xmin, stepx, xcoord, xmin+xcoord*stepx, ymin, stepy, ycoord,ymin+ycoord*stepy);
            
            
            //V(x,y)=x^2 + y^2
            MATRIX(latticemat,xcoord,ycoord)=(xmin+xcoord*stepx)*(xmin+xcoord*stepx) + (ymin+ycoord*stepy)*(ymin+ycoord*stepy)  ;
            
            //V(x,y)=cos(x*2*PI)+y^2;
            
            MATRIX(latticemat,xcoord,ycoord)=3*cos((xmin+xcoord*stepx)*2*PI)+(ymin+ycoord*stepy)*(ymin+ycoord*stepy);
            
            
            //MATRIX(latticemat,xcoord,ycoord)=0;
            
            //printf(" RES: %f\n",MATRIX(latticemat,xcoord,ycoord) );
        }
        
        //printf("\n");
        
    }
    
    
    
    
    printf("\n LATTICEMAT \n");
    print_matrix(&latticemat,stdout);
    printf("\n");
    
    //set p_ij
    igraph_vector_init(&neis, 0);
    
    for(int knot=0; knot<nodesnumber; ++knot){
        
        //pick its neighbours
        igraph_neighbors(&graph, &neis, knot, IGRAPH_OUT);
        
        //for each neighbour
        for(int neigh=0;neigh<igraph_vector_size(&neis);++neigh){
            //set admatrix
            int knotx, knoty, neighx, neighy;
            
            knotx=knot%mylatticeside;
            knoty=knot/mylatticeside;
            neighx=((int)VECTOR(neis)[neigh])%mylatticeside;
            neighy=((int)VECTOR(neis)[neigh])/mylatticeside;
            
            pot1=MATRIX(latticemat,knotx,knoty);
            pot2=MATRIX(latticemat,neighx,neighy);
            
            int a; int b;
            a=knot; b=VECTOR(neis)[neigh];
            MATRIX(admatrix,a,b)= exp((pot1-pot2)/2);
            
            //printf ("a:%i (lattice:%i %i) b:%i (lattice:%i %i) pot1=%f pot2=%f  admat(a,b)=%f\n",a,knotx,knoty, b, neighx, neighy, pot1,pot2,MATRIX(admatrix,a,b));
            
            
            
        }
        
    }
    
    
    //----- steady state
    
    igraph_vector_init(&statstate,nodesnumber);
    
    for (int knot=0; knot<nodesnumber; ++knot) {
        
        int knotx, knoty, neighx, neighy;
        
        knotx=knot%mylatticeside;
        knoty=knot/mylatticeside;
        
        pot1=MATRIX(latticemat,knotx,knoty);
        VECTOR(statstate)[knot]=exp(-pot1);
    }
    
    igraph_vector_scale(&statstate,1./igraph_vector_sum(&statstate));
    
    usesteady=1;
    
    
    //
    
    
    
    igraph_vector_destroy(&neis);
    
    igraph_matrix_destroy(&latticemat);
    
    //----------- prob matrix-----------------
    //ROW NORMALIZATION:
    
    igraph_vector_init(&row,nodesnumber);
    for (int i=0; i<nodesnumber; ++i) {
        igraph_matrix_get_row(&admatrix,&row,i);
        
        
        
        igraph_vector_scale(&row, 1./(igraph_vector_sum(&row)));
        igraph_matrix_set_row(&admatrix, &row,i);
        
    }
    igraph_vector_destroy(&row);
    
    
    /* printf("\n ADMATRIX NORMALIZED \n");
     print_matrix_ur(&admatrix,stdout);
     printf("\n");
     */
    
    //print_matrix_ur(&admatrix,stdout);
    
    
    
    isdissipating=0;
    sourcenode=0;
    
    
    if(diss==1 && drate>0){
        isdissipating=1;
        sourcenode=dnode;
        
        //add the dissipation node to admatrix:
        igraph_matrix_add_rows(&admatrix, 1);
        for(int i=0; i<igraph_matrix_ncol(&admatrix);++i){MATRIX(admatrix,igraph_matrix_nrow(&admatrix)-1,i)=0;}
        igraph_matrix_add_cols(&admatrix, 1);
        for(int i=0; i<igraph_matrix_nrow(&admatrix);++i){MATRIX(admatrix,i,igraph_matrix_ncol(&admatrix)-1)=0;}
        
        printf("\n clean version \n");
        //print_matrix_ur(&admatrix,stdout);printf("\n\n\n\n");
        
        //connect the sink
        for(int i=0;i<nodesnumber-1; ++i){
            
          /*  double totsum=0;
            for(int j=0;j<igraph_matrix_ncol(&admatrix);++j){
                totsum=totsum+MATRIX(admatrix,i,j);
            }
            
            //rescale the previous connections
            for(int j=0; j<igraph_matrix_ncol(&admatrix)-1;++j){
                double tempval;
                tempval=MATRIX(admatrix,i,j);
                tempval=tempval-(tempval*drate);
                MATRIX(admatrix,i,j)=tempval;
            }
            
            //add dissipation
            */
            MATRIX(admatrix,i,nodesnumber)=drate;
            
        }
        
        //connect the source (dissipation node---->cols of the source)
        //connect the sink (rows of the sink cluster--->dissipation node)
        if(dnode>nodesnumber){MATRIX(admatrix,nodesnumber,((nodesnumber/2)))=1;}
        else{MATRIX(admatrix,nodesnumber,dnode)=1;}
        
        
        printf("\n connected version \n");
        //print_matrix_ur(&admatrix,stdout);printf("\n\n\n\n");
        
        
        //generate graph from admatrix
        igraph_weighted_adjacency(&graph, &admatrix,IGRAPH_ADJ_UPPER, "w",0);
        
        
        nodesnumber=igraph_matrix_nrow(&admatrix);
        
        
        logDebug("\n\nnodesnumber=%i\n\n", nodesnumber);
        
        //print_matrix(&admatrix, stdout);
        
        igraph_copy(&sgraph,&graph);
        igraph_simplify(&sgraph, 1,1,0);
        igraph_to_undirected(&sgraph, IGRAPH_TO_UNDIRECTED_COLLAPSE,/*edge_comb=*/ 0);
        

        
        //----------- prob matrix-----------------
        //ROW NORMALIZATION:
        igraph_vector_init(&row,nodesnumber);
        for (int i=0; i<nodesnumber; ++i) {
            igraph_matrix_get_row(&admatrix,&row,i);
            
            igraph_vector_scale(&row, 1./(igraph_vector_sum(&row)));
            igraph_matrix_set_row(&admatrix, &row,i);
            
        }
        igraph_vector_destroy(&row);
        
        
        //print_matrix_ur(&admatrix,stdout);
        
        
        
        
    }
    
    
    
    
    //layout
    
    igraph_matrix_init(&layout, 0, 0);
    igraph_layout_grid(&graph, &layout, mylatticeside);
    
    adjustlayout();
    
    error=0;
    
    haveloadedstate=0;
    
    
    
}










//--------------------------------------------
int generaterandom1(int n, float probability, int diss, double drate, int dnode) {
    
    igraph_t graphtemp;
    
    int central;
    
    igraph_vector_t row;
    igraph_bool_t connected;
    igraph_vector_t res;
    
    
    for(int i=0;i<1000;++i){
        //generate graph
        igraph_erdos_renyi_game(&graphtemp, IGRAPH_ERDOS_RENYI_GNP, n, probability, 1,0);
        //check if it's connected
        igraph_is_connected(&graphtemp, &connected,IGRAPH_STRONG);
        
        if(connected)break;
    }
    
    if(connected){
        igraph_copy(&graph, &graphtemp);
        logDebug("\nThe graph is connected.\n");
        graphisloaded=1;
    }
    
    else
    {
        logDebug("\nERROR: UNABLE TO GENERATE A CONNECTED GRAPH.\n");
        graphisloaded=0; error=1;  rewrite=1;
        sprintf(errorstring,"ERROR!\nNON\nCONNECTED\nGRAPH!\nGRAPH\nNOT\nLOADED\ntry again");
        return 0;
    }
    
    
    //generate adjacece matrix
    igraph_matrix_init(&admatrix, 0, 0);
    igraph_get_adjacency(&graph, &admatrix,IGRAPH_GET_ADJACENCY_BOTH, 0);
    
    igraph_copy(&sgraph,&graph);
    igraph_simplify(&sgraph, 1,1,0);
    igraph_to_undirected(&sgraph, IGRAPH_TO_UNDIRECTED_COLLAPSE,/*edge_comb=*/ 0);
    
    nodesnumber=igraph_matrix_nrow(&admatrix);
    
    logDebug("\nnodesnumber=%i\n\n", nodesnumber);
    
    
    //----------- prob matrix-----------------
    //ROW NORMALIZATION:
    igraph_vector_init(&row,nodesnumber);
    for (int i=0; i<nodesnumber; ++i) {
        igraph_matrix_get_row(&admatrix,&row,i);
        igraph_vector_scale(&row, 1./(igraph_vector_sum(&row)));
        igraph_matrix_set_row(&admatrix, &row,i);
    }
    igraph_vector_destroy(&row);
    
    
    
    
    isdissipating=0;
    sourcenode=0;
    
    
    if(diss==1 && drate>0){
        
        isdissipating=1;
        sourcenode=dnode;
        
        //add the dissipation node to admatrix:
        igraph_matrix_add_rows(&admatrix, 1);
        for(int i=0; i<igraph_matrix_ncol(&admatrix);++i){MATRIX(admatrix,igraph_matrix_nrow(&admatrix)-1,i)=0;}
        igraph_matrix_add_cols(&admatrix, 1);
        for(int i=0; i<igraph_matrix_nrow(&admatrix);++i){MATRIX(admatrix,i,igraph_matrix_ncol(&admatrix)-1)=0;}
        
        printf("\n clean version \n");
        //print_matrix_ur(&admatrix,stdout);printf("\n\n\n\n");
        
        //connect the sink
        for(int i=0;i<nodesnumber-1; ++i){
            
            /*double totsum=0;
            for(int j=0;j<igraph_matrix_ncol(&admatrix);++j){
                totsum=totsum+MATRIX(admatrix,i,j);
            }
            
            //rescale the previous connections
            for(int j=0; j<igraph_matrix_ncol(&admatrix)-1;++j){
                double tempval;
                tempval=MATRIX(admatrix,i,j);
                tempval=tempval-(tempval*drate);
                MATRIX(admatrix,i,j)=tempval;
            }
            
            //add dissipation
            */
            MATRIX(admatrix,i,nodesnumber)=drate;
            
        }
        
        //connect the source (dissipation node---->cols of the source)
        //connect the sink (rows of the sink cluster--->dissipation node)
        if(dnode>nodesnumber){MATRIX(admatrix,nodesnumber,((nodesnumber/2)))=1;}
        else{MATRIX(admatrix,nodesnumber,dnode)=1;}
        
        
        printf("\n connected version \n");
        //print_matrix_ur(&admatrix,stdout);printf("\n\n\n\n");
        
        
        //generate graph from admatrix
        igraph_weighted_adjacency(&graph, &admatrix,IGRAPH_ADJ_UPPER, "w",0);
        
        
        nodesnumber=igraph_matrix_nrow(&admatrix);
        
        
        logDebug("\n\nnodesnumber=%i\n\n", nodesnumber);
        
        //print_matrix(&admatrix, stdout);
        
        igraph_copy(&sgraph,&graph);
        igraph_simplify(&sgraph, 1,1,0);
        igraph_to_undirected(&sgraph, IGRAPH_TO_UNDIRECTED_COLLAPSE,/*edge_comb=*/ 0);
        
        

        
        //----------- prob matrix-----------------
        //ROW NORMALIZATION:
        igraph_vector_init(&row,nodesnumber);
        for (int i=0; i<nodesnumber; ++i) {
            igraph_matrix_get_row(&admatrix,&row,i);
            
            igraph_vector_scale(&row, 1./(igraph_vector_sum(&row)));
            igraph_matrix_set_row(&admatrix, &row,i);
            
        }
        igraph_vector_destroy(&row);
        
        
        //print_matrix_ur(&admatrix,stdout);
        
        
        
        
    }
    
    
    
    
    
    //generate layout
    igraph_matrix_init(&layout, 0, 0);
    
    
    
    igraph_vector_init(&res,0);
    igraph_degree(&graph, &res, igraph_vss_all(), IGRAPH_ALL, IGRAPH_NO_LOOPS);
    central=igraph_vector_which_max(&res);
    igraph_vector_destroy(&res);
    igraph_layout_star(&graph, &layout,central,0);
    
    
    adjustlayout();
    
    igraph_destroy(&graphtemp);
    
    
    error=0;
    haveloadedstate=0;
    
    return central;
    
}




void generateclustsym(int clusnumber, int clusdim, double interconn, double intraconn, double wall){
    
    //i work on the ipper triangular diagonal.
    //i mainly work on the admatrix.
    //at the end, i symmetrize.
    
    init_genrand(0);
    
    igraph_matrix_t adtemp;
    igraph_t g;
    igraph_bool_t connected;
    
        igraph_vector_t row;
    
    igraph_matrix_init(&adtemp,clusnumber*clusdim,clusnumber*clusdim);
    igraph_matrix_null(&adtemp);
    
    
    //generate clusters
    for(int cluster=0; cluster<clusnumber; ++cluster){
        
        igraph_t gtemp;
        igraph_matrix_t adclus;
        
        igraph_erdos_renyi_game(&gtemp, IGRAPH_ERDOS_RENYI_GNP, clusdim, interconn, 0, 0);
        
        
        //ad clusters to adtemp
        igraph_matrix_init(&adclus,0,0);
        igraph_get_adjacency(&gtemp, &adclus, IGRAPH_GET_ADJACENCY_UPPER, 0);
        
        printf("\n cluster %i \n", cluster); //print_matrix_ur(&adclus,stdout); printf("\n\n");
        
        for (int j=0; j<clusdim; ++j) {
            for(int k=0; k<clusdim; ++k) {
                MATRIX(adtemp,(clusdim*cluster)+j,(clusdim*cluster)+k) = MATRIX(adclus,j,k);
                
            }
        }
        igraph_matrix_destroy(&adclus);
        igraph_destroy(&gtemp);
        
    }
    
    
    //link clusters
    
    //FROM the cluster:
    for(int clusfrom=0;clusfrom<clusnumber;++clusfrom){
        //TO every other cluster, upper triangular matrix
        for(int clusto=clusfrom+1;clusto<clusnumber;++clusto){
            
            //for each node in clusfrom
            for (int i=0; i<clusdim; ++i){
                //to each node in clusto
                for(int j=0; j<clusdim;++j){
                    
                    //throw coin
                    double coin;
                    coin=genrand_real1();
                    //if success
                    if (coin<intraconn) {
                        MATRIX(adtemp,(clusdim*clusfrom)+i,(clusdim*clusto)+j)=1./(wall+1);}
                    
                }}
            
        }
    }
    
    //print_matrix_ur(&adtemp,stdout);
    
    fflush(stdout);
    
    //generate graph from admatrix
    
    igraph_weighted_adjacency(&g, &adtemp,IGRAPH_ADJ_UPPER, "w",0);
    
    //check if it's connected
    igraph_is_connected(&g, &connected,IGRAPH_STRONG);
    
    if(connected){
        logDebug("\nThe network is connected.\n");
        graphisloaded=1;
    }
    
    else {
        logDebug("\nERROR: THE GRAPH GENERATED IS NOT CONNECTED.\n");
        graphisloaded=0; rewrite=1;
        sprintf(errorstring,"ERROR!\nNON\nCONNECTED\nGRAPH!\nGRAPH\nNOT\nLOADED\ntry again");
        error=1;
        igraph_destroy(&g);
        igraph_matrix_destroy(&adtemp);
        return;
        
        
    }
    
    igraph_destroy(&g);
    igraph_destroy(&graph);
    igraph_destroy(&sgraph);
    
    //generate definitive graph;
    
    igraph_matrix_update(&admatrix,&adtemp);
    igraph_matrix_destroy(&adtemp);
    
    //generate graph from admatrix
    igraph_weighted_adjacency(&graph, &admatrix,IGRAPH_ADJ_UPPER, "w",0);
    
    //symmetrize admatrix
    /*for (int i=0; i<clusdim*clusnumber; ++i) {
        for(int j=i; j<clusdim*clusnumber;++j)
        {
            MATRIX(admatrix,j,i)=MATRIX(admatrix,i,j);
        }
    }*/
    
    
    nodesnumber=igraph_matrix_nrow(&admatrix);
    
    
    logDebug("\n\nnodesnumber=%i\n\n", nodesnumber);
    
    //print_matrix(&admatrix, stdout);
    
    igraph_copy(&sgraph,&graph);
    igraph_simplify(&sgraph, 1,1,0);
    igraph_to_undirected(&sgraph, IGRAPH_TO_UNDIRECTED_COLLAPSE,/*edge_comb=*/ 0);
    
    

    
    //----------- prob matrix-----------------
    //ROW NORMALIZATION:
    igraph_vector_init(&row,nodesnumber);
    for (int i=0; i<nodesnumber; ++i) {
        igraph_matrix_get_row(&admatrix,&row,i);
        
        igraph_vector_scale(&row, 1./(igraph_vector_sum(&row)));
        igraph_matrix_set_row(&admatrix, &row,i);
        
    }
    igraph_vector_destroy(&row);
    
    
    //print_matrix_ur(&admatrix,stdout);
    
    //layout
    
    printf("-layout ...-");
    igraph_matrix_null(&layout);
    igraph_matrix_init(&layout, 0, 0);
    
    igraph_t tempgraph;
    igraph_vector_t weight;
    
    igraph_weighted_adjacency(&tempgraph, &admatrix,IGRAPH_ADJ_UNDIRECTED, "w",0);
    
    igraph_vector_init(&weight,0);
    
    EANV(&tempgraph,"w",&weight);
    
    
    
    igraph_destroy(&tempgraph);
    
    igraph_vector_scale(&weight,10.);
    
    igraph_layout_fruchterman_reingold(&graph, &layout,500, (float)nodesnumber,
                                       sqrt((double)nodesnumber), 1.5,
                                       (float)(nodesnumber*nodesnumber), 0,
                                       &weight,
                                       NULL,NULL,NULL,NULL);
    
    igraph_vector_destroy(&weight);
    
    adjustlayout();
    
    printf("-layout done-");
    
    fflush(stdout);
    
    
    error=0;
    
    haveloadedstate=0;
    
    
    
}





void generateclustsympots(int clusnumber, int clusdim, int cnodes, double interconn, double intraconn, double wall, int dissipate, double dissipation, int dnode){
    
    //i mainly work on the admatrix.
    
    init_genrand(0);
    
    igraph_matrix_t pots; //potentials of the clusters
    igraph_matrix_t wpots; // potentials of the wall
    
    igraph_matrix_t adtemp;
    igraph_t g;
    igraph_bool_t connected;
    
    int source;
    int totnodes;
    
    double actualinter=0, actualintra=0;
    
    
    igraph_vector_t row;
    
    
    //genreate potentials for CLUSTERS
    igraph_matrix_init(&pots,clusnumber,clusdim);
    igraph_matrix_null(&pots);
    
    for (int i=0; i<clusnumber; ++i) {
        for (int j=0; j<clusdim; ++j) {
            //random potentials:
            //MATRIX(pots,i,j)=i*wall+genrand_real1()*j;
            
            //fixed potentials
            MATRIX(pots,i,j)=0;
            
        }
    }
    
    init_genrand(time(NULL));
    
    
    //calculate number of bridges ("bridges")
    int nbridges=0;
    for (int i=2; i<=clusnumber; ++i) {
        nbridges=nbridges+(i-1);
    }
    
    printf("\n nbridges=%i", nbridges);
    
    
    igraph_matrix_init(&adtemp,(clusnumber*clusdim)+(cnodes*nbridges),(clusnumber*clusdim)+(cnodes*nbridges));
    igraph_matrix_null(&adtemp);
    
    //generate clusters
    for(int cluster=0; cluster<clusnumber; ++cluster){
        
        igraph_t gtemp;
        igraph_matrix_t adclus;
        
        igraph_erdos_renyi_game(&gtemp, IGRAPH_ERDOS_RENYI_GNP, clusdim, interconn, 0, 0);
        actualinter=actualinter+igraph_ecount(&gtemp);
        
        //ad clusters to adtemp
        igraph_matrix_init(&adclus,0,0);
        igraph_get_adjacency(&gtemp, &adclus, IGRAPH_GET_ADJACENCY_BOTH, 0);
        
        printf("\n cluster %i \n", cluster); //print_matrix_ur(&adclus,stdout); printf("\n\n");
        
        for (int j=0; j<clusdim; ++j) {
            for(int k=0; k<clusdim; ++k) {
                //the link is weighet wrt potential
                //one direction
                MATRIX(adtemp,(clusdim*cluster)+j,(clusdim*cluster)+k) =
                MATRIX(adclus,j,k)*exp((MATRIX(pots,cluster,j)-MATRIX(pots,cluster,k))/2);
                //other direction!
                
                MATRIX(adtemp,(clusdim*cluster)+k,(clusdim*cluster)+j) =
                MATRIX(adclus,k,j)*exp((MATRIX(pots,cluster,k)-MATRIX(pots,cluster,j))/2);
                
            }
        }
        igraph_matrix_destroy(&adclus);
        igraph_destroy(&gtemp);
        
    }
    
    actualinter=actualinter/((clusdim-1)*clusdim);
    
    
    
    igraph_matrix_init(&wpots,cnodes,nbridges); igraph_matrix_null(&wpots);
    
    //set potentials
    for (int i=0; i<cnodes; ++i) {
        for (int j=0; j<nbridges; ++j) {
            double ii=(double)i;
            double cn=(double)cnodes;
            //quadratic
            //MATRIX(wpots,i,j)= (-1)* wall * (( (ii-cn/2)/cn)*((ii-cn/2)/cn)  -1 );
            //constant
            MATRIX(wpots,i,j)=wall; //  /(clusdim-1+clusnumber-1);
        }
    }
    
    //link the nodes of a bridge (on a line)
    for (int i=0; i<(cnodes-1); ++i) {
        for (int j=0; j<nbridges; ++j) {
            
            printf("\nBRIDGE %i, linking inner nodes (from %i    to %i)\n", j,i,i+1);
            //in one direction
            MATRIX(adtemp,(clusdim*clusnumber)+(cnodes*j)+i, (clusdim*clusnumber)+(cnodes*j)+(i+1))=
            exp( - (MATRIX(wpots,i,j)-MATRIX(wpots,i+1,j))/2);
            //and in the other one!
            MATRIX(adtemp,(clusdim*clusnumber)+(cnodes*j)+(i+1), (clusdim*clusnumber)+(cnodes*j)+(i))=
            exp( - (MATRIX(wpots,i+1,j)-MATRIX(wpots,i,j))/2);
            
            
            
        }
    }
    
    //link first and last node to clusters.
    int contator=0;
    for (int clusfrom=0; clusfrom<clusnumber; ++clusfrom) {
        for (int clusto=0; clusto<clusfrom; ++clusto) {
            
            {
                printf("\n from=%i to=%i contator=%i  clusdim=%i  clusnumber=%i\n", clusfrom, clusto,contator,clusdim, clusnumber);
                //FROM: first node of the bridge
                //for all the possible nodes of CLUSFROM
                for (int i=0; i<clusdim; ++i) {
                    
                    double coin;
                    coin=genrand_real1();
                    if(coin<intraconn){ actualintra=actualintra+1;
                        //in one direction
                        MATRIX(adtemp,(clusdim*clusfrom)+i, (clusdim*clusnumber)+(contator*cnodes))=
                        exp( - (MATRIX(wpots,0,contator)-MATRIX(pots,clusfrom,i))/2);
                        
                        //and in the other one!
                        MATRIX(adtemp,(clusdim*clusnumber)+(contator*cnodes),(clusdim*clusfrom)+i)=
                        exp( - (MATRIX(pots,clusfrom,i)-MATRIX(wpots,0,contator))/2);
                        
                    }
                }
                
            }
            
            //TO: last node of the bridge
            //for all the possible nodes of CLUSTO
            for (int i=0; i<clusdim; ++i) {
                double coin;
                coin=genrand_real1();
                if(coin<intraconn){actualintra=actualintra+1;
                    //in one direction
                    MATRIX(adtemp,   (clusdim*clusto)+i  ,     (clusdim*clusnumber)+(contator*cnodes)+cnodes-1   )=
                    exp( - (MATRIX(wpots,cnodes-1,contator)-MATRIX(pots,clusto,i))/2);
                    //and in the other one!
                    MATRIX(adtemp,   (clusdim*clusnumber)+(contator*cnodes)+cnodes-1  ,    (clusdim*clusto)+i    )=
                    exp( - (MATRIX(pots,clusto,i)-MATRIX(wpots,cnodes-1,contator))/2);
                    
                    
                    
                }
                
            }
            
            ++contator;
        }
    }
    
    
    actualintra=actualintra/(clusdim*clusnumber);
    
    
    igraph_matrix_destroy(&pots);
    igraph_matrix_destroy(&wpots);
    printf("\n\n\n\n");
    //print_matrix_ur(&adtemp,stdout);
    printf("\n\n\n\n");
    
    
    totnodes=igraph_matrix_nrow(&adtemp);
    
    
    //----------- prob matrix-----------------
    //ROW NORMALIZATION:
    igraph_vector_init(&row,nodesnumber);
    for (int i=0; i<totnodes; ++i) {
        igraph_matrix_get_row(&adtemp,&row,i);
        
        igraph_vector_scale(&row, 1./(igraph_vector_sum(&row)));
        igraph_matrix_set_row(&adtemp, &row,i);
        
    }
    igraph_vector_destroy(&row);
    
    //include dissipation.
    isdissipating=0;
    sourcenode=0;
    
    
    
    if(dissipate==1 && dissipation>0){
        
        isdissipating=1;
        sourcenode=dnode;
        
        
        //add the dissipation node to admatrix:
        igraph_matrix_add_rows(&adtemp, 1);
        for(int i=0; i<igraph_matrix_ncol(&adtemp);++i)
        {MATRIX(adtemp,totnodes,i)=0;}
        
        igraph_matrix_add_cols(&adtemp, 1);
        for(int i=0; i<igraph_matrix_nrow(&adtemp);++i)
        {MATRIX(adtemp,i,totnodes)=0;}
        
        printf("\n --- ADDING DISSIPATION: adtemp, clean version (added new row and new column, all 0) \n");
        //print_matrix_ur(&adtemp,stdout);printf("\n\n\n\n");
        
        //connect the sink
        for(int i=0;i<totnodes-(cnodes*nbridges); ++i){
            
            MATRIX(adtemp,i,totnodes)=dissipation;
            
        }
        
        //connect the source (dissipation node---->cols of the source)
        //connect the sink (rows of the sink cluster--->dissipation node)
        if(dnode>totnodes){MATRIX(adtemp,totnodes,((totnodes/2)))=1;}
        else{MATRIX(adtemp,totnodes,dnode)=1;}
        
        
        printf("\n connected version \n");
        //print_matrix_ur(&adtemp,stdout);printf("\n\n\n\n");
        
        
    }
    
    
    //generate graph from admatrix
    
    igraph_weighted_adjacency(&g, &adtemp,IGRAPH_ADJ_UNDIRECTED, "w",0);
    
    
    //check if it's connected
    igraph_is_connected(&g, &connected,IGRAPH_STRONG);
    
    if(connected){
        logDebug("\nThe network is connected.\n");
        graphisloaded=1;
    }
    
    else {
        logDebug("\nERROR: THE GRAPH GENERATED IS NOT CONNECTED.\n");
        graphisloaded=0; rewrite=1;
        sprintf(errorstring,"ERROR!\nNON\nCONNECTED\nGRAPH!\nGRAPH\nNOT\nLOADED\ntry again");
        error=1;
        igraph_destroy(&g);
        igraph_matrix_destroy(&adtemp);
        
        return;
        
        
    }
    
    igraph_destroy(&g);
    igraph_destroy(&graph);
    igraph_destroy(&sgraph);
    
    //generate definitive graph;
    
    igraph_matrix_update(&admatrix,&adtemp);
    igraph_matrix_destroy(&adtemp);
    
    //generate graph from admatrix
    igraph_weighted_adjacency(&graph, &admatrix,IGRAPH_ADJ_DIRECTED, "w",0);
    
    printf("\n \n ---------------------- NON SYMMETRIZED ADMATRIX!! --------------------------------\n \n ");
    print_matrix(&admatrix, stdout);
    
    
    //symmetrize admatrix
    /*for (int i=0; i<clusdim*clusnumber; ++i) {
        for(int j=i; j<clusdim*clusnumber;++j)
        {
            MATRIX(admatrix,j,i)=MATRIX(admatrix,i,j);
        }
    }
    */
    
    nodesnumber=igraph_matrix_nrow(&admatrix);
    
    
    logDebug("\n\nnodesnumber=%i\n\n", nodesnumber);
    
    printf("\n \n ---------------------- SYMMETRIZED ADMATRIX!! --------------------------------\n \n ");
    print_matrix(&admatrix, stdout);
    
    igraph_copy(&sgraph,&graph);
    igraph_simplify(&sgraph, 1,1,0);
    igraph_to_undirected(&sgraph, IGRAPH_TO_UNDIRECTED_COLLAPSE,/*edge_comb=*/ 0);
    
    
    //GET BRIDGES LINKS
    
    igraph_matrix_init(&bridgeslinks,nbridges,cnodes-1);
    igraph_matrix_null(&bridgeslinks);
    printf("\n BRIDGES: nbrigdes=%i, cnodes=%i\n", nbridges,cnodes);
    for(int i=0; i<nbridges; ++i){
        for(int j=0; j<cnodes-1; ++j){
            int ttt;
            igraph_get_eid(&sgraph, &ttt,(clusnumber*clusdim)+(nbridges*cnodes*i)+j, (clusnumber*clusdim)+(nbridges*cnodes*i)+j+1,0,0);
            MATRIX(bridgeslinks,i,j)=ttt;
            
        }
        
    }
    
    
    
    
    
    
    
    
    printf("\n   LINKS OF THE BRIDGES:   \n"); //print_matrix_ur(&bridgeslinks,stdout); printf("\n \n");
    
    

    
    //----------- prob matrix-----------------
    //ROW NORMALIZATION:
    igraph_vector_init(&row,nodesnumber);
    for (int i=0; i<nodesnumber; ++i) {
        igraph_matrix_get_row(&admatrix,&row,i);
        
        igraph_vector_scale(&row, 1./(igraph_vector_sum(&row)));
        igraph_matrix_set_row(&admatrix, &row,i);
        
    }
    igraph_vector_destroy(&row);
    
    
    //print_matrix_ur(&admatrix,stdout);
    
    //layout
    
    printf("\n-layout ...-");
    igraph_matrix_null(&layout);
    igraph_matrix_init(&layout, 0, 0);
    
    igraph_t tempgraph;
    igraph_vector_t weight;
    
    igraph_weighted_adjacency(&tempgraph, &admatrix,IGRAPH_ADJ_DIRECTED, "w",0);
    
    igraph_vector_init(&weight,0);
    
    EANV(&tempgraph,"w",&weight);
    
    
    
    igraph_destroy(&tempgraph);
    
    igraph_vector_scale(&weight,10.);
    
    igraph_layout_fruchterman_reingold(&graph, &layout,500, (float)nodesnumber,
                                       sqrt((double)nodesnumber), 1.5,
                                       (float)(nodesnumber*nodesnumber), 0,
                                       &weight,
                                       NULL,NULL,NULL,NULL);
    
    igraph_vector_destroy(&weight);
    
    adjustlayout();
    
    printf("-layout done-\n \n");
    
    printf(" \n \n \n ---- actual INTER = %f    actual INTRA =  %f   ------- \n \n \n", actualinter, actualintra);
    
    
    fflush(stdout);
    
    
    error=0;
    
    haveloadedstate=0;
    
    
    
}




void generateclusger2(int clusnumber1, int clusnumber2, int clusdim1, int cnodes1,
                      int cnodes2, double interconn, double intraconn1, double intraconn2,
                      double wall1, double wall2, int dissipate, double dissipation, int dnode){
    
    //i mainly work on the admatrix.
    
    init_genrand(0);
    
    igraph_matrix_t pots; //potentials of the clusters
    igraph_matrix_t wpots; // potentials of the wall
    
    igraph_matrix_t adtemp;
    igraph_t g;
    igraph_bool_t connected;
    igraph_vector_t row;
    
    int sink;
    
    int tempnodesnumber; int totnodes;
    
    int contator;
    
    double actualinter=0, actualintra1=0, actualintra2=0;
    
    
    //genreate potentials for CLUSTERS
    igraph_matrix_init(&pots,clusnumber1,clusdim1);
    igraph_matrix_null(&pots);
    
    for (int i=0; i<clusnumber1; ++i) {
        for (int j=0; j<clusdim1; ++j) {
            //random potentials:
            //MATRIX(pots,i,j)=i*wall+genrand_real1()*j;
            
            //fixed potentials
            MATRIX(pots,i,j)=0;
            
        }
    }
    
    init_genrand(time(NULL));
    
    
    //calculate number of bridges
    int nbridges1=0;
    for (int i=2; i<=clusnumber1; ++i) {
        nbridges1=nbridges1+(i-1);
    }
    
    int nbridges2=0;
    for (int i=2; i<=clusnumber2; ++i) {
        nbridges2=nbridges2+(i-1);
    }
    
    printf("\n nbridges1=%i    nbridges2=%i  ", nbridges1, nbridges2);
    
    
    int bigclusdim;
    bigclusdim= (clusnumber1*clusdim1) + (cnodes1*nbridges1);
    
    int adtempdim;
    adtempdim= (clusnumber2*bigclusdim) +  (cnodes2*nbridges2);
    
    
    igraph_matrix_init(&adtemp,adtempdim,adtempdim);
    igraph_matrix_null(&adtemp);
    
    
    
    //--------------------------- GENERATE INTERN CLUSTERS. (LAYER 1=INTERN LAYER) --------
    // i generate 1 external cluster made up of intern clusters. (the prorotype)
    //then, the extern clusters will be a copy of the prototype.
    
    
    for(int cluster=0; cluster<clusnumber1; ++cluster){
        
        igraph_t gtemp;
        igraph_matrix_t adclus;
        
        igraph_erdos_renyi_game(&gtemp, IGRAPH_ERDOS_RENYI_GNP, clusdim1, interconn, 0, 0);
        actualinter=actualinter+igraph_ecount(&gtemp);
        
        //ad clusters to adtemp
        igraph_matrix_init(&adclus,0,0);
        igraph_get_adjacency(&gtemp, &adclus, IGRAPH_GET_ADJACENCY_BOTH, 0);
        
        printf("\n cluster %i \n", cluster); //print_matrix_ur(&adclus,stdout); printf("\n\n");
        
        for (int j=0; j<clusdim1; ++j) {
            for(int k=0; k<clusdim1; ++k) {
                //the link is weighet wrt potential
                //one direction
                MATRIX(adtemp,(clusdim1*cluster)+j,(clusdim1*cluster)+k) =
                MATRIX(adclus,j,k)*exp((MATRIX(pots,cluster,j)-MATRIX(pots,cluster,k))/2);
                //other direction!
                
                MATRIX(adtemp,(clusdim1*cluster)+k,(clusdim1*cluster)+j) =
                MATRIX(adclus,k,j)*exp((MATRIX(pots,cluster,k)-MATRIX(pots,cluster,j))/2);
                
            }
        }
        igraph_matrix_destroy(&adclus);
        igraph_destroy(&gtemp);
        
    }
    
    actualinter=actualinter/((clusdim1-1)*clusdim1);
    
    
    //generate WALLS (L1)
    igraph_matrix_init(&wpots,cnodes1,nbridges1); igraph_matrix_null(&wpots);
    
    //set potentials
    for (int i=0; i<cnodes1; ++i) {
        for (int j=0; j<nbridges1; ++j) {
            double ii=(double)i;
            double cn=(double)cnodes1;
            //quadratic
            //MATRIX(wpots,i,j)= (-1)* wall * (( (ii-cn/2)/cn)*((ii-cn/2)/cn)  -1 );
            //constant
            MATRIX(wpots,i,j)=wall1; // /(clusdim1-1+clusnumber1-1+clusnumber2-1);
        }
    }
    
    //link the nodes of a bridge (on a line)
    for (int i=0; i<(cnodes1-1); ++i) {
        for (int j=0; j<nbridges1; ++j) {
            
            printf("\nBRIDGE %i, linking inner nodes (from %i    to %i)\n", j,i,i+1);
            //in one direction
            MATRIX(adtemp,(clusdim1*clusnumber1)+(cnodes1*j)+i, (clusdim1*clusnumber1)+(cnodes1*j)+(i+1))=
            exp( - (MATRIX(wpots,i,j)-MATRIX(wpots,i+1,j))/2);
            //and in the other one!
            MATRIX(adtemp,(clusdim1*clusnumber1)+(cnodes1*j)+(i+1), (clusdim1*clusnumber1)+(cnodes1*j)+(i))=
            exp( - (MATRIX(wpots,i+1,j)-MATRIX(wpots,i,j))/2);
            
            
            
        }
    }
    
    //link first and last node to clusters.
    contator=0;
    for (int clusfrom=0; clusfrom<clusnumber1; ++clusfrom) {
        for (int clusto=0; clusto<clusfrom; ++clusto) {
            
            
            printf("\n from=%i to=%i contator=%i  clusdim=%i  clusnumber=%i\n", clusfrom, clusto,contator,clusdim1, clusnumber1);
            //FROM: first node of the bridge
            //for all the possible nodes of CLUSFROM
            for (int i=0; i<clusdim1; ++i) {
                
                double coin;
                coin=genrand_real1();
                if(coin<intraconn1){ actualintra1=actualintra1+1;
                    //in one direction
                    MATRIX(adtemp,(clusdim1*clusfrom)+i, (clusdim1*clusnumber1)+(contator*cnodes1))=
                    exp( - (MATRIX(wpots,0,contator)-MATRIX(pots,clusfrom,i))/2);
                    
                    //and in the other one!
                    MATRIX(adtemp,(clusdim1*clusnumber1)+(contator*cnodes1),(clusdim1*clusfrom)+i)=
                    exp( - (MATRIX(pots,clusfrom,i)-MATRIX(wpots,0,contator))/2);
                    
                }
            }
            
            
            
            //TO: last node of the bridge
            //for all the possible nodes of CLUSTO
            for (int i=0; i<clusdim1; ++i) {
                double coin;
                coin=genrand_real1();
                if(coin<intraconn1){actualintra1=actualintra1+1;
                    //in one direction
                    MATRIX(adtemp,   (clusdim1*clusto)+i  ,     (clusdim1*clusnumber1)+(contator*cnodes1)+cnodes1-1   )=
                    exp( - (MATRIX(wpots,cnodes1-1,contator)-MATRIX(pots,clusto,i))/2);
                    //and in the other one!
                    MATRIX(adtemp,   (clusdim1*clusnumber1)+(contator*cnodes1)+cnodes1-1  ,    (clusdim1*clusto)+i    )=
                    exp( - (MATRIX(pots,clusto,i)-MATRIX(wpots,cnodes1-1,contator))/2);
                    
                    
                    
                }
                
            }
            
            ++contator;
        }
    }
    
    
    actualintra1=actualintra1/(clusdim1*clusnumber1);
    
    igraph_matrix_destroy(&wpots);
    
    
    
    //--------------------------------------------- GENERATE EXTERN CLUSTERS
    
    for (int bigclus=1; bigclus<clusnumber2; ++bigclus) {
        for (int i=0; i<bigclusdim; ++i) {
            for (int j=0; j<bigclusdim; ++j) {
                
                MATRIX(adtemp,(bigclus*bigclusdim)+i,(bigclus*bigclusdim)+j)=MATRIX(adtemp,i,j);
                
                
            }
            
            
        }
    }
    
    
    
    printf("\n\n                NON CONNECTED BIG CLUSTERS              \n\n");
    //print_matrix_ur(&adtemp,stdout);
    printf("\n\n\n\n");
    
    
    
    
    //--------------------------------------------- GENERATE EXTERN BRIDGES
    
    
    
    
    
    //generate WALLS (L2)
    igraph_matrix_init(&wpots,cnodes2,nbridges2);
    igraph_matrix_null(&wpots);
    
    //set potentials
    for (int i=0; i<cnodes2; ++i) {
        for (int j=0; j<nbridges2; ++j) {
            double ii=(double)i;
            double cn=(double)cnodes2;
            //quadratic
            //MATRIX(wpots,i,j)= (-1)* wall * (( (ii-cn/2)/cn)*((ii-cn/2)/cn)  -1 );
            //constant
            MATRIX(wpots,i,j)=wall2; // /(((clusdim1-1)+clusnumber1-1+clusnumber2-1));
        }
    }
    
    //link the nodes of a bridge (on a line)
    for (int i=0; i<(cnodes2-1); ++i) {
        for (int j=0; j<nbridges2; ++j) {
            
            printf("\nBRIDGE %i, linking inner nodes (from %i    to %i)\n", j,i,i+1);
            //in one direction
            MATRIX(adtemp,(bigclusdim*clusnumber2)+(cnodes2*j)+i, (bigclusdim*clusnumber2)+(cnodes2*j)+(i+1))=
            exp( - (MATRIX(wpots,i,j)-MATRIX(wpots,i+1,j))/2);
            //and in the other one!
            MATRIX(adtemp,(bigclusdim*clusnumber2)+(cnodes2*j)+(i+1), (bigclusdim*clusnumber2)+(cnodes2*j)+(i))=
            exp( - (MATRIX(wpots,i+1,j)-MATRIX(wpots,i,j))/2);
            
            
            
        }
    }
    
    //link first and last node to clusters.
    contator=0;
    for (int clusfrom=0; clusfrom<clusnumber2; ++clusfrom) {
        for (int clusto=0; clusto<clusfrom; ++clusto) {
            
            printf("\n from=%i to=%i contator=%i  clusdim1=%i  clusnumber1=%i\n", clusfrom, clusto,contator,clusdim1, clusnumber1);
            
            //FROM: first node of the bridge
            //for all the possible nodes of CLUSFROM
            for (int i=0; i<(clusnumber1*clusdim1); ++i) {
                
                double coin;
                coin=genrand_real1();
                if(coin<intraconn2){
                    actualintra2=actualintra2+1;
                    
                    int cn, bn; //cluster node(to connect), bridge node (toconnect)
                    cn=(bigclusdim*clusfrom)+i;
                    bn=(bigclusdim*clusnumber2)+(contator*cnodes2);
                    
                    double cnp, bnp; //potentials of cluster node, bridge node.
                    cnp=MATRIX(pots,i/clusdim1,i%clusdim1);
                    bnp=MATRIX(wpots,0,contator);
                    
                    printf("\n [i=%i] cnp=%f  bnp=%f   \n\n pots:\n\n",i, cnp, bnp);
                    //print_matrix_ur(&pots,stdout); printf("\n  wpots: \n"); print_matrix_ur(&wpots,stdout);
                    
                    //in one direction
                    MATRIX(adtemp,  cn, bn)=    exp( -(bnp-cnp)/2);
                    
                    //and in the other one!
                    MATRIX(adtemp,  bn, cn)=    exp( -(cnp-bnp)/2);
                    
                }
            }
            
            
            
            //TO: last node of the bridge
            //for all the possible nodes of CLUSTO
            for (int i=0; i<(clusnumber1*clusdim1); ++i) {
                double coin;
                coin=genrand_real1();
                if(coin<intraconn2){
                    actualintra2=actualintra2+1;
                    
                    int cn, bn; //cluster node(to connect), bridge node (toconnect)
                    cn=(bigclusdim*clusto)+i ;
                    bn=(bigclusdim*clusnumber2)+(contator*cnodes2)+cnodes2-1 ;
                    
                    //in one direction
                    MATRIX(adtemp, cn,bn)=
                    exp( - (MATRIX(wpots,cnodes2-1,contator)-MATRIX(pots,i/clusdim1,i%clusdim1))/2);
                    
                    //and in the other one!
                    MATRIX(adtemp,bn,cn)=
                    exp( - (MATRIX(pots,i/clusdim1,i%clusdim1)-MATRIX(wpots,cnodes2-1,contator))/2);
                    
                    
                    
                }
                
            }
            
            ++contator;
        }
    }
    
    
    actualintra2=actualintra2/(bigclusdim*clusnumber2);
    
    
    igraph_matrix_destroy(&pots);
    igraph_matrix_destroy(&wpots);
    
    
    
    
    
    
    printf("\n\n  THIS IS ADTEMP:   \n\n");
    //print_matrix_ur(&adtemp,stdout);
    printf("\n\n\n\n");
    
    
    //check if i have some NAN.
    for (int i=0; i<igraph_matrix_ncol(&adtemp); ++i) {
        for (int j=0; j<igraph_matrix_ncol(&adtemp); ++j) {
            
            if(MATRIX(adtemp,i,j)!=MATRIX(adtemp,i,j)){
                
                logDebug("\nERROR: ADTEMP HAS 'NAN' VALUES\n");
                graphisloaded=0; rewrite=1;
                sprintf(errorstring,"ERROR!\nNAN\nVALUES\n!\nGRAPH\nNOT\nLOADED\ntry again");
                error=1;
                igraph_matrix_destroy(&adtemp);
                
                return;
                
                
            }
        }
    }
    
    
    
    
    
    tempnodesnumber=igraph_matrix_nrow(&adtemp);
    
    
    
    //(OLD)new version of dissipation
   /* isdissipating=0;
    if(dissipate==1 && dissipation>0){
        isdissipating=1;
        
        //select source *node*
        if(source>tempnodesnumber-1){
            printf("\n INVALID SOURCE. I PICKED ONE RANDOMLY. (sourse was:%i, but we have just %i nodes!)\n\n",
                    source,tempnodesnumber);
            source=genrand_int31()%tempnodesnumber;}
        
        printf("\n\n\n\n -------------------------- dissipate=%i, ------>  source=%i\n\n\n", dissipate, source);
        
        
        //add the dissipation node to admatrix:
        igraph_matrix_add_rows(&adtemp, 1);
        for(int i=0; i<igraph_matrix_ncol(&adtemp);++i){MATRIX(adtemp,igraph_matrix_nrow(&adtemp)-1,i)=0;}
        igraph_matrix_add_cols(&adtemp, 1);
        for(int i=0; i<igraph_matrix_nrow(&adtemp);++i){MATRIX(adtemp,i,igraph_matrix_ncol(&adtemp)-1)=0;}
        
        printf("\n clean version \n");
        print_matrix_ur(&adtemp,stdout);printf("\n\n\n\n");
        
        //connect the sink (rows of the sink cluster--->dissipation node)
        for(sink=0;sink<clusnumber2;++sink){
            for(int i=0;i<(bigclusdim-(cnodes1*nbridges1));++i){
            
            MATRIX(adtemp,(sink*bigclusdim)+i,igraph_matrix_ncol(&adtemp)-1)=dissipation;
            
        }}
        
        //connect the source (dissipation node---->cols of the source)
        //connect the sink (rows of the sink cluster--->dissipation node)
        
        MATRIX(adtemp,igraph_matrix_nrow(&adtemp)-1,source)=1;
         
            
        
        printf("\n connected version \n");
        print_matrix_ur(&adtemp,stdout);printf("\n\n\n\n");
        
    }
    
    */
    
    
    
    
    totnodes=igraph_matrix_nrow(&adtemp);
    
    
    //----------- prob matrix-----------------
    //ROW NORMALIZATION:
    igraph_vector_init(&row,nodesnumber);
    for (int i=0; i<totnodes; ++i) {
        igraph_matrix_get_row(&adtemp,&row,i);
        
        igraph_vector_scale(&row, 1./(igraph_vector_sum(&row)));
        igraph_matrix_set_row(&adtemp, &row,i);
        
    }
    igraph_vector_destroy(&row);
    
    //include dissipation.
    isdissipating=0;
    sourcenode=0;
    
    
    if(dissipate==1 && dissipation>0){
        
        isdissipating=1;
        sourcenode=dnode;
        
        
        //add the dissipation node to admatrix:
        igraph_matrix_add_rows(&adtemp, 1);
        for(int i=0; i<igraph_matrix_ncol(&adtemp);++i)
        {MATRIX(adtemp,totnodes,i)=0;}
        
        igraph_matrix_add_cols(&adtemp, 1);
        for(int i=0; i<igraph_matrix_nrow(&adtemp);++i)
        {MATRIX(adtemp,i,totnodes)=0;}
        
        printf("\n --- ADDING DISSIPATION: adtemp, clean version (added new row and new column, all 0) \n");
        //print_matrix_ur(&adtemp,stdout);printf("\n\n\n\n");
        
        //connect the sink
        for(int k=0;k<clusnumber2;++k){
            for(int i=0;i<bigclusdim-(cnodes1*nbridges1); ++i){
    
            MATRIX(adtemp,(bigclusdim*k)+i,totnodes)=dissipation;
            
        }
        }
        
        //connect the source (dissipation node---->cols of the source)
        //connect the sink (rows of the sink cluster--->dissipation node)
        if(dnode>totnodes){MATRIX(adtemp,totnodes,((totnodes/2)))=1;}
        else{MATRIX(adtemp,totnodes,dnode)=1;}
        
        
        printf("\n connected version \n");
        //print_matrix_ur(&adtemp,stdout);printf("\n\n\n\n");
        
        
    }
    
    
    
    
    //generate graph from admatrix
    
    igraph_weighted_adjacency(&g, &adtemp,IGRAPH_ADJ_UNDIRECTED, "w",0);
    
    
    //check if it's connected
    igraph_is_connected(&g, &connected,IGRAPH_STRONG);
    
    if(connected){
        logDebug("\nThe network is connected.\n");
        graphisloaded=1;
    }
    
    else {
        logDebug("\nERROR: THE GRAPH GENERATED IS NOT CONNECTED.\n");
        graphisloaded=0; rewrite=1;
        sprintf(errorstring,"ERROR!\nNON\nCONNECTED\nGRAPH!\nGRAPH\nNOT\nLOADED\ntry again");
        error=1;
        igraph_destroy(&g);
        igraph_matrix_destroy(&adtemp);
        
        return;
        
        
    }
    
    igraph_destroy(&g);
    igraph_destroy(&graph);
    igraph_destroy(&sgraph);
    
    //generate definitive graph;
    logDebug("    GENERATE DEFINITIVE GRAPH \n\n");
    igraph_matrix_update(&admatrix,&adtemp);
    igraph_matrix_destroy(&adtemp);
    
    //generate graph from admatrix
    
    logDebug("    FROM ADMATRIX \n\n");
    igraph_weighted_adjacency(&graph, &admatrix,IGRAPH_ADJ_DIRECTED, "w",0);
    
    
     logDebug("    SYMMETRIXE MATRIX \n\n");
     //symmetrize admatrix
    /* for (int i=0; i<igraph_matrix_nrow(&admatrix); ++i) {
     for(int j=i; j<igraph_matrix_nrow(&admatrix);++j)
     {
     MATRIX(admatrix,j,i)=MATRIX(admatrix,i,j);
     }
     }*/
    
    
    nodesnumber=igraph_matrix_nrow(&admatrix);
    
    
    logDebug("\n\nnodesnumber=%i\n\n", nodesnumber);
    
    //print_matrix(&admatrix, stdout);
    
    igraph_copy(&sgraph,&graph);
    igraph_simplify(&sgraph, 1,1,0);
    igraph_to_undirected(&sgraph, IGRAPH_TO_UNDIRECTED_COLLAPSE,/*edge_comb=*/ 0);
    
    logDebug("    SIMPLIFIED AND COLLAPSED.   NOW GETTING BRIDGES LINKS \n\n");
    
    //----------------GET BRIDGES LINKS
    
    int maxbridgesnodes=0;
    maxbridgesnodes=cnodes1; if(cnodes2>cnodes1) maxbridgesnodes=cnodes2;
    
    logDebug("  GENERATING BRIDGESLINKS MATRIX \n\n");
    igraph_matrix_init(&bridgeslinks,(nbridges1*clusnumber2)+nbridges2,maxbridgesnodes-1);
    igraph_matrix_null(&bridgeslinks);
    printf("\n BRIDGES : nbrigdes=%i, cnodes=%i\n", (nbridges1*clusnumber2)+nbridges2,maxbridgesnodes);
    // layer 1
    logDebug("  BRIDGESLINKS MATRIX   --- GETTING LAYER 1\n\n");
    for (int cluin=0; cluin<clusnumber2; ++cluin) {
        
        for(int i=0; i<nbridges1; ++i){
            for(int j=0; j<cnodes1-1; ++j){
                int ttt;
                igraph_get_eid(&sgraph, &ttt,
                               (cluin*bigclusdim)+((clusnumber1*clusdim1)+(nbridges1*cnodes1*i)+j),
                               (cluin*bigclusdim)+ (clusnumber1*clusdim1)+(nbridges1*cnodes1*i)+j+1,
                               0,0);
                MATRIX(bridgeslinks,(cluin*nbridges1)+i,j)=ttt;
                
            }
            
        }
    }
    // layer 2
    logDebug("  BRIDGESLINKS MATRIX   --- GETTING LAYER 2\n\n");
    
    for(int i=0; i<nbridges2; ++i){
        for(int j=0; j<cnodes2-1; ++j){
            int ttt;
            igraph_get_eid(&sgraph, &ttt,(clusnumber2*bigclusdim)+(nbridges2*cnodes2*i)+j, (clusnumber2*bigclusdim)+(nbridges2*cnodes2*i)+j+1,0,0);
            MATRIX(bridgeslinks,i,j)=ttt;
            
        }
        
    }
    
    
    
    
    
    
    
    printf("\n   LINKS OF THE BRIDGES:   \n"); //print_matrix_ur(&bridgeslinks,stdout); printf("\n \n");
    
    

    
    //----------- prob matrix-----------------
    //ROW NORMALIZATION:
    igraph_vector_init(&row,nodesnumber);
    for (int i=0; i<nodesnumber; ++i) {
        igraph_matrix_get_row(&admatrix,&row,i);
        igraph_vector_scale(&row, 1./(igraph_vector_sum(&row)));
        igraph_matrix_set_row(&admatrix, &row,i);
        
    }
    igraph_vector_destroy(&row);
    
    
    //print_matrix_ur(&admatrix,stdout);
    
    
    
    //check if i have some NAN (AGAIN).
    for (int i=0; i<igraph_matrix_ncol(&admatrix); ++i) {
        for (int j=0; j<igraph_matrix_ncol(&admatrix); ++j) {
            
            if(MATRIX(admatrix,i,j)!=MATRIX(admatrix,i,j)){
                
                logDebug("\nERROR: ADTEMP HAS 'NAN' VALUES\n");
                graphisloaded=0; rewrite=1;
                sprintf(errorstring,"ERROR!\nNAN\nVALUES\n!\nGRAPH\nNOT\nLOADED\ntry again");
                error=1;
                
                return;
                
                
            }
        }
    }
    
    
    
    
    
    //layout
    
    printf("\n--- DOING layout ...-\n\n"); fflush(stdout);
    
    igraph_matrix_null(&layout);
    igraph_matrix_init(&layout, 0, 0);
    
    igraph_t tempgraph;
    igraph_vector_t weight;
    
    igraph_weighted_adjacency(&tempgraph, &admatrix,IGRAPH_ADJ_DIRECTED, "w",0);
    
    
    igraph_vector_init(&weight,0);
    
    EANV(&tempgraph,"w",&weight);
    
    
    igraph_vector_scale(&weight,10.);
    
    
    igraph_layout_fruchterman_reingold(&tempgraph, &layout,500, (float)nodesnumber,
                                       sqrt((double)nodesnumber), 1.5,
                                       (float)(nodesnumber*nodesnumber), 0,
                                       &weight,
                                       NULL,NULL,NULL,NULL);
    
    igraph_destroy(&tempgraph);
    
    
    
    
    igraph_vector_destroy(&weight);
    fflush(stdout);
    
    fflush(stdout);
    

    
    adjustlayout();
    
        if(dissipate==1){MATRIX(layout,nodesnumber-1,0)=0;MATRIX(layout,nodesnumber-1,1)=0;}
    
    printf("-layout done-\n \n");
    
    printf(" \n \n \n ---- actual INTER = %f    actual INTRA1 =  %f  actual INTRA2 =  %f   ------- \n \n \n", actualinter, actualintra1,actualintra2);
    
    
    fflush(stdout);
    
    
    error=0;
    
    haveloadedstate=0;
    
    
}










void StationaryState() {
    
    
    //--------------------------------------------------- EIGENVALUES & EIGENVECTORS WITH LAPACK (included in igraph library)
    
    
    
    igraph_vector_t valuesreal, valuesimag;
    igraph_matrix_t vectorsleft;
    int info;
    double tempreal;
    int evec;
    int i;
    
    
    igraph_vector_init(&valuesreal,0);
    igraph_vector_init(&valuesimag,0);
    igraph_matrix_init(&vectorsleft, 0, 0);
    
    igraph_lapack_dgeev(&admatrix, &valuesreal, &valuesimag,&vectorsleft,NULL,&info);
    logDebug("\nLAPACK DGEEV INFO = %i \n",info);
    
    
    //print_matrix(&admatrix,stdout);
    
    printf("\n EIGENVALUES \n");
    for (i=0; i<nodesnumber; ++i) {printf("%i ::::: %f + i*%f \n", i, VECTOR(valuesreal)[i],VECTOR(valuesimag)[i]);}
    /*
     printf("\n EIGENVECTORS LEFT \n");
     print_matrix(&vectorsleft,stdout);
     */
    
    //get steady state
    evec=-1;
    tempreal=0;
    for (i=0; i<nodesnumber; ++i) {if(VECTOR(valuesreal)[i]>tempreal){ evec=i; tempreal=VECTOR(valuesreal)[i];}}
    
    logDebug("\nEigenvector ''EVEC=%i''  has eigenvalue ''EVAL=%f''.\n ",evec,tempreal);
    
    igraph_vector_init(&statstate,0);
    igraph_matrix_get_col(&vectorsleft,&statstate,evec);
    
    
    print_vector_line(&statstate,stdout);
    //PROBABILITY NORMALIZATION
    igraph_vector_scale(&statstate, 1./(igraph_vector_sum(&statstate)));
    logDebug("\nStationary State normalized.\n");
    print_vector_line(&statstate,stdout);
    
    printf("\n \n \n ");
    //get all the other states
    igraph_matrix_init(&estates,nodesnumber,nodesnumber);
    
    
    for (int j=0; j<nodesnumber; ++j) {
        
        igraph_vector_t tempvect;
        igraph_vector_init(&tempvect,0);
        evec=-1;
        tempreal=-1;
        for (i=0; i<nodesnumber; ++i) {if(VECTOR(valuesreal)[i]>tempreal){ evec=i; tempreal=VECTOR(valuesreal)[i];}}
        
        igraph_matrix_get_col(&vectorsleft,&tempvect,evec);
        
        logDebug("\nEigenvector ''EVEC=%i''  has eigenvalue ''EVAL=%f''    -- NORM %f.\n ",evec,tempreal, (igraph_vector_sum(&tempvect)));
        // print_vector_line(&tempvect,stdout);
        
        VECTOR(valuesreal)[evec]=-1000;
        
        /*
         //PROBABILITY NORMALIZATION
         igraph_vector_scale(&tempvect, 1./(igraph_vector_sum(&tempvect)));
         logDebug("\n Evec normalized.\n");
         print_vector_line(&tempvect,stdout);
         */
        igraph_matrix_set_col(&estates,&tempvect, j);
        
        igraph_vector_destroy(&tempvect);
        
    }
    
    
    printf("\n \n eigenvectors:");
   // print_matrix_ur(&estates,stdout);
    printf("\n\n");
    
    
    
    //_______________________________________ DESTROY DESTROY DESTROY DESTROY DESTROY DESTROY DESTROY DESTROY
    igraph_vector_destroy(&valuesimag);
    igraph_vector_destroy(&valuesreal);
    igraph_matrix_destroy(&vectorsleft);
    
    
}










//--------------------------------------------
void loadnetwork(char *mypath) {
    
    int mynodesnumber;
    int mytotrun;
    int myusesteady=0;
    int havestate=0;
    
    double mytempd;
    int mytempi;
    
    igraph_t graphtemp;
    igraph_matrix_t adtemp;
    igraph_vector_t steadytemp;
    igraph_matrix_t statetemp;
    
    igraph_vector_t row;
    igraph_bool_t connected;
    
    FILE *instream;
    instream=fopen(mypath,"r");
    
    char s[100];
    
    
    //load graph:
    
    //read head
    fgets(s, 100, instream);
    if(s[0]!='#'||s[1]!='l'||s[2]!='f'||s[3]!='\n'){
        logDebug("\nERROR: NOT A LAURAFILE\n");
        error=1; rewrite=1;
        sprintf(errorstring,"ERROR!\nNOT\nLAURA\nFILE!");
        return;}
    
    fscanf(instream,"%i",&mytotrun); logDebug("\nmytotrun = %i\n",mytotrun);
    fscanf(instream,"%i",&mynodesnumber); logDebug("\nmynodesnumber = %i\n",mynodesnumber);
    
    fgets(s, 100, instream);
    if(s[0]!='\n'){
        logDebug("\nERROR: CORRUPTED LAURAFILE (1)\n");
        error=1;  rewrite=1;
        sprintf(errorstring,"ERROR!\nCORRU\nPTED\nLAURA\nFILE!");
        return;}
    
    
    
    igraph_matrix_init(&adtemp,mynodesnumber,mynodesnumber);
    igraph_vector_init(&steadytemp,mynodesnumber);
    igraph_matrix_init(&statetemp,mynodesnumber,mytotrun);
    igraph_matrix_init(&loadedstate,mynodesnumber,mytotrun);
    
    
    //read PART I
    fgets(s, 100, instream);
    if(s[0]!='#'){
        logDebug("\nERROR: CORRUPTED LAURAFILE (2)\n");
        error=1;  rewrite=1;
        sprintf(errorstring,"ERROR!\nCORRU\nPTED\nLAURA\nFILE!");
        return;}
    
    
    for(int i=0;i<mynodesnumber;++i){
        for(int j=0;j<mynodesnumber;++j){
            fscanf(instream,"%lf",&mytempd);
            MATRIX(adtemp,i,j)=mytempd;}
    }
    
    fgets(s, 100, instream);
    
    //read PART II
    fgets(s, 100, instream);
    if((s[0]!='#'|| s[1]!='#')&&(s[0]!='#'||s[1]!='e'||s[2]!='l'||s[3]!='f')){
        logDebug("\nERROR: CORRUPTED LAURAFILE (3)\n");
        error=1; rewrite=1;
        sprintf(errorstring,"ERROR!\nCORRU\nPTED\nLAURA\nFILE!");
        return;}
    
    if (s[2]=='\n') {
        
        for(int i=0;i<mynodesnumber;++i){
            fscanf(instream,"%lf",&mytempd);
            VECTOR(steadytemp)[i]=mytempd;
        }
        
        myusesteady=1;
        
        fgets(s, 100, instream);
        fgets(s, 100, instream);
        
    }
    
    
    //read PART III
    if ((s[0]!='#'||s[1]!='#'||s[2]!='#')&&(s[0]!='#'||s[1]!='e'||s[2]!='l'||s[3]!='f')){
        logDebug("\nERROR: CORRUPTED LAURAFILE (4)\n");
        error=1;  rewrite=1;
        sprintf(errorstring,"ERROR!\nCORRU\nPTED\nLAURA\nFILE!");
        return;}
    
    if (s[3]=='\n'){
        
        for(int i=0;i<mynodesnumber;++i){
            for(int j=0;j<mytotrun;++j){
                fscanf(instream,"%lf",&mytempd);
                MATRIX(statetemp,i,j)=(int)mytempd;
                MATRIX(loadedstate,i,j)=(int)mytempd;
            }
        }
        
        havestate=1;
        
        fgets(s, 100, instream);
        fgets(s, 100, instream);
        
    }
    
    //read END
    if (s[0]!='#'||s[1]!='e'||s[2]!='l'||s[3]!='f'){
        logDebug("\nERROR: CORRUPTED LAURAFILE (5)\n");
        error=1;  rewrite=1;
        sprintf(errorstring,"ERROR!\nCORRU\nPTED\nLAURA\nFILE!");
        return;}
    
    logDebug("\nLaurafile correctly read.\n");
    
    
    fclose(instream);
    
    //logDebug("\nADJ MATRIX: ");
    //print_matrix(&adtemp,stdout);
    logDebug("\nState:\n");
    //print_matrix_ur(&statetemp,stdout);
    
    if(myusesteady==1){    logDebug("\nSteady State: ");
        print_vector_line(&steadytemp,stdout);}
    else logDebug("\nNo steady state. ");
    
    
    
    
    //SET ADMATRIX
    igraph_weighted_adjacency(&graphtemp, &adtemp,IGRAPH_ADJ_DIRECTED, 0,0);
    //check if it's connected
    igraph_is_connected(&graphtemp, &connected,IGRAPH_STRONG);
    
    if(connected){
        logDebug("\nThe network is connected.\n");
        graphisloaded=1;
    }
    
    else {
        logDebug("\nERROR: THE GRAPH LOADED IS NOT CONNECTED.\n");
        graphisloaded=0; rewrite=1;
        sprintf(errorstring,"ERROR!\nNON\nCONNECTED\nGRAPH!\nGRAPH\nNOT\nLOADED\ntry again");
        error=1;
        return;
    }
    
    
    //----------- prob matrix-----------------
    //ROW NORMALIZATION:
    igraph_vector_init(&row,mynodesnumber);
    for (int i=0; i<mynodesnumber; ++i) {
        igraph_matrix_get_row(&adtemp,&row,i);
        igraph_vector_scale(&row, 1./(igraph_vector_sum(&row)));
        igraph_matrix_set_row(&adtemp, &row,i);
    }
    igraph_vector_destroy(&row);
    
    logDebug("\nAdtemp is now a stochastic matrix.\n");
    fflush(stdout);
    
    
    //generate Layout
    igraph_matrix_init(&layout, 0, 0);
    
    
    igraph_layout_star(&graphtemp, &layout,0,0);
    
    
    adjustlayout();
    
    
    
    error=0;
    
    logDebug("\nLayout done.\n");
    
    fflush(stdout);
    
    //update everything
    igraph_matrix_update(&admatrix,&adtemp);
    igraph_copy(&graph, &graphtemp);
    igraph_copy(&sgraph,&graph);
    igraph_simplify(&sgraph, 1,1,0);
    igraph_to_undirected(&sgraph, IGRAPH_TO_UNDIRECTED_COLLAPSE,/*edge_comb=*/ 0);
    
    
    if(havestate=1){
        igraph_matrix_update(&state,&statetemp);
    }
    
    if(myusesteady==1){
        igraph_vector_init(&statstate,0);
        igraph_vector_update(&statstate,&steadytemp);}
    
    
    usesteady=myusesteady;
    haveloadedstate=havestate;
    nodesnumber=mynodesnumber;
    totrun=mytotrun;
    graphisloaded=1;
    error=0;
    
    
    igraph_matrix_destroy(&adtemp);
    igraph_matrix_destroy(&statetemp);
    igraph_vector_destroy(&steadytemp);
    
    igraph_destroy(&graphtemp);
}






//--------------------------------------------

void savenetwork(char *mypath) {
    
    
    FILE *outstream;
    outstream=fopen(mypath,"w");
    
    logDebug("\n SAVING... \n");
    
    //open file where i'll save
    if( (outstream=fopen(mypath,"w")) ==NULL) {
        logDebug("\nERROR SAVING\n");
        sprintf(errorstring,"ERROR\nSAVING\nCANT\nOPEN\n%s",mypath);
        error=1;  rewrite=1;
        return;
    }
    
    
    //----------------- WRITE A "LAURA" FILE::::::::::::::::
    
    //head
    fprintf(outstream,"#lf\n");
    fprintf(outstream,"%i\n",totrun);
    
    //PART I: NETWORK
    fprintf(outstream,"%i\n",nodesnumber); //number of nodes
    fprintf(outstream,"#\n");
    //print ADJACENCE MATRIX
    print_matrix(&admatrix,outstream);
    
    //PART II: STEADY STATE
    if(usesteady==1){
        fprintf(outstream,"##\n");
        print_vector_line(&statstate,outstream);
    }
    
    //PART III: CURRENT STATE
    fprintf(outstream,"###\n");
    print_matrix(&state, outstream);
    
    
    fprintf(outstream,"#elf\n");
    
    fclose(outstream);
    
    
    havepath=1;
    
    error=0;
    
    
    
}















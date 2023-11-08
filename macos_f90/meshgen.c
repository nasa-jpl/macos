/**
 *  meshgen.c
 */

#include "qhull_a.h"

int get_numfacets(facetT *facetlist, setT *facets,
                   boolT printall) 
{
  int numfacets, numsimplicial, numridges, totneighbors, 
      numcoplanars, numtricoplanars;
  int i, num;
  facetT *facet, **facetp;

  qh printoutnum= 0;
  qh_countfacets (facetlist, facets, printall, &numfacets, &numsimplicial,
      &totneighbors, &numridges, &numcoplanars, &numtricoplanars);
  return numfacets;
} /* get_numfacets */

#if 0
void printfacet3vertex(FILE *fp, facetT *facet, int format)
{
  vertexT *vertex, **vertexp;
  setT *vertices;

  /*printf("In io.c::qh_printfacet3vertex() ...\n");*/

  vertices= qh_facet3vertex (facet);

  FOREACHvertex_(vertices) {
    fprintf (fp, "%i ", qh_pointid(vertex->point));
  }
  fprintf (fp, "\n");
  qh_settempfree(&vertices);
} /* printfacet3vertex */
#endif


boolT prt_strt;

void printafacet(FILE *fp, int format, facetT *facet, int *fid_p,
                 boolT printall) 
{
  realT color[4], offset, dist, outerplane, innerplane;
  boolT zerodiv;
  coordT *point, *normp, *coordp, **pointp, *feasiblep;
  int k;
  vertexT *vertex, **vertexp;
  facetT *neighbor, **neighborp;
  setT *vertices;
  static int *fid_p_loc;

  if (!printall && qh_skipfacet (facet))
    return;
  if (facet->visible && qh NEWfacets && format != qh_PRINTfacets)
    return;
  qh printoutnum++;

  if (prt_strt) { 
    fid_p_loc=fid_p;
    prt_strt=False;
  }

  vertices= qh_facet3vertex (facet);

  FOREACHvertex_(vertices) {
    *fid_p_loc=qh_pointid(vertex->point); 
    /*fprintf (fp, "%i ", *fid_p_loc);*/ 
    ++fid_p_loc;
  }
  /*fprintf (fp, "\n");*/
  qh_settempfree(&vertices);
} /* printafacet */


int gen_trimesh(int, coordT*, int*, int**);

#define DIM 2     /* dimension of points, must be < 31 for SIZEcube */
/*#define SIZEcube (1<<DIM)*/
/*#define SIZEdiamond (2*DIM)*/
/*#define TOTpoints (SIZEcube + SIZEdiamond)*/
/*#define TOTpoints 7*/


/* *** coordT is now set to "double" *** */

/*
 * main interface routine (to F90 MACOS)
 */

void meshgen_(int *n_mesh_pts, double *mesh_coords,
              int *n_mesh_elts, int *mesh_elt_nodes_out)
{
  int *mesh_elt_nodes;

  printf("** meshgen_: number of mesh points = %i\n", *n_mesh_pts);

  /* tri-mesh generation */
  gen_trimesh(*n_mesh_pts, mesh_coords, n_mesh_elts, &mesh_elt_nodes);

  if (*n_mesh_elts > 3*(*n_mesh_pts)) {
    printf("** meshgen_: number mesh elements exceeds mesh_elt_nodes buffer\n");
    exit(0);
   }
  else {
    memcpy(mesh_elt_nodes_out,mesh_elt_nodes,
           3*(*n_mesh_elts)*sizeof(int));
  }
 
  /* Free buffer allocated in gen_trimesh() */ 
  free(mesh_elt_nodes);
} // meshgen_


#if 0
void meshgen_()
{
  FILE *fp;
  int n_mesh_pts, n_mesh_elts, ir, ie;
  coordT *mesh_points, d; 
  int *mesh_elt_nodes;

  /* Read in mesh point coords */
  fp = fopen("SurfMap011.txt","r");
  if (!fp) {
    printf("open file %s failed!\n", "SurfMap011.txt");
    exit(0);
  }

  fscanf(fp,"%i",&n_mesh_pts);
  printf("** number of mesh points = %i\n", n_mesh_pts);

  mesh_points = (coordT*) calloc(2*n_mesh_pts,sizeof(coordT));

  /* Read in mesh point coords  */
  for (ir=0; ir<n_mesh_pts; ++ir) {
    fscanf(fp,"%lf%lf%lf%lf",
      &mesh_points[2*ir],&mesh_points[2*ir+1],&d,&d);
    /*printf("%e %e\n",mesh_points[2*ir],mesh_points[2*ir+1]);*/ 
  }
  fclose(fp);

  /* tri-mesh generation */
  gen_trimesh(n_mesh_pts, mesh_points, &n_mesh_elts, &mesh_elt_nodes);

  /* Write out mesh triangle node ids */
  fp = fopen("SurfMap011_trimesh.txt","w");
  if (!fp) {
    printf("open file %s failed!\n", "SurfMap011_trimesh.txt");
    exit(0);
  }
  fprintf(fp,"%i\n", n_mesh_elts);
  for (ie=0; ie<n_mesh_elts; ++ie) {
    fprintf(fp,"%i %i %i\n", mesh_elt_nodes[3*ie],
            mesh_elt_nodes[3*ie+1],mesh_elt_nodes[3*ie+2]);
  }

  free(mesh_points);
  free(mesh_elt_nodes);
} /* meshgen_ */
#endif


int gen_trimesh(int TOTpoints, coordT *points_in,
                int *TOTelts, int **elts) 
{
  int dim= DIM;	          /* dimension of points */
  int numpoints;          /* number of points */
  /*coordT points[(DIM+1)*TOTpoints]; */ /* array of coordinates for each point */
  coordT *points;
  /*coordT *rows[TOTpoints];*/
  coordT **rows;
  boolT ismalloc= False;    /* True if qhull should free points in qh_freeqhull() 
			       or reallocation */
  char flags[250];          /* option flags for qhull, see qh_opt.htm */
  FILE *outfile= stdout;    /* output from qh_produce_output()
			       use NULL to skip qh_produce_output() */
  FILE *errfile= stderr;    /* error messages from qhull code */
  int exitcode;             /* 0 if no error from qhull */
  facetT *facet;	    /* set by FORALLfacets */
  int curlong, totlong;	    /* memory remaining after qh_memfreeshort */
  boolT printall;
  int i, ifacet, num_facets;
  int *facet_id_arr;

  /* Local storage allocation and initialization */
  points = (coordT*) calloc((DIM+1)*TOTpoints,sizeof(coordT));
  rows = (coordT**) calloc(TOTpoints,sizeof(coordT*)); 
  for (i=0; i<2*TOTpoints; ++i) points[i]=points_in[i];
     

  /*
    Run 2: Delaunay triangulation
  */

  printall=True;

  printf( "\ncompute %d-d Delaunay triangulation\n", dim);
  /*sprintf (flags, "qhull s i d Qu");*/
  sprintf (flags, "qhull s i d");

  numpoints= TOTpoints;
  printf("numpoints= %i \n",numpoints);

  /* mesh points now defined in main() */
  /*makeDelaunay(points, numpoints, dim, time(NULL));*/

  for (i=numpoints; i--; )
    rows[i]= points+dim*i;

  /* print out input, now turned off */
  /* qh_printmatrix (outfile, "input", rows, numpoints, dim); */

  /* Build new qhull data structure, defined in user.c */
  exitcode= qh_new_qhull(dim, numpoints, points, ismalloc,
                         flags, outfile, errfile);


  if (!exitcode) {                  /* if no error */
    /* 'qh facet_list' contains the convex hull */
    /* If you want a Voronoi diagram ('v') and do not request output (i.e., outfile=NULL), 
       call qh_setvoronoi_all() after qh_new_qhull(). */

    
    num_facets = get_numfacets(qh facet_list, NULL, !qh_ALL);
    /*printf("\nNumber of triangles in mesh = %d\n", num_facets);*/
    *TOTelts = num_facets;

    facet_id_arr = (int*) calloc(num_facets*3, sizeof(int));
    
    /*printf("Begin ForAll printafacet ...\n");*/
    prt_strt=True;
    FORALLfacet_(qh facet_list) 
      printafacet(stdout, qh_PRINTincidences, facet, facet_id_arr,
                  !qh_ALL);
   
    /* 
    printf("facet_id_arr = ");
    for (ifacet=0; ifacet<num_facets*3; ifacet++)
      printf("%i ", facet_id_arr[ifacet]);
    printf("\n");*/

    /*free(facet_id_arr);*/

    /*printf("End ForAll printafacet ...\n");*/
  }  /* if !exitcode */

  free(points);
  free(rows);

  qh_freeqhull(!qh_ALL);   /* free long memory */
  qh_memfreeshort (&curlong, &totlong);  /* free short memory and memory allocator */
  if (curlong || totlong) 
    fprintf (errfile, "qhull internal warning (user_eg, #2): \
     did not free %d bytes of long memory (%d pieces)\n", totlong, curlong);

  *elts = facet_id_arr;
  return exitcode;
} /* gen_trimesh */


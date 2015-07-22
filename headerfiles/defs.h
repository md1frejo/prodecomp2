  gsl_matrix *defs1=gsl_matrix_calloc(4,2);
  gsl_matrix *defs2=gsl_matrix_calloc(10,3);
  gsl_matrix *defs3=gsl_matrix_calloc(13,3);

  gsl_matrix_set(defs1,0,0,1);
  gsl_matrix_set(defs1,1,1,1);
  gsl_matrix_set(defs1,2,0,1);
  gsl_matrix_set(defs1,2,1,1);
  gsl_matrix_set(defs1,3,0,1);
  gsl_matrix_set(defs1,3,1,-1);  

  /*****************************/

  gsl_matrix_set(defs2,0,0,1);

  gsl_matrix_set(defs2,1,1,1);

  gsl_matrix_set(defs2,2,0,1);
  gsl_matrix_set(defs2,2,1,1);

  gsl_matrix_set(defs2,3,0,1);
  gsl_matrix_set(defs2,3,1,-1);

  gsl_matrix_set(defs2,4,0,1);
  gsl_matrix_set(defs2,4,1,0);  
  gsl_matrix_set(defs2,4,2,1);

  gsl_matrix_set(defs2,5,0,1);
  gsl_matrix_set(defs2,5,1,0);  
  gsl_matrix_set(defs2,5,2,-1);

  gsl_matrix_set(defs2,6,0,1);
  gsl_matrix_set(defs2,6,1,1);  
  gsl_matrix_set(defs2,6,2,1);

  gsl_matrix_set(defs2,7,0,1);
  gsl_matrix_set(defs2,7,1,-1); 
  gsl_matrix_set(defs2,7,2,1);

  gsl_matrix_set(defs2,8,0,1);
  gsl_matrix_set(defs2,8,1,1);  
  gsl_matrix_set(defs2,8,2,-1);

  gsl_matrix_set(defs2,9,0,1); 
  gsl_matrix_set(defs2,9,1,-1); 
  gsl_matrix_set(defs2,9,2,-1); 

/**********************************/

  gsl_matrix_set(defs3,0,0,1);

  gsl_matrix_set(defs3,0,1,1);

  gsl_matrix_set(defs3,1,1,1);
  gsl_matrix_set(defs3,1,2,1);

  gsl_matrix_set(defs3,2,1,1);
  gsl_matrix_set(defs3,2,2,-1);

  gsl_matrix_set(defs3,3,0,1);
  gsl_matrix_set(defs3,3,1,0);  
  gsl_matrix_set(defs3,3,2,0);

  gsl_matrix_set(defs3,4,0,1);
  gsl_matrix_set(defs3,4,1,0);  
  gsl_matrix_set(defs3,4,2,1);

  gsl_matrix_set(defs3,5,0,1);
  gsl_matrix_set(defs3,5,1,0);  
  gsl_matrix_set(defs3,5,2,-1);

  gsl_matrix_set(defs3,6,0,1);
  gsl_matrix_set(defs3,6,1,1); 
  gsl_matrix_set(defs3,6,2,0);

  gsl_matrix_set(defs3,7,0,1);
  gsl_matrix_set(defs3,7,1,-1);  
  gsl_matrix_set(defs3,7,2,0);

  gsl_matrix_set(defs3,8,0,1); 
  gsl_matrix_set(defs3,8,1,1); 
  gsl_matrix_set(defs3,8,2,1); 

  gsl_matrix_set(defs3,9,0,1); 
  gsl_matrix_set(defs3,9,1,-1); 
  gsl_matrix_set(defs3,9,2,1); 

  gsl_matrix_set(defs3,10,0,1);  
  gsl_matrix_set(defs3,10,1,1);  
  gsl_matrix_set(defs3,10,2,-1); 

  gsl_matrix_set(defs3,11,0,1);  
  gsl_matrix_set(defs3,11,1,-1);  
  gsl_matrix_set(defs3,11,2,-1);  

/*                 N       CO      Hall                     */

/* 5000            0       0       1 */
/* 5001            0       1       0 */
/* 5002            0       1       1 */
/* 5003            0       1       -1 */
/* 5004            1       0       0 */
/* 5005            1       0       1 */
/* 5006            1       0       -1 */
/* 5007            1       1       0 */
/* 5008            1       -1      0 */
/* 5009            1       1       1 */
/* 5010            1       -1      1 */
/* 5011            1       1       -1 */
/* 5012            1       -1      -1 */


/********************************/

/*   gsl_matrix_set(defs2,0,0,1); */

/*   gsl_matrix_set(defs2,1,1,1); */

/*   gsl_matrix_set(defs2,2,0,1); */
/*   gsl_matrix_set(defs2,2,1,1); */

/*   gsl_matrix_set(defs2,3,0,1); */
/*   gsl_matrix_set(defs2,3,1,-1); */

/*   gsl_matrix_set(defs2,4,0,1); */
/*   gsl_matrix_set(defs2,4,1,0);   */
/*   gsl_matrix_set(defs2,4,2,1); */

/*   gsl_matrix_set(defs2,5,0,1); */
/*   gsl_matrix_set(defs2,5,1,0);   */
/*   gsl_matrix_set(defs2,5,2,-1); */

/*   gsl_matrix_set(defs2,6,0,1); */
/*   gsl_matrix_set(defs2,6,1,1);   */
/*   gsl_matrix_set(defs2,6,2,0); */
/*   gsl_matrix_set(defs2,6,3,1); */

/*   gsl_matrix_set(defs2,7,0,1); */
/*   gsl_matrix_set(defs2,7,1,-1);  */
/*   gsl_matrix_set(defs2,7,2,0);    */
/*   gsl_matrix_set(defs2,7,2,1); */

/*   gsl_matrix_set(defs2,8,0,1); */
/*   gsl_matrix_set(defs2,8,1,1);   */
/*   gsl_matrix_set(defs2,8,2,0 ); */
/*   gsl_matrix_set(defs2,8,2,-1); */

/*   gsl_matrix_set(defs2,9,0,1);  */
/*   gsl_matrix_set(defs2,9,1,-1);  */
/*   gsl_matrix_set(defs2,9,2,0);   */
/*   gsl_matrix_set(defs2,9,2,-1);  */

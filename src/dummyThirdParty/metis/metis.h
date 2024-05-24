#ifndef METIS_H
#define METIS_H 1

/* *** DUMMY VERSION of metis.h - this file should not be included if you have metis
 *     installed in the correct position in $WM_THIRD_PARTY_DIR - see
 *     decompositionMethods/metis/Make/options
 */

#warning "Dummy metis.h - gets included since it cannot find metis installation."

#define IDXTYPEWIDTH 32

/*------------------------------------------------------------------------
* Undefine the following #define in order to use short idxtype as the idxtype 
*-------------------------------------------------------------------------*/
#if IDXTYPEWIDTH == 32
  #define SCNIDX  SCNd32
  #define PRIIDX  PRId32

  typedef int32_t idxtype;
#elif IDXTYPEWIDTH == 64
  #define SCNIDX  SCNd64
  #define PRIIDX  PRId64

  typedef int64_t idxtype;
#else
  #error "Incorrect user-supplied value fo IDXTYPEWIDTH"
#endif


void METIS_WPartGraphRecursive(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt,
                   idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, float *tpwgts, 
                   idxtype *options, idxtype *edgecut, idxtype *part);
void METIS_PartGraphRecursive(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt,
                   idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, idxtype *options, 
                   idxtype *edgecut, idxtype *part);
void METIS_WPartGraphKway(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
                   idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, float *tpwgts, 
                   idxtype *options, idxtype *edgecut, idxtype *part); 
void METIS_PartGraphKway(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
                   idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, idxtype *options, 
                   idxtype *edgecut, idxtype *part); 



#endif

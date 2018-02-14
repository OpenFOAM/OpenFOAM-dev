/*----------------------------------
 *     User Defined Reader - checker
 *----------------------------------*/

/********************************************************************
 *
 *      ****************************************
 *       Copyright 2004 Computational Engineering International, Inc.
 *       All Rights Reserved.
 *
 *              Restricted Rights Legend
 *
 *       Use, duplication, or disclosure of this
 *       software and its documentation by the
 *       Government is subject to restrictions as
 *       set forth in subdivision [(b)(3)(ii)] of
 *       the Rights in Technical Data and Computer
 *       Software clause at 52.227-7013.
 *******************************************************************/

/*----------------------------------------------------------------------
 * MAJOR ROUTINES ACCESS:     (VERSION 2.00)    Gold_Userd   API
 *
 *        Get the name of the reader
 *        ==========================
 *        USERD_get_name_of_reader
 *        USERD_get_extra_gui_numbers    (optional)
 *        USERD_get_extra_gui_defaults   (optional)
 *
 *        Get the reader version
 *        ======================
 *        USERD_get_reader_version
 *
 *        Set the filenames, gather timeset and time info
 *        ===============================================
 *        USERD_set_extra_gui_data       (optional)
 *        USERD_set_server_number
 *        USERD_set_filenames
 *        USERD_get_number_of_timesets
 *        USERD_get_geom_timeset_number
 *
 *        for each timeset:
 *          USERD_get_timeset_description
 *          USERD_get_num_of_time_steps
 *          USERD_get_sol_times
 *
 *        USERD_set_time_set_and_step
 *
 *
 *        Gather variable and time info
 *        =============================
 *        USERD_get_number_of_variables
 *        USERD_get_gold_variable_info
 *
 *
 *        Get initial part building info
 *        ==============================
 *        USERD_set_time_set_and_step
 *        USERD_get_changing_geometry_status
 *        USERD_get_node_label_status
 *        USERD_get_element_label_status
 *        USERD_get_number_of_files_in_dataset
 *        USERD_get_dataset_query_file_info
 *        USERD_get_descrip_lines                 (geometry)
 *        USERD_get_number_of_model_parts
 *        USERD_get_gold_part_build_info
 *        USERD_get_ghosts_in_model_flag
 *        USERD_get_maxsize_info
 *        USERD_get_ghosts_in_block_flag          (if any ghost cells in model)
 *        USERD_get_model_extents    **OR**
 *             USERD_get_part_coords  **AND/OR**
 *             USERD_get_block_coords_by_component
 *
 *
 *
 *         Part Builder
 *         ============
 *
 *         both unstructured and structured
 *         --------------------------------
 *         USERD_set_time_set_and_step
 *
 *         if unstructured
 *         ---------------
 *         USERD_get_part_element_ids_by_type
 *         USERD_get_part_elements_by_type
 *
 *         if any nsided elements:
 *           USERD_get_nsided_conn
 *
 *         if any nfaced elements:
 *           USERD_get_nfaced_nodes_per_face
 *           USERD_get_nfaced_conn
 *
 *         USERD_get_part_coords
 *         USERD_get_part_node_ids
 *
 *         else if structured
 *         ------------------
 *         USERD_get_block_iblanking
 *         USERD_get_block_coords_by_component
 *         USERD_get_block_ghost_flags
 *         USERD_get_part_node_ids              (If node ids given)
 *         USERD_get_part_element_ids_by_type   (If element ids given)
 *
 *         both again
 *         ----------
 *         USERD_get_border_availability        (If border representation
 *         USERD_get_border_elements_by_type     is selected)
 *
 *         USERD_stop_part_building
 *
 *
 *         Changing geometry
 *         =================
 *
 *         changing coords only
 *         --------------------
 *         USERD_set_time_set_and_step
 *         USERD_get_descrip_lines
 *         USERD_get_part_coords
 *         USERD_get_block_coords_by_component
 *
 *         changing connectivity
 *         ---------------------
 *         USERD_set_time_set_and_step
 *         USERD_get_descrip_lines
 *         USERD_get_number_of_model_parts
 *         USERD_get_gold_part_build_info
 *         USERD_get_ghosts_in_model_flag
 *         USERD_get_ghosts_in_block_flag       (If any ghost cells in model)
 *         USERD_get_model_extents     **OR**
 *             USERD_get_part_coords  **AND/OR**
 *             USERD_get_block_coords_by_component
 *         USERD_get_part_element_ids_by_type
 *         USERD_get_part_elements_by_type
 *         USERD_get_part_coords
 *         USERD_get_part_node_ids
 *         USERD_get_block_iblanking
 *         USERD_get_block_coords_by_component
 *         USERD_get_block_ghost_flags          (If ghost cells in part)
 *         USERD_get_part_node_ids              (If node ids given)
 *         USERD_get_part_element_ids_by_type   (If element ids given)
 *
 *         USERD_get_border_availability        (If border representation
 *         USERD_get_border_elements_by_type     is selected)
 *
 *
 *         Loading Variables
 *         ==================
 *
 *         constants
 *         ---------
 *         USERD_set_time_set_and_step
 *         USERD_get_constant_val
 *
 *         scalars/vectors/tensors
 *         -----------------------
 *         USERD_get_description_lines
 *         USERD_set_time_set_and_step
 *         USERD_get_var_by_component
 *
 *
 *         Node or Element queries over time
 *         =================================
 *         USERD_get_var_value_at_specific
 *
 *
 *  At 2.03, added:
 *  ---------------
 *
 *         Materials
 *         =========
 *         USERD_get_number_of_material_sets
 *         USERD_get_matf_set_info
 *
 *         If any material sets in the model (calls once per material set)
 *           USERD_get_number_of_materials
 *           USERD_get_matf_var_info
 *
 *         For each element type of each part containing material ids
 *           USERD_size_matf_data
 *           USERD_load_matf_data
 *
 *
 *  At 2.04, added:
 *  ---------------
 *          USERD_get_uns_failed_params  - Sets params used in element failure
 *
 *
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------
 * MAJOR ROUTINES ACCESS:       (Version 1.00)  original API
 *
 *        Get the name of the reader
 *        ==========================
 *        USERD_get_name_of_reader
 *        USERD_get_extra_gui_numbers     (optional #define _EGS)
 *        USERD_get_extra_gui_defaults    (optional #define _EGS)
 *
 *        Set the filenames
 *        =================
 *        USERD_set_extra_gui_data        (optional #define _EGS)
 *        USERD_set_filenames
 *        USERD_get_number_of_time_steps
 *        USERD_get_solution_times
 *        USERD_set_time_step
 *
 *
 *        Gather variable and time info
 *        =============================
 *        USERD_get_number_of_variables
 *        USERD_get_variable_info
 *
 *
 *        Get initial part building info
 *        ==============================
 *        USERD_set_time_step
 *        USERD_get_changing_geometry_status
 *        USERD_get_node_label_status
 *        USERD_get_element_label_status
 *        USERD_get_number_of_files_in_dataset
 *        USERD_get_dataset_query_file_info
 *        USERD_get_description_lines                     (geometry)
 *        USERD_get_number_of_model_parts
 *        USERD_get_part_build_info
 *        USERD_get_number_global_nodes
 *        USERD_get_global_coords
 *        USERD_get_block_coords_by_component
 *
 *        Failure Info
 *        ============
 *        USERD_get_uns_failed_params
 *
 *
 *        Part Builder
 *        ============
 *        USERD_set_time_step
 *        USERD_get_global_coords
 *        USERD_get_global_node_ids
 *        USERD_get_element_connectivities_for_part
 *        USERD_get_element_ids_for_part
 *        USERD_get_block_iblanking
 *        USERD_get_block_coords_by_component
 *
 *        USERD_stop_part_building
 *
 *
 *        Changing geometry
 *        =================
 *
 *        changing coords only
 *        --------------------
 *        USERD_set_time_step
 *        USERD_get_global_coords
 *        USERD_get_block_coords_by_component
 *
 *        changing connectivity
 *        ---------------------
 *        USERD_set_time_step
 *        USERD_get_number_of_model_parts
 *        USERD_get_part_build_info
 *        USERD_get_number_global_nodes
 *        USERD_get_global_coords
 *        USERD_get_global_node_ids
 *        USERD_get_element_connectivities_for_part
 *        USERD_get_element_ids_for_part
 *        USERD_get_block_iblanking
 *        USERD_get_block_coords_by_component
 *
 *        Loading Variables
 *        =================
 *
 *        constants:
 *        ----------
 *        USERD_set_time_step
 *        USERD_get_constant_value
 *
 *        scalars:
 *        --------
 *        USERD_get_description_lines
 *        USERD_set_time_step
 *        USERD_get_scalar_values
 *        USERD_get_block_scalar_values
 *
 *        vectors:
 *        --------
 *        USERD_get_description_lines
 *        USERD_set_time_step
 *        USERD_get_vector_values
 *        USERD_get_block_vector_values_by_component
 *
 *
 *        Node or Element queries over time
 *        =================================
 *        USERD_get_variable_value_at_specific
 *
 *----------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#ifndef WIN32
#include <unistd.h>
#endif
#include <string.h>
#include <ctype.h>
#include <sys/types.h>
#ifndef CONVEX
#include <malloc.h>
#endif
#include <math.h>

#include "global_extern.h"

/* -----------------------------------------------
 *  If you wish to test out the Extra GUI stuff
 *     you need to uncomment this section
 *
 *   #define _EGS
 *
 * ----------------------------------------------- */

#if (!defined USERD_API_100)
#define GT_USERD_API_100
#endif

#if (!defined USERD_API_100 && !defined USERD_API_200)
#define GT_USERD_API_200
#endif

#if (!defined USERD_API_100 && !defined USERD_API_200 && !defined USERD_API_201 && !defined USERD_API_202)
#define GT_USERD_API_202
#endif

#if (!defined USERD_API_100 && !defined USERD_API_200 && !defined USERD_API_201 && !defined USERD_API_202 && !defined USERD_API_203)
#define GT_USERD_API_203
#endif


#define EOS '\0'

typedef struct {
  int  id;                /* part_id                                      */
  char desc[Z_BUFL];      /* description given in the file                */
  int  type;              /* Z_UNSTRUCTURED, Z_STRUCTURED, Z_IBLANKED     */
  int  ne[Z_MAXTYPE];     /* Number of elements per type (Z_UNSTRUCTURED) */
                          /* or ne[0] = I dimension (Z_STRUCTURED)        */
                          /*    ne[1] = J dimension                       */
                          /*    ne[2] = K dimension                       */
  int  nn;                /* Num of unstructured nodes (All_Local only)   */
  int  ghosts;            /* TRUE if ghost cells, FALSE otherwise         */
}BUILDINFO;

typedef struct {
  char  description[Z_BUFL];  /* description */
  char  filename[Z_BUFL];     /* real filename */
  char  ifilename[Z_BUFL];    /* imaginary filename */
  int   type;
  int   classify;
  int   complex;
  float freq;
  int   contran;
  int   timeset;
}VARINFO;


typedef struct {
  char name[12];          /* Z_POINT, Z_G_POINT, Z_BAR02, ... */
  int  con_len;           /* Number of nodes per element      */
}EINFO;


/* Global variables
 *-----------------*/
int Geom_status;
int Node_labels;
int Element_labels;
int Ghosts_in_model;
int Num_parts;
int Num_vars;
int Num_materials_sets;
BUILDINFO *Pbuild;
VARINFO *Varinfo;
char Version_number[Z_MAX_USERD_NAME];
int Num_time_sets;
int Num_time_steps;

/* ------------------
 * Extra GUI stuff
 * ------------------ */
int Num_toggles = 0;
int Num_pulldowns = 0;
int Num_fields = 0;
char **Toggle_title;
int *Toggle_default_status;
char **Pulldown_title;
int *Pulldown_number_in_list;
int *Pulldown_default_selection;
char ***Pulldown_item_strings;
char **Field_title;
char **Field_user_string;
int *Toggle_choice; /* user choice */
int *Pulldown_choice; /* user choice */

/* ---------------------------
 * Failed elements (API 2.04)
 * --------------------------- */
int Any_uns_failed_model_elems = FALSE;


#if (defined USERD_API_100 || defined USERD_API_200)
EINFO Elem_info[Z_MAXTYPE] = {"Z_POINT",1,
                              "Z_BAR02",2,
                              "Z_BAR03",3,
                              "Z_TRI03",3,
                              "Z_TRI06",6,
                              "Z_QUA04",4,
                              "Z_QUA08",8,
                              "Z_TET04",4,
                              "Z_TET10",10,
                              "Z_PYR05",5,
                              "Z_PYR13",13,
                              "Z_PEN06",6,
                              "Z_PEN15",15,
                              "Z_HEX08",8,
                              "Z_HEX20",20};
#elif defined USERD_API_201
EINFO Elem_info[Z_MAXTYPE] = {"Z_POINT",  1,
                              "Z_G_POINT",1,
                              "Z_BAR02",  2,
                              "Z_G_BAR02",2,
                              "Z_BAR03",  3,
                              "Z_G_BAR03",3,
                              "Z_TRI03",  3,
                              "Z_G_TRI03",3,
                              "Z_TRI06",  6,
                              "Z_G_TRI06",6,
                              "Z_QUA04",  4,
                              "Z_G_QUA04",4,
                              "Z_QUA08",  8,
                              "Z_G_QUA08",8,
                              "Z_TET04",  4,
                              "Z_G_TET04",4,
                              "Z_TET10",  10,
                              "Z_G_TET10",10,
                              "Z_PYR05",  5,
                              "Z_G_PYR05",5,
                              "Z_PYR13",  13,
                              "Z_G_PYR13",13,
                              "Z_PEN06",  6,
                              "Z_G_PEN06",6,
                              "Z_PEN15",  15,
                              "Z_G_PEN15",15,
                              "Z_HEX08",  8,
                              "Z_G_HEX08",8,
                              "Z_HEX20",  20,
                              "Z_G_HEX20",20};
#else
EINFO Elem_info[Z_MAXTYPE] = {"Z_POINT",  1,
                              "Z_G_POINT",1,
                              "Z_BAR02",  2,
                              "Z_G_BAR02",2,
                              "Z_BAR03",  3,
                              "Z_G_BAR03",3,
                              "Z_TRI03",  3,
                              "Z_G_TRI03",3,
                              "Z_TRI06",  6,
                              "Z_G_TRI06",6,
                              "Z_QUA04",  4,
                              "Z_G_QUA04",4,
                              "Z_QUA08",  8,
                              "Z_G_QUA08",8,
                              "Z_TET04",  4,
                              "Z_G_TET04",4,
                              "Z_TET10",  10,
                              "Z_G_TET10",10,
                              "Z_PYR05",  5,
                              "Z_G_PYR05",5,
                              "Z_PYR13",  13,
                              "Z_G_PYR13",13,
                              "Z_PEN06",  6,
                              "Z_G_PEN06",6,
                              "Z_PEN15",  15,
                              "Z_G_PEN15",15,
                              "Z_HEX08",  8,
                              "Z_G_HEX08",8,
                              "Z_HEX20",  20,
                              "Z_G_HEX20",20,
                              "Z_NSIDED",  1,   /* Not yet implemented */
                              "Z_G_NSIDED",1,   /* Not yet implemented */
                              "Z_NFACED",  1,   /* Not yet implemented */
                              "Z_G_NFACED",1};  /* Not yet implemented */
#endif


/* Prototypes
 *-----------*/
static int load_fail_defaults(void);
static int prelim_info(int *two_fields, int *any_extra_gui);
static int get_input(int set_server_number,
                     int use_playfile,
                     char playfile[Z_MAXFILENP],
                     int two_fields,
                     int any_extra_gui,
                     int *swapbytes);
static int time_info( void );
static int part_build_info(int geom_time_step);
static int variable_info( void );

#if (defined GT_USERD_API_100)
static int gold_part_builder(int geom_time_step);
static int gold_var_loader(int var_time_step);
#else
static int part_builder(int geom_time_step);
static int var_loader(int var_time_step);
#endif

#if (defined GT_USERD_API_100)
static int materials_info( void );
static int gold_materials_loader(int geom_time_step);
#endif

static int entity_querys(int var_time_step);
static int exercise_bkup( void );
static void usage( void );


/*=============
 * Main Routine
 *=============*/
#ifdef WIN32
int main(int argc, char *argv[])
#else
int main(int argc, char *argv[])
#endif
{
  /* Command line option variables
   *------------------------------*/
  int set_server_number = FALSE;
  int use_playfile = FALSE;
  char playfile[Z_MAXFILENP];
  FILE *fplay;
  int geom_time_step = 0;
  int var_time_step = 0;

  /* Other local variables
   *----------------------*/
  int i, j;
  int err;
  int two_fields;
  int any_extra_gui = FALSE;
  int swapbytes;
  int indx;

  /*----------------------------
   * Command argument processing
   *----------------------------*/
  fprintf(stderr,"\n");
  fprintf(stderr,"\n");
  fprintf(stderr,"********************************************\n");
  fprintf(stderr,"*  EnSight User Defined Reader Debug Tool  *\n");
  fprintf(stderr,"********************************************\n");
  fprintf(stderr,"\n");

  indx = 1;
  while(indx < argc) {

    if(!strcmp("-h",argv[indx])) {
      usage();
    }
    else if(!strcmp("-help",argv[indx])) {
      usage();
    }

    /* if you want to test the server number routines
     *
     * Use:
     * > checker -server_number
     *
     * You will then be prompted for the current and total
     * number of servers
     *-----------------------------------------------*/
    else if(!strcmp("-server_number",argv[indx])) {
      set_server_number = TRUE;
    }


    /* if you want to use a "playfile" instead of being prompted
     * for the data loader information
     *
     * Use:
     * > checker -p <playfile>
     *
     * This playfile should have 3 [or 4] lines:
     * line 1:   the path
     * line 2:   filename_1
     * line 3:   [filename_2]   (if two_fields is TRUE)
     * line 4:   0 or 1, for swapytes (0 is FALSE, 1 is TRUE)
     *
     * example (two_fields is FALSE, so only 3 lines):
     *
     *  /usr/scratch/stealth/bjn/dat/sp_gold/binary
     *  simple.case
     *  1
     *
     *------------------------------------------------------*/
    else if(!strcmp("-p",argv[indx])) {
      indx++;
      if((indx < argc) && (argv[indx][0] != '-')) {
        use_playfile = TRUE;
        memset(playfile,EOS,Z_MAXFILENP);
        strcpy(playfile,argv[indx]);
      }
      else {
        usage();
      }
    }

    /* if you want to specify the geometry timestep to test (default is step 0)
     *
     * Use:
     * > checker -gts #
     *
     * Where # is the step number (zero based)
     *-------------------------------------------------------------------------*/
    else if(!strcmp("-gts",argv[indx])) {
      indx++;
      if((indx < argc) && (argv[indx][0] != '-')) {
        geom_time_step = atoi(argv[indx]);
      }
      else {
        usage();
      }
    }

    /* if you want to specify the variable timestep to test (default is step 0)
     * (will use this step for the appropriate timeset of each variable)
     *
     * Use:
     * > checker -vts #
     *
     * Where # is the step number (zero based)
     *-------------------------------------------------------------------------*/
    else if(!strcmp("-vts",argv[indx])) {
      indx++;
      if((indx < argc) && (argv[indx][0] != '-')) {
        var_time_step = atoi(argv[indx]);
      }
      else {
        usage();
      }
    }
    else {
      usage();
    }

    indx++;
  }


  /*-------------------------------------------------------------
   *
   * Now start exercising EnSight
   *
   *--------------------------------------------------------------*/

  /*-----------------
   * Preliminary info
   *-----------------*/
  err = prelim_info(&two_fields,&any_extra_gui);
  if(err == Z_ERR) {
    fprintf(stderr,"Stopping because of error in prelim_info\n");
    exit(1);
  }


  /*------------------
   * User input needed
   *------------------*/
  err = get_input(set_server_number,
                  use_playfile,
                  playfile,
                  two_fields,
                  any_extra_gui,
                  &swapbytes);
  if(err == Z_ERR) {
    fprintf(stderr,"Stopping because of error in get_input\n");
    exit(1);
  }


  /*----------
   * Time info
   *----------*/
  err = time_info();
  if(err == Z_ERR) {
    fprintf(stderr,"Stopping because of error in time_info\n");
    exit(1);
  }


  /*----------------
   * Part build info
   *----------------*/
  err = part_build_info(geom_time_step);
  if(err == Z_ERR) {
    fprintf(stderr,"Stopping because of error in part_build_info\n");
    exit(1);
  }


  /*------------------
   * Get Variable Info
   *------------------*/
  err = variable_info();
  if(err == Z_ERR) {
    fprintf(stderr,"Stopping because of error in variable_info\n");
    exit(1);
  }


#if (defined GT_USERD_API_202)
  /*-------------------
   * Get Materials Info
   *-------------------*/
  err = materials_info();
  if(err == Z_ERR) {
    fprintf(stderr,"Stopping because of error in materials_info\n");
    exit(1);
  }
#endif

#if (defined GT_USERD_API_203)
  if (Z_ERR == load_fail_defaults()) {
    fprintf(stderr,"Stopping due to error in failed element flag routine\n");
    exit(1);
  }
#endif

  /*------------------------
   * Act like building parts
   *------------------------*/
  if(Num_parts > 0) {

#if (defined GT_USERD_API_100)
    err = gold_part_builder(geom_time_step);
#else
    err = part_builder(geom_time_step);
#endif
    if(err == Z_ERR) {
      fprintf(stderr,"Stopping because of error in part_builder\n");
      exit(1);
    }
    else {
      USERD_stop_part_building();
    }
  }


  /*---------------------------
   * Act like loading variables
   *---------------------------*/
  if(Num_vars > 0) {

#if (defined GT_USERD_API_100)
    err = gold_var_loader(var_time_step);
#else
    err = var_loader(var_time_step);
#endif
    if(err == Z_ERR) {
      fprintf(stderr,"Stopping because of error in var_loader\n");
      exit(1);
    }
  }

#if (defined GT_USERD_API_202)
  /*---------------------------
   * Act like loading materials
   *---------------------------*/
  if(Num_materials_sets > 0) {
    err = gold_materials_loader(geom_time_step);
    if(err == Z_ERR) {
      fprintf(stderr,"Stopping because of error in materials_loader\n");
      exit(1);
    }
  }
#endif



  /*----------------------------------------------------
   * See if can do node and/or element queries over time
   *----------------------------------------------------*/
  if(Num_parts > 0 &&
     Num_vars > 0) {
    err = entity_querys(var_time_step);
    if(err == Z_ERR) {
      fprintf(stderr,"Stopping because of error in entity_querys\n");
      exit(1);
    }
  }

  /*----------------------------------------
   * Call the bkup file once in save mode,
   * then again in restore mode - so someone
   * could debug if desired
   *----------------------------------------*/
  err = exercise_bkup();
  if(err == Z_ERR) {
    fprintf(stderr,"Stopping due to error in saving and/or restoring archive\n");
    exit(1);
  }

  /*-------------
   * Exit Routine
   *-------------*/
  fprintf(stderr,"\n----------------- exiting ---------------\n");

#if (defined GT_USERD_API_100)
  USERD_exit_routine();
#endif

  fprintf(stderr,"\n\n");
}

/*--------------
 * Usage routine
 *--------------*/
static void
usage( void )
{
  fprintf(stderr,"------------------------------------------------------------\n");
  fprintf(stderr,"USAGE: checker [-p pfile] [-server_number] [-gts #] [-vts #]\n");
  fprintf(stderr,"------------------------------------------------------------\n");
  fprintf(stderr,"  -h, -help       Prints out this USAGE text.\n");
  fprintf(stderr,"  -gts #          Specify the geometry times step to use.)\n");
  fprintf(stderr,"  -p pfile        Plays the checker playfile (pfile).\n");
  fprintf(stderr,"  -server_number  Cause servers numbers to be prompted for.\n");
  fprintf(stderr,"  -vts #          Specify the variable times step to use.)\n");
  fprintf(stderr,"\n");
  exit(1);
}




/*------------
 * prelim_info
 *------------*/
static int
prelim_info(int *two_fields, int *any_extra_gui)
{
  int err;
  char reader_name[Z_MAX_USERD_NAME];
  char release_number[Z_MAX_USERD_NAME];
  char description[Z_MAXFILENP];
  int i,j;

  *any_extra_gui = FALSE;

  /* Get the reader name
   *--------------------*/
  err = USERD_get_name_of_reader(reader_name,two_fields);
  if(err == Z_OK) {
    fprintf(stderr," Name of reader: %s\n",reader_name);
    if(*two_fields==1) {
      fprintf(stderr," two_fields:     TRUE\n");
    }
    else if(*two_fields==0){
      fprintf(stderr," two_fields:     FALSE\n");
    }
    else if(*two_fields < 0) {
      fprintf(stderr," two_fields:     -1 (optional string) \n");
    }
  }
  else {
    fprintf(stderr,"Error: Could not get name of reader\n");
    return(Z_ERR);
  }

  /* Get the Extra GUI stuff (optional)
   * ---------------------------------------------------------- */
#ifdef _EGS

  /* Get the Extra GUI numbers of toggles, pulldowns, & fields
   * ---------------------------------------------------------- */

 USERD_get_extra_gui_numbers(      &Num_toggles,
                                  &Num_pulldowns,
                                  &Num_fields );

 if ( Num_toggles > 0 || Num_pulldowns > 0 || Num_fields > 0 ) {


   *any_extra_gui = TRUE;

   if (Num_toggles>0) {
     Toggle_title = (char **) calloc(Num_toggles,sizeof(char*));
     if (Toggle_title == (char **)NULL) return(Z_ERR);
     for (i=0; i<Num_toggles; i++) {
       Toggle_title[i] = (char *) calloc(Z_LEN_GUI_TITLE_STR,sizeof(char));
       if ( Toggle_title[i] == (char *)NULL ) return(Z_ERR);
     }
     Toggle_default_status = (int *) calloc(Num_toggles,sizeof(int));
     Toggle_choice = (int *) calloc(Num_toggles,sizeof(int));
   }

   if (Num_pulldowns > 0) {
     Pulldown_title = (char **) calloc( Num_pulldowns , sizeof(char*) );
     if (Pulldown_title == (char **)NULL) return(Z_ERR);

     Pulldown_item_strings = (char ***) calloc( Num_pulldowns , sizeof(char**) );
     if (Pulldown_item_strings == (char ***)NULL) return(Z_ERR);

     for (i=0; i<Num_pulldowns; i++) {
       Pulldown_title[i] = (char *) calloc( Z_LEN_GUI_TITLE_STR , sizeof(char) );
       if ( Pulldown_title[i] == (char *)NULL ) return(Z_ERR);

       Pulldown_item_strings[i] = (char **) calloc( Z_MAX_NUM_GUI_PULL_ITEMS , sizeof(char *) );
       if (Pulldown_item_strings[i] == (char **)NULL) return(Z_ERR);

       for(j = 0; j < Z_MAX_NUM_GUI_PULL_ITEMS; j++) {
         Pulldown_item_strings[i][j] = (char *) calloc( Z_LEN_GUI_PULL_STR , sizeof(char) );
         if ( Pulldown_item_strings[i][j] == (char *)NULL ) return(Z_ERR);
       }
     }
     Pulldown_number_in_list = (int *) calloc(Num_pulldowns,sizeof(int));
     Pulldown_default_selection = (int *) calloc(Num_pulldowns,sizeof(int));
     Pulldown_choice = (int *) calloc(Num_pulldowns,sizeof(int));
   }

   if (Num_fields > 0) {
     Field_title = (char **) calloc(Num_fields,sizeof(char*));
     Field_user_string = (char **) calloc(Num_fields,sizeof(char*));
     if (Field_title == (char **) NULL) return(Z_ERR);
     for (i=0; i<Num_fields; i++) {
       Field_title[i] = (char *) calloc(Z_LEN_GUI_TITLE_STR,sizeof(char));
       if ( Field_title[i] == (char *)NULL) return(Z_ERR);
       Field_user_string[i] = (char *) calloc(Z_LEN_GUI_FIELD_STR,sizeof(char));
       if ( Field_user_string[i] == (char *)NULL) return(Z_ERR);
     }
   }


   err = USERD_get_extra_gui_defaults(
                                    Toggle_title,                 /* [num_toggles][Z_LEN_GUI_TITLE_STR] */
                                    Toggle_default_status,        /* [num_toggles] */
                                    Pulldown_title,               /* [num_pulldowns][Z_LEN_GUI_TITLE_STR] */
                                    Pulldown_number_in_list,      /* [num_pulldowns] */
                                    Pulldown_default_selection,   /* [num_pulldowns] */
                                    Pulldown_item_strings,        /* [num_pulldowns][Z_MAX_NUM_GUI_PULL_ITEMS][Z_LEN_GUI_PULL_STR] */
                                    Field_title,                  /* [num_fields][Z_LEN_GUI_TITLE_STR] */
                                    Field_user_string              /* [num_fields][Z_LEN_GUI_FIELD_STR] */
                                    );
   if (Z_ERR == err) return(Z_ERR);

   fprintf(stderr,"\n**********************************************\n");
   fprintf(stderr,"****          Extra GUI Information        ***\n");
   fprintf(stderr,"**********************************************\n\n");

   fprintf(stderr,"\nTOGGLE INFO: %d active toggles\n",Num_toggles);
   for (i=0; i<Num_toggles; i++) {
     fprintf(stderr,"Toggle Title %d : %s\n",i,Toggle_title[i]);
     fprintf(stderr,"Default status = %d \n",Toggle_default_status[i]);
   }

   fprintf(stderr,"\nPULLDOWN INFO: %d active pulldowns\n",Num_pulldowns);
   for (i=0; i<Num_pulldowns; i++) {
     fprintf(stderr,"Pulldown Title %d : %s\n", i , Pulldown_title[i] );
     for (j=0; j<Z_MAX_NUM_GUI_PULL_ITEMS; j++) {
       fprintf(stderr,"Pulldown_item %d : %s\n",j,Pulldown_item_strings[i][j]);
       if (strlen(Pulldown_item_strings[i][j]) == 0) {
         Pulldown_number_in_list[i] = j;
         break;
       }
     }
     fprintf(stderr,"Number of items in list: %d\n",Pulldown_number_in_list[i]);
     fprintf(stderr,"Default selection: %d\n\n",Pulldown_default_selection[i]);
   }
   fprintf(stderr,"\n");

   fprintf(stderr,"\nFIELDINFO: %d active fields\n",Num_fields);
   for (i=0; i<Num_fields; i++) {
     fprintf(stderr,"Field Title %d : %s\n",i,Field_title[i]);
     fprintf(stderr,"Field string %d: %s\n",i,Field_user_string[i]);
   }
   fprintf(stderr,"\n\n\n");
 }

#endif


#if (defined GT_USERD_API_100)

  /* Get the reader api used
   *------------------------*/
  err = USERD_get_reader_version(Version_number);
  if(err == Z_OK) {
    fprintf(stderr," API version:    %s\n",Version_number);
  }
  else {
    fprintf(stderr,"Error: Could not get reader api version\n");
    return(Z_ERR);
  }

  /* Get the reader release
   *-----------------------*/
  err = USERD_get_reader_release(release_number);
  if(err == Z_OK) {
    fprintf(stderr," Release:        %s\n",release_number);
  }
  else {
    fprintf(stderr,"Error: Could not get reader release\n");
    return(Z_ERR);
  }
#else
  fprintf(stderr," API version:    1.00\n");
#endif


#if 0
  /* Get the reader description
   *---------------------------*/
  err = USERD_get_reader_descrip(description);
  if(err == Z_OK) {
    fprintf(stderr," Description:\n\n");
    fprintf(stderr,"%s\n\n",description);
  }
  else {
    fprintf(stderr,"Error: Could not get reader description\n");
    return(Z_ERR);
  }
#else
  fprintf(stderr,"  Note: Not currently calling USERD_get_reader_descrip\n");
  fprintf(stderr,"        because it is optional.\n");
#endif

  return(Z_OK);
}


/*----------
 * get_input
 *----------*/
static int
get_input(int set_server_number,
          int use_playfile,
          char playfile[Z_MAXFILENP],
          int two_fields,
          int any_extra_gui,
          int *swapbytes)
{
  FILE *fplay;

  int i, j;
  int err;
  int tot_servers;
  int cur_server;
  char the_path[Z_MAXFILENP];
  char file1[Z_MAXFILENP];
  char file2[Z_MAXFILENP];
  char filename_1[Z_MAXFILENP];
  char filename_2[Z_MAXFILENP];


  fprintf(stderr,"\n-------------- get_input ----------------\n");

  /*-----------------------------------------------------
   * Prompt for the two input values, as the client would
   * And set this info for the reader
   *-----------------------------------------------------*/

#if (defined GT_USERD_API_100)

  /* Set the server number - if command line option to do so
   *--------------------------------------------------------*/
  if(set_server_number) {
    fprintf(stderr,"     Enter total number of servers: ");
    scanf("%d",&tot_servers);

    fprintf(stderr,"     Enter current server number: ");
    scanf("%d",&cur_server);

    fprintf(stderr," Setting %d of %d for server number\n",cur_server,tot_servers);
    USERD_set_server_number(cur_server,tot_servers);
  }
#endif

  /* Set the filenames
   *------------------*/
  memset(the_path,EOS,Z_MAXFILENP);
  memset(file1,EOS,Z_MAXFILENP);
  memset(file2,EOS,Z_MAXFILENP);
  memset(filename_1,EOS,Z_MAXFILENP);
  memset(filename_2,EOS,Z_MAXFILENP);


  if(!use_playfile) {
    fprintf(stderr,"     Enter the path: ");
    scanf("%s",the_path);


    fprintf(stderr,"     Enter filename_1: ");
    scanf("%s",file1);

    if(two_fields == TRUE) {
      fprintf(stderr,"     Enter filename_2: ");
      scanf("%s",file2);
    }

    fprintf(stderr,"     Enter Swapbytes (0 if FALSE, 1 if TRUE): ");
    scanf("%d",swapbytes);

    if (TRUE == any_extra_gui ) {
      fprintf(stderr,"\n**********************************************\n");
      fprintf(stderr,"****          Extra GUI INPUT                ***\n");
      fprintf(stderr,"**********************************************\n\n");

      fprintf(stderr, "\n      TOGGLE INPUT \n");
      for (i=0; i<Num_toggles; i++) {
        fprintf(stderr, "      Enter Toggle Value for '%s' (1=toggle on, 0=toggle off)\n",Toggle_title[i]);
        scanf("%d",&Toggle_choice[i]);
      }
      fprintf(stderr, "\n      PULLDOWN INPUT \n");
      for (i=0; i<Num_pulldowns; i++) {
        fprintf(stderr, "\n      PULLDOWN # %d \n",i);
        for (j = 0; j<Pulldown_number_in_list[i]; j++) {
          fprintf(stderr, "              %d %s\n",j,Pulldown_item_strings[i][j]);
        }
        fprintf(stderr, "              Enter Pulldown Value for '%s' (0 to %d)\n",Pulldown_title[i],Pulldown_number_in_list[i]-1);
        scanf("%d",&Pulldown_choice[i]);
      }
      fprintf(stderr, "\n      FIELD INPUT \n");
      for (i=0; i<Num_fields; i++) {
        fprintf(stderr, "Enter string for field %d '%s'\n",i,Field_title[i]);
        scanf("%s",Field_user_string[i]);
      }

    }                /* end if there is any extra gui stuff */
  }                  /* end if not using playfile */
  else {
    fplay = fopen(playfile,"rb");
    if(fplay == (FILE *)NULL) {
      fprintf(stderr,"Error: Opening the playfile %s\n",playfile);
      return(Z_ERR);
    }
    else {
      fscanf(fplay,"%s",the_path);
      fscanf(fplay,"%s",file1);
      if(two_fields == TRUE) {
        fscanf(fplay,"%s",file2);
      }
      fscanf(fplay,"%d",swapbytes);

      /* ---------------------
       * Extra GUI stuff
       * --------------------- */
      if (TRUE == any_extra_gui) {

        for (i=0; i<Num_toggles; i++) {
          fscanf(fplay,"%d",&Toggle_choice[i]);
        }

        for (i=0; i<Num_pulldowns; i++) {
          fscanf(fplay,"%d",&Pulldown_choice[i]);
        }

        for (i=0; i<Num_fields; i++) {
          fscanf(fplay,"%s",Field_user_string[i]);
        }
      }
      fclose(fplay);
    }
  }

#ifdef _EGS
  /* -------------------------------------------
   * set the user choices here and run the code
   * ------------------------------------------- */

  /* set your choices here
     Toggle_choice[0..Num_toggles]
     Pulldown_choice[0..Num_pulldowns]
     Field_user_string[Num_fields][0..Numfields]
     amd then send your choices into this routine */

  USERD_set_extra_gui_data(
                  Toggle_choice,            /* [num_toggle] */
                  Pulldown_choice,          /* [num_pulldown] */
                  Field_user_string  );    /* [num_fields][Z_LEN_GUI_FIELD_STR] */

  for (i=0; i<Num_toggles; i++) {
    fprintf(stderr,"Toggle Title %d : %s\n",i,Toggle_title[i]);
    fprintf(stderr,"User selection = %d \n",Toggle_choice[i]);
  }
  fprintf(stderr,"\n\n");

  for (i=0; i<Num_pulldowns; i++) {
    fprintf(stderr,"Pulldown Title %d : %s\n", i , Pulldown_title[i] );
    fprintf(stderr,"Pulldown selection is # %d : %s\n",Pulldown_choice[i],Pulldown_item_strings[i][Pulldown_choice[i]]);
  }

  for (i=0; i<Num_fields; i++) {
    fprintf(stderr,"Field Title %d : %s\n",i,Field_title[i]);
    fprintf(stderr,"Field string %d: %s\n",i,Field_user_string[i]);

  }


#endif

  if(strncmp(file1,"/",1)) {
    strcpy(filename_1,the_path);
    strcat(filename_1,"/");
    strcat(filename_1,file1);
  }
  if(two_fields == TRUE) {
    if(strncmp(file2,"/",1)) {
      strcpy(filename_2,the_path);
      strcat(filename_2,"/");
      strcat(filename_2,file2);
    }
  }
  if(*swapbytes == 0) {
    *swapbytes = FALSE;
  }
  else {
    *swapbytes = TRUE;
  }

  /* Feedback
   *---------*/
  fprintf(stderr," path: %s\n",the_path);
  fprintf(stderr," filename_1: %s\n",filename_1);
  fprintf(stderr," filename_2: %s\n",filename_2);
  if(*swapbytes) {
    fprintf(stderr," Swapbytes:    TRUE\n");
  }
  else {
    fprintf(stderr," Swapbytes:    FALSE\n");
  }

  err = USERD_set_filenames(filename_1,filename_2,the_path,*swapbytes);
  if(err == Z_ERR) {
    fprintf(stderr,"Error: Trouble setting the filenames\n");
    return(Z_ERR);
  }

  return(Z_OK);
}


/*----------
 * time_info
 *----------*/
static int
time_info( void )
{
  int i;
  int err;
  int geom_time_set;
  int ts;
  float *sol_times;
  char ts_desc[Z_BUFL];

  fprintf(stderr,"\n-------------- time_info ----------------\n");

#if (defined GT_USERD_API_100)

  /* Get the number of timesets
   *---------------------------*/
  Num_time_sets = USERD_get_number_of_timesets();
  fprintf(stderr," number of timesets: %d\n",Num_time_sets);
  if(Num_time_sets == 0) {
    fprintf(stderr," So, static geometry and variables\n");
    return(Z_OK);
  }

  /* Get the timeset used for the geometry
   *--------------------------------------*/
  geom_time_set = USERD_get_geom_timeset_number();

  fprintf(stderr," geom timeset number: %d\n",geom_time_set);
  if(geom_time_set < 1 && Num_time_sets > 0) {
    fprintf(stderr,"Error: timeset numbers must be 1 or greater\n");
    fprintf(stderr,"       (unless Num_time_sets is zero also)\n");
  }


  /* For each timeset
   *-----------------*/
  for(ts=1; ts<=Num_time_sets; ++ts) {

    fprintf(stderr," Timeset %d:\n",ts);

    /* Get the timeset descriptions
     *-----------------------------*/
    err = USERD_get_timeset_description(ts,ts_desc);
    if(err == Z_ERR) {
      fprintf(stderr,"Error: getting timeset description\n");
      return(Z_ERR);
    }
    else {
      fprintf(stderr,"   description: %s\n",ts_desc);
    }

    /* Get the number of time steps
     *-----------------------------*/
    Num_time_steps = USERD_get_num_of_time_steps(ts);
    fprintf(stderr,"   number of time steps: %d\n",Num_time_steps);
    if(Num_time_steps < 1) {
      fprintf(stderr," Error: Number of time steps returned: %d\n",Num_time_steps);
      fprintf(stderr," (Must be >0 to be okay)\n");
      return(Z_ERR);
    }


    /* Get the solution times
     *-----------------------*/
    if(Num_time_steps > 0) {
      sol_times = (float *) calloc(Num_time_steps,sizeof(float));
      if(sol_times == (float *)NULL) {
        fprintf(stderr,"Error: allocating for solution times\n");
        return(Z_ERR);
      }
      else {
        err = USERD_get_sol_times(ts,sol_times);
        if(err == Z_ERR) {
          fprintf(stderr,"Error: getting solution times\n");
          return(Z_ERR);
        }
        else {
          for(i=0; i<Num_time_steps; ++i) {
            fprintf(stderr,"   At step %d, time = %f\n",i,sol_times[i]);
          }
        }
      }
      free(sol_times);
    }
  }

#else


  /* Get the number of time steps
   *-----------------------------*/
  Num_time_steps = USERD_get_number_of_time_steps();
  fprintf(stderr," Nnumber of time steps: %d\n",Num_time_steps);
  if(Num_time_steps < 1) {
    fprintf(stderr," Error: Number of time steps returned: %d\n",Num_time_steps);
    fprintf(stderr," (Must be >0 to be okay)\n");
    return(Z_ERR);
  }


  /* Get the solution times
   *-----------------------*/
  if(Num_time_steps > 0) {
    sol_times = (float *) calloc(Num_time_steps,sizeof(float));
    if(sol_times == (float *)NULL) {
      fprintf(stderr,"Error: allocating for solution times\n");
      return(Z_ERR);
    }
    else {
      err = USERD_get_solution_times(sol_times);
      if(err == Z_ERR) {
        fprintf(stderr,"Error: getting solution times\n");
        return(Z_ERR);
      }
      else {
        for(i=0; i<Num_time_steps; ++i) {
          fprintf(stderr,"   At step %d, time = %f\n",i,sol_times[i]);
        }
      }
    }
    free(sol_times);
  }

#endif

  return(Z_OK);
}



/*----------------
 * part_build_info
 *----------------*/
static int
part_build_info(int geom_time_step)
{
  int i, j;
  int fn;
  int err;
  int num_dataset_files;
  int geom_time_set;
  Z_QFILES *qfiles;
  char line1[Z_BUFL];
  char line2[Z_BUFL];

  int *part_ids;
  int *part_types;
  int *number_of_nodes;
  int **num_elems;
  int **ijk_dimensions;
  int **iblanking_options;
  char **part_descriptions;

  int ghosts_in_block;

  int *max_num_nodes;
  int **max_num_elems;
  int **max_ijk_dimensions;
  float extents[6];


  fprintf(stderr,"\n------------ part_build_info ------------\n");

#if (defined GT_USERD_API_100)

  /* Get the timeset used for the geometry
   *--------------------------------------*/
  geom_time_set = USERD_get_geom_timeset_number();

  /* Set the timeset and step - to first step
   *-----------------------------------------*/

  USERD_set_time_set_and_step(geom_time_set,geom_time_step);

#else

  /* Set the time step - to first step
   *----------------------------------*/
  USERD_set_time_step(geom_time_step);

#endif

  /* Get the changing geometry status
   *---------------------------------*/
  Geom_status = USERD_get_changing_geometry_status();
  if(Geom_status == Z_STATIC) {
    fprintf(stderr," Geom changing status: Z_STATIC\n");
  }
  else if(Geom_status == Z_CHANGE_COORDS) {
    fprintf(stderr," Geom changing status: Z_CHANGE_COORDS\n");
  }
  else if(Geom_status == Z_CHANGE_CONN) {
    fprintf(stderr," Geom changing status: Z_CHANGE_CONN\n");
  }
  else {
    fprintf(stderr," Invalid Geom changing status!!\n");
  }


  /* Get the node label status
   *--------------------------*/
  Node_labels = USERD_get_node_label_status();
  if(Node_labels) {
    fprintf(stderr," Node labels will be provided\n");
  }
  else {
    fprintf(stderr," Node labels will NOT be provided\n");
  }

  /* Get the element label status
   *-----------------------------*/
  Element_labels = USERD_get_element_label_status();
  if(Element_labels) {
    fprintf(stderr," Element labels will be provided\n");
  }
  else {
    fprintf(stderr," Element labels will NOT be provided\n");
  }

  fprintf(stderr,"\n");

  /* Get the number of files in the dataset
   *---------------------------------------*/
  num_dataset_files = USERD_get_number_of_files_in_dataset();
  fprintf(stderr," Number of files in dataset: %d\n",num_dataset_files);


  /* Get the dataset query file info
   *--------------------------------*/
  if(num_dataset_files > 0) {

    qfiles = (Z_QFILES *) calloc(num_dataset_files,sizeof(Z_QFILES));
    if(qfiles == (Z_QFILES *)NULL) {
      fprintf(stderr,"Error: allocating for dataset query files\n");
      return(Z_ERR);
    }
    else {

      for(i=0; i<num_dataset_files; ++i) {
        qfiles[i].f_desc = (char **) calloc(10,sizeof(char *));
        if(qfiles[i].f_desc == (char **)NULL) {
          fprintf(stderr,"Error: allocating for dataset query descrip lines\n");
          return(Z_ERR);
        }
        else {
          for(j=0; j<10; ++j) {
            qfiles[i].f_desc[j] = (char *) calloc(Z_MAXFILENP,sizeof(char));
            if(qfiles[i].f_desc[j] == (char *)NULL) {
              fprintf(stderr,"Error: allocating for dataset query descrip lines\n");
              return(Z_ERR);
            }
          }
        }
      }

      err = USERD_get_dataset_query_file_info(qfiles);
      if(err == Z_OK) {
        for(fn=0; fn<num_dataset_files; ++fn) {
          fprintf(stderr," For dataset file %d:\n",fn);

          fprintf(stderr,"   name:           %s\n",qfiles[fn].name);
          fprintf(stderr,"   size:           %d\n",qfiles[fn].sizeb);
          fprintf(stderr,"   time:           %s\n",qfiles[fn].timemod);
          fprintf(stderr,"   num desc lines: %d\n",qfiles[fn].num_d_lines);
          for(i=0; i<qfiles[fn].num_d_lines; ++i) {
            fprintf(stderr,"    desc line %d: %s\n",i,qfiles[fn].f_desc[i]);
          }
        }
      }
      else {
        fprintf(stderr,"Error: getting dataset query info\n");
        return(Z_ERR);
      }
    }

    /* Free allocated memory
     *----------------------*/
    for(i=0; i<num_dataset_files; ++i) {
      for(j=0; j<10; ++j) {
        free(qfiles[i].f_desc[j]);
      }
      free(qfiles[i].f_desc);
    }
    free(qfiles);
  }

  fprintf(stderr,"\n-----------------------------------------\n");

#if (defined GT_USERD_API_100)

  /* Get the geometry description lines
   *-----------------------------------*/
  err = USERD_get_descrip_lines(Z_GEOM,0,FALSE,line1,line2);
  if(err == Z_OK) {
    fprintf(stderr," Geom Desc line1: %s\n",line1);
    fprintf(stderr," Geom Desc line2: %s\n",line2);
  }
  else {
    fprintf(stderr,"Error: getting geom description lines\n");
    return(Z_ERR);
  }

#else

  /* Get the geometry description lines
   *-----------------------------------*/
  err = USERD_get_description_lines(Z_GEOM,0,line1,line2);
  if(err == Z_OK) {
    fprintf(stderr," Geom Desc line1: %s\n",line1);
    fprintf(stderr," Geom Desc line2: %s\n",line2);
  }
  else {
    fprintf(stderr,"Error: getting geom description lines\n");
    return(Z_ERR);
  }

#endif

  /* Get the number of model parts
   *------------------------------*/
  Num_parts = USERD_get_number_of_model_parts();
  if(Num_parts > 0) {
    fprintf(stderr," Number of parts: %d\n",Num_parts);
  }
  else {
    fprintf(stderr," Problems getting number of parts\n");
    return(Z_ERR);
  }



  /* Get the gold part build info
   *-----------------------------*/
  Pbuild = (BUILDINFO *) calloc(Num_parts,sizeof(BUILDINFO));
  if(Pbuild == (BUILDINFO *)NULL) {
    fprintf(stderr," Problems allocating for Pbuild structure\n");
    return(Z_ERR);
  }


  part_ids = (int *) calloc(Num_parts,sizeof(int));
  if(part_ids == (int *)NULL) {
    fprintf(stderr," Problems allocating for part ids\n");
    return(Z_ERR);
  }

  part_types = (int *) calloc(Num_parts,sizeof(int));
  if(part_types == (int *)NULL) {
    fprintf(stderr," Problems allocating for part types\n");
    return(Z_ERR);
  }

  part_descriptions = (char **) calloc(Num_parts,sizeof(char *));
  if(part_descriptions == (char **)NULL) {
    fprintf(stderr," Problems allocating for part descriptions\n");
    return(Z_ERR);
  }
  else {
    for(i=0; i<Num_parts; ++i) {
      part_descriptions[i] = (char *) calloc(Z_BUFL,sizeof(char));
      if(part_descriptions[i] == (char *)NULL) {
        fprintf(stderr," Problems allocating for part descriptions\n");
        return(Z_ERR);
      }
    }
  }

  number_of_nodes = (int *) calloc(Num_parts,sizeof(int));
  if(number_of_nodes == (int *)NULL) {
    fprintf(stderr," Problems allocating for part number of nodes\n");
    return(Z_ERR);
  }

  num_elems = (int **) calloc(Num_parts,sizeof(int *));
  if(num_elems == (int **)NULL) {
    fprintf(stderr," Problems allocating for part number of elements\n");
    return(Z_ERR);
  }
  else {
    for(i=0; i<Num_parts; ++i) {
      num_elems[i] = (int *) calloc(Z_MAXTYPE,sizeof(int));
      if(num_elems[i] == (int *)NULL) {
        fprintf(stderr," Problems allocating for part number of elements\n");
        return(Z_ERR);
      }
    }
  }

  ijk_dimensions = (int **) calloc(Num_parts,sizeof(int *));
  if(ijk_dimensions == (int **)NULL) {
    fprintf(stderr," Problems allocating for part ijk dimensions\n");
    return(Z_ERR);
  }
  else {

#if (defined GT_USERD_API_202)
    for(i=0; i<Num_parts; ++i) {
      ijk_dimensions[i] = (int *) calloc(9,sizeof(int));
      if(ijk_dimensions[i] == (int *)NULL) {
        fprintf(stderr," Problems allocating for part ijk dimensions\n");
        return(Z_ERR);
      }
      else {
        for(j=0; j<9; ++j) {
          ijk_dimensions[i][j] = -1;
        }
      }
    }
#else
    for(i=0; i<Num_parts; ++i) {
      ijk_dimensions[i] = (int *) calloc(3,sizeof(int));
      if(ijk_dimensions[i] == (int *)NULL) {
        fprintf(stderr," Problems allocating for part ijk dimensions\n");
        return(Z_ERR);
      }
    }
#endif
  }

  iblanking_options = (int **) calloc(Num_parts,sizeof(int *));
  if(iblanking_options == (int **)NULL) {
    fprintf(stderr," Problems allocating for part iblanking options\n");
    return(Z_ERR);
  }
  else {
    for(i=0; i<Num_parts; ++i) {
      iblanking_options[i] = (int *) calloc(6,sizeof(int));
      if(iblanking_options[i] == (int *)NULL) {
        fprintf(stderr," Problems allocating for part iblanking options\n");
        return(Z_ERR);
      }
    }
  }


#if (defined GT_USERD_API_100)

  err = USERD_get_gold_part_build_info(part_ids,
                                       part_types,
                                       part_descriptions,
                                       number_of_nodes,
                                       num_elems,
                                       ijk_dimensions,
                                       iblanking_options);
#else

  err = USERD_get_part_build_info(part_ids,
                                  part_types,
                                  part_descriptions,
                                  num_elems,
                                  ijk_dimensions,
                                  iblanking_options);

#endif

  if(err == Z_ERR) {
    fprintf(stderr," Problems getting part build info\n");
    return(Z_ERR);
  }
  else {
    for(i=0; i<Num_parts; ++i) {
      fprintf(stderr," For part %d:\n",i+1);

      fprintf(stderr,"   part id:   %d\n",part_ids[i]);
      Pbuild[i].id = part_ids[i];

      if(part_types[i] == Z_UNSTRUCTURED) {
        fprintf(stderr,"   part type: Z_UNSTRUCTURED\n");
      }
      else if(part_types[i] == Z_STRUCTURED) {
        fprintf(stderr,"   part type: Z_STRUCTURED\n");
      }
      else if(part_types[i] == Z_IBLANKED) {
        fprintf(stderr,"   part type: Z_IBLANKED\n");
      }
      else {
        fprintf(stderr,"   Invalid part type\n");
        return(Z_ERR);
      }
      Pbuild[i].type = part_types[i];

      fprintf(stderr,"   part desc:   %s\n",part_descriptions[i]);
      strncpy(Pbuild[i].desc,part_descriptions[i],Z_BUFL);

#if (defined GT_USERD_API_100)
      fprintf(stderr,"   number of nodes :   %d\n",number_of_nodes[i]);
      Pbuild[i].nn = number_of_nodes[i];
#else
      Pbuild[i].nn = USERD_get_number_of_global_nodes();
#endif

      for(j=0; j<Z_MAXTYPE; ++j) {
        if(num_elems[i][j] > 0) {
          fprintf(stderr,"   # %s elements:   %d\n",
                  Elem_info[j].name,num_elems[i][j]);
          Pbuild[i].ne[j] = num_elems[i][j];
        }
      }

      if(part_types[i] != Z_UNSTRUCTURED) {

        /* For this checker, we will place the following in the
         * Pbuild[].ne[] structure:
         *
         *   Note this is can be used for block size whether ranges or not
         *   -------------------------------------------------------------
         *   Pbuild[].ne[0] = i dim of current block (to the range selected)
         *   Pbuild[].ne[1] = j dim of current block (to the range selected)
         *   Pbuild[].ne[2] = k dim of current block (to the range selected)
         *
         *   Thus if ranges:
         *   ---------------
         *   Pbuild[].ne[3] = i min range          (-1 indicates no ranges)
         *   Pbuild[].ne[4] = i max range
         *   Pbuild[].ne[5] = i min range
         *   Pbuild[].ne[6] = i max range
         *   Pbuild[].ne[7] = i min range
         *   Pbuild[].ne[8] = i max range
         *
         *   Pbuild[].ne[9]  = i dim of total block (if ranges)
         *   Pbuild[].ne[10] = j dim of total block (if ranges)
         *   Pbuild[].ne[11] = k dim of total block (if ranges)
         *
         *   What comes back from the api is:
         *   --------------------------------
         *   before 2.03 (no ranges)
         *   -----------------------
         *   ijk_dimensions[][0] = i dim of block
         *   ijk_dimensions[][1] = j dim of block
         *   ijk_dimensions[][2] = k dim of block
         *
         *   at 2.03 (if no ranges)
         *   -------
         *   ijk_dimensions[][0] = i dim of block
         *   ijk_dimensions[][1] = j dim of block
         *   ijk_dimensions[][2] = k dim of block
         *   ijk_dimensions[][3] = -1
         *
         *   at 2.03 (if ranges)
         *   -------
         *   ijk_dimensions[][0] = i dim of total block
         *   ijk_dimensions[][1] = j dim of total block
         *   ijk_dimensions[][2] = k dim of total block
         *   ijk_dimensions[][3] = i min range
         *   ijk_dimensions[][4] = i max range
         *   ijk_dimensions[][5] = j min range
         *   ijk_dimensions[][6] = j max range
         *   ijk_dimensions[][7] = k min range
         *   ijk_dimensions[][8] = k max range
         *--------------------------------------------------------------*/

#if (defined GT_USERD_API_202)
        if(ijk_dimensions[i][3] == -1) {
          fprintf(stderr,"   ijk_dimensions: %d %d %d\n",
                  ijk_dimensions[i][0],
                  ijk_dimensions[i][1],
                  ijk_dimensions[i][2]);
          Pbuild[i].ne[0] = ijk_dimensions[i][0];
          Pbuild[i].ne[1] = ijk_dimensions[i][1];
          Pbuild[i].ne[2] = ijk_dimensions[i][2];
          Pbuild[i].ne[3] = ijk_dimensions[i][3];
        }
        else {

          /* If empty part
           *--------------*/
          if(ijk_dimensions[i][0] == 0 &&
             ijk_dimensions[i][1] == 0 &&
             ijk_dimensions[i][2] == 0) {
            fprintf(stderr,"   ijk_dimensions: %d %d %d\n",
                    ijk_dimensions[i][0],
                    ijk_dimensions[i][1],
                    ijk_dimensions[i][2]);
            Pbuild[i].ne[0] = ijk_dimensions[i][0];
            Pbuild[i].ne[1] = ijk_dimensions[i][1];
            Pbuild[i].ne[2] = ijk_dimensions[i][2];
            Pbuild[i].ne[3] = -1;
          }

          /* range part
           *-----------*/
          else {
            Pbuild[i].ne[0] = ijk_dimensions[i][4] - ijk_dimensions[i][3] + 1;
            Pbuild[i].ne[1] = ijk_dimensions[i][6] - ijk_dimensions[i][5] + 1;
            Pbuild[i].ne[2] = ijk_dimensions[i][8] - ijk_dimensions[i][7] + 1;

            Pbuild[i].ne[3] = ijk_dimensions[i][3];
            Pbuild[i].ne[4] = ijk_dimensions[i][4];
            Pbuild[i].ne[5] = ijk_dimensions[i][5];
            Pbuild[i].ne[6] = ijk_dimensions[i][6];
            Pbuild[i].ne[7] = ijk_dimensions[i][7];
            Pbuild[i].ne[8] = ijk_dimensions[i][8];

            Pbuild[i].ne[9] = ijk_dimensions[i][0];
            Pbuild[i].ne[10] = ijk_dimensions[i][1];
            Pbuild[i].ne[11] = ijk_dimensions[i][2];

            fprintf(stderr,"   Part has ranges:\n");
            fprintf(stderr,"   ijk dimensions of total block: %d %d %d\n",
                    Pbuild[i].ne[9],
                    Pbuild[i].ne[10],
                    Pbuild[i].ne[11]);
            fprintf(stderr,"     i range: %d  to  %d\n",
                    Pbuild[i].ne[3],
                    Pbuild[i].ne[4]);
            fprintf(stderr,"     j range: %d  to  %d\n",
                    Pbuild[i].ne[5],
                    Pbuild[i].ne[6]);
            fprintf(stderr,"     k range: %d  to  %d\n",
                    Pbuild[i].ne[7],
                    Pbuild[i].ne[8]);
            fprintf(stderr,"   ijk dimensions of range portion: %d %d %d\n",
                    Pbuild[i].ne[0],
                    Pbuild[i].ne[1],
                    Pbuild[i].ne[2]);
          }
        }
#else
        fprintf(stderr,"   ijk_dimensions: %d %d %d\n",
                ijk_dimensions[i][0],
                ijk_dimensions[i][1],
                ijk_dimensions[i][2]);
        Pbuild[i].ne[0] = ijk_dimensions[i][0];
        Pbuild[i].ne[1] = ijk_dimensions[i][1];
        Pbuild[i].ne[2] = ijk_dimensions[i][2];
        Pbuild[i].ne[3] = -1;
#endif
        if(part_types[i] == Z_IBLANKED) {
          fprintf(stderr,"   Ibanking options on:\n");
          if(iblanking_options[i][Z_EXT]) {
            fprintf(stderr,"     Z_EXT\n");
          }
          if(iblanking_options[i][Z_INT]) {
            fprintf(stderr,"     Z_INT\n");
          }
          if(iblanking_options[i][Z_BND]) {
            fprintf(stderr,"     Z_BND\n");
          }
          if(iblanking_options[i][Z_INTBND]) {
            fprintf(stderr,"     Z_INTBND\n");
          }
          if(iblanking_options[i][Z_SYM]) {
            fprintf(stderr,"     Z_SYM\n");
          }
        }
      }
    }
  }


#if (defined GT_USERD_API_200)

  /* Get ghosts in model flag
   *-------------------------*/
  Ghosts_in_model = USERD_get_ghosts_in_model_flag();
  if(Ghosts_in_model) {
    fprintf(stderr," Ghosts in Model:  TRUE\n");
  }
  else {
    fprintf(stderr," Ghosts in Model:  FALSE\n");
  }

  /* Get ghosts in block flag - if needed
   *-------------------------------------*/
  for(i=1; i<=Num_parts; ++i) {
    if(part_types[i-1] != Z_UNSTRUCTURED && Ghosts_in_model) {
      ghosts_in_block = USERD_get_ghosts_in_block_flag(i);
      Pbuild[i-1].ghosts = ghosts_in_block;
      if(ghosts_in_block) {
        fprintf(stderr," Ghosts in block part %d:  TRUE\n",i);
      }
      else {
        fprintf(stderr," Ghosts in block part %d:  FALSE\n",i);
      }
    }
  }

#endif


#if (defined GT_USERD_API_100)

  /* Get maxsize info
   *-----------------*/
  max_num_nodes = (int *) calloc(Num_parts,sizeof(int));
  if(max_num_nodes == (int *)NULL) {
    fprintf(stderr," Problems allocating for part max num of nodes\n");
    return(Z_ERR);
  }

  max_num_elems = (int **) calloc(Num_parts,sizeof(int *));
  if(max_num_elems == (int **)NULL) {
    fprintf(stderr," Problems allocating for part max num of elements\n");
    return(Z_ERR);
  }
  else {
    for(i=0; i<Num_parts; ++i) {
      max_num_elems[i] = (int *) calloc(Z_MAXTYPE,sizeof(int));
      if(max_num_elems[i] == (int *)NULL) {
        fprintf(stderr," Problems allocating for part max_num of elements\n");
        return(Z_ERR);
      }
    }
  }

  max_ijk_dimensions = (int **) calloc(Num_parts,sizeof(int *));
  if(max_ijk_dimensions == (int **)NULL) {
    fprintf(stderr," Problems allocating for part max ijk dimensions\n");
    return(Z_ERR);
  }
  else {
    for(i=0; i<Num_parts; ++i) {
      max_ijk_dimensions[i] = (int *) calloc(3,sizeof(int));
      if(max_ijk_dimensions[i] == (int *)NULL) {
        fprintf(stderr," Problems allocating for part max ijk dimensions\n");
        return(Z_ERR);
      }
    }
  }

  err = USERD_get_maxsize_info(max_num_nodes,
                               max_num_elems,
                               max_ijk_dimensions);
  if(err == Z_ERR) {
    fprintf(stderr," No maxsize info provided (or error getting them)\n");
  }
  else {

    for(i=0; i<Num_parts; ++i) {
      fprintf(stderr," For part %d:\n",i+1);

      fprintf(stderr,"   max number of nodes :   %d\n",max_num_nodes[i]);

      for(j=0; j<Z_MAXTYPE; ++j) {
        if(max_num_elems[i][j] > 0) {
          fprintf(stderr,"   max # %s elems:   %d\n",
                  Elem_info[j].name,max_num_elems[i][j]);
        }
      }

      if(part_types[i] != Z_UNSTRUCTURED) {
        fprintf(stderr,"   max_ijk_dimensions: %d %d %d\n",
                max_ijk_dimensions[i][0],
                max_ijk_dimensions[i][1],
                max_ijk_dimensions[i][2]);
      }
    }
  }

  /* Get model extents - if given
   *-----------------------------*/
  err = USERD_get_model_extents(extents);
  if(err == Z_ERR) {
    fprintf(stderr," No extents given\n");
  }
  else {
    fprintf(stderr," Min x: %g\n",extents[0]);
    fprintf(stderr," Max x: %g\n",extents[1]);
    fprintf(stderr," Min y: %g\n",extents[2]);
    fprintf(stderr," Max y: %g\n",extents[3]);
    fprintf(stderr," Min z: %g\n",extents[4]);
    fprintf(stderr," Max z: %g\n",extents[5]);
  }

#endif

  /* Free the allocated memory
   *--------------------------*/
  free(part_ids);
  free(part_types);
  free(number_of_nodes);

  for(i=0; i<Num_parts; ++i) {
    free(ijk_dimensions[i]);
    free(num_elems[i]);
    free(part_descriptions[i]);
  }
  free(ijk_dimensions);
  free(num_elems);
  free(iblanking_options);
  free(part_descriptions);

#if (defined GT_USERD_API_100)
  for(i=0; i<Num_parts; ++i) {
    free(max_ijk_dimensions[i]);
    free(max_num_elems[i]);
  }
  free(max_num_nodes);
  free(max_num_elems);
  free(max_ijk_dimensions);

#endif

  return(Z_OK);
}


/*--------------
 * variable_info
 *--------------*/
static int
variable_info( void )
{
  int i,j;
  int err;

  char **var_description;
  char **var_filename;
  int *var_type;
  int *var_classify;
  int *var_complex;
  char **var_ifilename;
  float *var_freq;
  int *var_contran;
  int *var_timeset;


  fprintf(stderr,"\n------------ variable_info --------------\n");

  /* Get the number of variables
   *----------------------------*/
  Num_vars = USERD_get_number_of_variables();
  if(Num_vars < 0) {
    fprintf(stderr,"Error: getting the number of variables\n");
  }
  else {
    fprintf(stderr," Number of variables: %d\n",Num_vars);
  }


  /* Get the gold variable info
   *---------------------------*/
  Varinfo = (VARINFO *) calloc(Num_vars,sizeof(VARINFO));
  if(Varinfo == (VARINFO *)NULL) {
    fprintf(stderr," Problems allocating for Varinfo structure\n");
    return(Z_ERR);
  }


  var_description = (char **) calloc(Num_vars,sizeof(char *));
  if(var_description == (char **)NULL) {
    fprintf(stderr," Problems allocating for var description\n");
    return(Z_ERR);
  }
  else {
    for(i=0; i<Num_vars; ++i) {
      var_description[i] = (char *) calloc(Z_BUFL,sizeof(char));
      if(var_description[i] == (char *)NULL) {
        fprintf(stderr," Problems allocating for var description\n");
        return(Z_ERR);
      }
    }
  }

  var_filename = (char **) calloc(Num_vars,sizeof(char *));
  if(var_filename == (char **)NULL) {
    fprintf(stderr," Problems allocating for var filename\n");
    return(Z_ERR);
  }
  else {
    for(i=0; i<Num_vars; ++i) {
      var_filename[i] = (char *) calloc(Z_BUFL,sizeof(char));
      if(var_filename[i] == (char *)NULL) {
        fprintf(stderr," Problems allocating for var filename\n");
        return(Z_ERR);
      }
    }
  }

  var_type = (int *) calloc(Num_vars,sizeof(int));
  if(var_type == (int *)NULL) {
    fprintf(stderr," Problems allocating for var type\n");
    return(Z_ERR);
  }

  var_classify = (int *) calloc(Num_vars,sizeof(int));
  if(var_classify == (int *)NULL) {
    fprintf(stderr," Problems allocating for var classify\n");
    return(Z_ERR);
  }

  var_complex = (int *) calloc(Num_vars,sizeof(int));
  if(var_complex == (int *)NULL) {
    fprintf(stderr," Problems allocating for var complex\n");
    return(Z_ERR);
  }


  var_ifilename = (char **) calloc(Num_vars,sizeof(char *));
  if(var_ifilename == (char **)NULL) {
    fprintf(stderr," Problems allocating for var ifilename\n");
    return(Z_ERR);
  }
  else {
    for(i=0; i<Num_vars; ++i) {
      var_ifilename[i] = (char *) calloc(Z_BUFL,sizeof(char));
      if(var_ifilename[i] == (char *)NULL) {
        fprintf(stderr," Problems allocating for var ifilename\n");
        return(Z_ERR);
      }
    }
  }

  var_freq = (float *) calloc(Num_vars,sizeof(float));
  if(var_freq == (float *)NULL) {
    fprintf(stderr," Problems allocating for var freq\n");
    return(Z_ERR);
  }

  var_contran = (int *) calloc(Num_vars,sizeof(int));
  if(var_contran == (int *)NULL) {
    fprintf(stderr," Problems allocating for var contran\n");
    return(Z_ERR);
  }

  var_timeset = (int *) calloc(Num_vars,sizeof(int));
  if(var_timeset == (int *)NULL) {
    fprintf(stderr," Problems allocating for var timeset\n");
    return(Z_ERR);
  }

#if (defined GT_USERD_API_100)

  err = USERD_get_gold_variable_info(var_description,
                                     var_filename,
                                     var_type,
                                     var_classify,
                                     var_complex,
                                     var_ifilename,
                                     var_freq,
                                     var_contran,
                                     var_timeset);
#else

  err = USERD_get_variable_info(var_description,
                                var_filename,
                                var_type,
                                var_classify);

#endif

  if(err == Z_ERR) {
    fprintf(stderr,"Error: getting variable info\n");
  }
  else {
    for(i=0; i<Num_vars; ++i) {

      /* Loading the global
       * (for use in other routines)
       *----------------------------*/
      strncpy(Varinfo[i].description,var_description[i],Z_BUFL);
      strncpy(Varinfo[i].filename,var_filename[i],Z_BUFL);
      strncpy(Varinfo[i].ifilename,var_ifilename[i],Z_BUFL);
      Varinfo[i].type     = var_type[i];
      Varinfo[i].classify = var_classify[i];
      Varinfo[i].complex  = var_complex[i];
      Varinfo[i].freq     = var_freq[i];
      Varinfo[i].contran  = var_contran[i];
      Varinfo[i].timeset  = var_timeset[i];

      /* Echo some feedback
       *-------------------*/
      fprintf(stderr," For Variable %d:\n",i+1);

      fprintf(stderr,"   var desc:      %s\n",var_description[i]);
      fprintf(stderr,"   var filename:  %s\n",var_filename[i]);

#if (defined GT_USERD_API_100)
      if(var_complex[i]) {
        fprintf(stderr,"   var complex:   TRUE\n");
        fprintf(stderr,"   var ifilename: %s\n",var_ifilename[i]);
        fprintf(stderr,"   var freq:      %g\n",var_freq[i]);
      }
      else {
        fprintf(stderr,"   var complex:   FALSE\n");
      }
#endif

      if(var_type[i] == Z_CONSTANT) {
        fprintf(stderr,"   var type:      Z_CONSTANT\n");

#if (defined GT_USERD_API_100)
        if(var_contran[i]) {
          fprintf(stderr,"   var contran:  TRUE\n");
        }
        else {
          fprintf(stderr,"   var contran:  FALSE\n");
        }
#endif

      }
      else if(var_type[i] == Z_SCALAR) {
        fprintf(stderr,"   var type:      Z_SCALAR\n");
      }
      else if(var_type[i] == Z_VECTOR) {
        fprintf(stderr,"   var type:      Z_VECTOR\n");
      }
      else if(var_type[i] == Z_TENSOR) {
        fprintf(stderr,"   var type:      Z_TENSOR\n");
      }
      else if(var_type[i] == Z_TENSOR9) {
        fprintf(stderr,"   var type:      Z_TENSOR9\n");
      }
      else {
        fprintf(stderr,"   Invalid var type\n");
        return(Z_ERR);
      }

      if(var_classify[i] == Z_PER_NODE) {
        fprintf(stderr,"   var classify:  Z_PER_NODE\n");
      }
      else if(var_classify[i] == Z_PER_ELEM) {
        fprintf(stderr,"   var classify:  Z_PER_ELEM\n");
      }
      else if(var_classify[i] != Z_CONSTANT) {
        fprintf(stderr,"   Invalid var classify\n");
        return(Z_ERR);
      }

#if (defined GT_USERD_API_100)
      fprintf(stderr,"   var timeset:   %d\n",var_timeset[i]);
#endif
    }
  }

  /* Free the allocated memory
   *--------------------------*/
  for(i=0; i<Num_vars; ++i) {
    free(var_description[i]);
    free(var_filename[i]);
    free(var_ifilename[i]);
  }
  free(var_description);
  free(var_filename);
  free(var_ifilename);
  free(var_type);
  free(var_classify);
  free(var_complex);
  free(var_freq);
  free(var_contran);
  free(var_timeset);

  return(Z_OK);
}


#if (defined GT_USERD_API_100)
/*------------------
 * gold_part_builder
 *------------------*/
static int
gold_part_builder(int geom_time_step)
{
  int i, j, k, jj, kk;
  int err;
  int geom_timeset;
  int p, pn;
  int et, ne;
  int *elemids;
  int **conns;
  int nn;
  int comp;
  int bdim[3];
  int ib[5];
  int num_ghosts;
  int num_dims;
  int cell_type;
  float mm[6];
  float **coords;
  int *nodeids;
  int *iblanking;
  int *ghost_flag;
  short *parent_type;
  int *parent_num;
  int num_elems[Z_MAXTYPE];
  CRD *crds;
  int bd1,bd2,bd3;
  int empty_part;
  int *pdata;
  int nsid_len;
  int *nsid_con;
  int nface_len;
  int *nface_con;
  int npf_len;
  int *npf_con;
  int maxcheck;
  int num_failed = 0;
  int *fail_flags = (int *) NULL;

  fprintf(stderr,"\n------------- part_builder --------------\n");

  if(Num_time_sets > 0) {
    /* Get the timeset used for the geometry
     *--------------------------------------*/
    geom_timeset = USERD_get_geom_timeset_number();

    /* Get the number of time steps for this timeset
     *----------------------------------------------*/
    Num_time_steps = USERD_get_num_of_time_steps(geom_timeset);
    if(Num_time_steps < 1) {
      fprintf(stderr," Error: Number of time steps returned: %d\n",Num_time_steps);
      fprintf(stderr," (Must be >0 to be okay)\n");
      return(Z_ERR);
    }
    if(geom_time_step > (Num_time_steps - 1)) {
      geom_time_step = Num_time_steps - 1;
    }

    /* Set the timeset and step - to first step by default, but
     * can set it at others using -gts command argument
     *---------------------------------------------------------*/
    USERD_set_time_set_and_step(geom_timeset,geom_time_step);

    fprintf(stderr," Using timeset:   %d  (step range is %d through %d)\n",
            geom_timeset,0,Num_time_steps-1);
    fprintf(stderr," Using time step: %d\n",geom_time_step);
  }

  for(p=0; p<Num_parts; ++p) {
    pn = p+1;

    fprintf(stderr,"\n\n----------------------------------------");
    fprintf(stderr," Part %d:\n",pn);

    /*-----------------------
     * For unstructured parts
     *-----------------------*/
    if(Pbuild[p].type == Z_UNSTRUCTURED) {

      for(et=0; et<Z_MAXTYPE; ++et) {
        ne = Pbuild[p].ne[et];

        if(ne > 0) {

          pdata = (int *)calloc(ne*Elem_info[et].con_len,sizeof(int));
          if(pdata == (int *) NULL) {
            fprintf(stderr,"Error: allocating conns array\n");
            return(Z_ERR);
          }
          else {
            conns = (int **) calloc(ne,sizeof(int *));
            if(conns == (int **) NULL) {
              fprintf(stderr,"Error: allocating conns array\n");
              return(Z_ERR);
            }
            for(i=0; i<ne; ++i) {
              conns[i] = pdata;
              pdata += Elem_info[et].con_len;
            }
          }


          /* Get the elements
           *-----------------*/
          err = USERD_get_part_elements_by_type(pn,et,conns);
          if(err == Z_ERR) {
            fprintf(stderr,"Error: getting element connectivities\n");
            return(Z_ERR);
          }

          if(Element_labels) {
            elemids = (int *) calloc(ne,sizeof(int));
            if(elemids == (int *) NULL) {
              fprintf(stderr,"Error: allocating elemids array\n");
              return(Z_ERR);
            }
          }

          /* Get the element ids - if any
           *-----------------------------*/
          if(Element_labels) {
            err = USERD_get_part_element_ids_by_type(pn,et,elemids);
            if(err == Z_ERR) {
              fprintf(stderr,"Error: getting element ids\n");
              return(Z_ERR);
            }
          }

          /* Echo "some" info
           *-----------------*/

#if (defined GT_USERD_API_202)

          maxcheck = Z_NSIDED;

          /* Nsided elements, if any
           *------------------------*/
          if(et == Z_NSIDED ||
             et == Z_G_NSIDED) {

            nsid_len = 0;
            for(i=0; i<ne; ++i) {
              nsid_len += conns[i][0];
            }

            nsid_con = (int *) calloc(nsid_len,sizeof(int));
            if(nsid_con == (int *) NULL) {
              fprintf(stderr,"Error: allocating nsided conn array\n");
              return(Z_ERR);
            }

            err = USERD_get_nsided_conn(pn,nsid_con);
            if(err == Z_ERR) {
              fprintf(stderr,"Error: getting nsided conn array\n");
              return(Z_ERR);
            }

            /* First element of the type
             *--------------------------*/
            i = 0;
            fprintf(stderr,"   %s Element %d of %d:\n",Elem_info[et].name,i+1,ne);
            if(Element_labels) {
              fprintf(stderr,"      id: %d\n",elemids[i]);
            }
            fprintf(stderr,"      connectivity:");
            for(j=0; j<conns[i][0]; ++j) {
              fprintf(stderr," %d",nsid_con[j]);
            }
            fprintf(stderr,"\n");

            /* Last element of the type
             *-------------------------*/
            i = ne - 1;
            if(i > 0) {
              fprintf(stderr,"   %s Element %d of %d:\n",Elem_info[et].name,i+1,ne);
              if(Element_labels) {
                fprintf(stderr,"      id: %d\n",elemids[i]);
              }
              fprintf(stderr,"      connectivity:");

              for(j=nsid_len-conns[i][0]; j<nsid_len; ++j) {
                fprintf(stderr," %d",nsid_con[j]);
              }
              fprintf(stderr,"\n");
            }
          }

          /* Nfaced elements if any
           *-----------------------*/
          if(et == Z_NFACED ||
             et == Z_G_NFACED) {

            nface_len = 0;
            for(i=0; i<ne; ++i) {
              nface_len += conns[i][0];
            }

            nface_con = (int *) calloc(nface_len,sizeof(int));
            if(nface_con == (int *) NULL) {
              fprintf(stderr,"Error: allocating nfaced face array\n");
              return(Z_ERR);
            }

            err = USERD_get_nfaced_nodes_per_face(pn,nface_con);
            if(err == Z_ERR) {
              fprintf(stderr,"Error: getting nfaced face array array\n");
              return(Z_ERR);
            }

            npf_len = 0;
            for(i=0; i<nface_len; ++i) {
              npf_len += nface_con[i];
            }

            npf_con = (int *) calloc(npf_len,sizeof(int));
            if(npf_con == (int *) NULL) {
              fprintf(stderr,"Error: allocating nfaced npf array\n");
              return(Z_ERR);
            }

            err = USERD_get_nfaced_conn(pn,npf_con);
            if(err == Z_ERR) {
              fprintf(stderr,"Error: getting nfaced conn array\n");
              return(Z_ERR);
            }

            /* First element of the type
             *--------------------------*/
            jj = 0;
            kk = 0;
            for(i=0; i<ne; ++i) {

              if(i == 0 ||
                 i == ne-1) {
                fprintf(stderr,"   %s Element %d of %d:\n",Elem_info[et].name,
                        i+1,ne);
                if(Element_labels) {
                  fprintf(stderr,"      id: %d\n",elemids[i]);
                }
                for(j=0; j<conns[i][0]; ++j) {
                  fprintf(stderr,"      face %d connectivity:",j+1);
                  for(k=0; k<nface_con[jj]; ++k) {
                    fprintf(stderr," %d",npf_con[kk]);
                    ++kk;
                  }
                  fprintf(stderr,"\n");
                  ++jj;
                }
              }
              else {
                for(j=0; j<conns[i][0]; ++j) {
                  for(k=0; k<nface_con[jj]; ++k) {
                    ++kk;
                  }
                  ++jj;
                }
              }
            }
          }
#else

          maxcheck = Z_MAXTYPE;

#endif

          /* Regular elements
           *-----------------*/
          if(et < maxcheck) {

            /* First element of the type
             *--------------------------*/
            i = 0;
            fprintf(stderr,"   %s Element %d of %d:\n",Elem_info[et].name,i+1,ne);
            if(Element_labels) {
              fprintf(stderr,"      id: %d\n",elemids[i]);
            }
            fprintf(stderr,"      connectivity:");
            for(j=0; j<Elem_info[et].con_len; ++j) {
              fprintf(stderr," %d",conns[i][j]);
            }
            fprintf(stderr,"\n");

            /* check the connectivity for negative numbers
             * -------------------------------------------*/
#if defined GT_USERD_API_100
            for (i=0;i<ne;i++){
              for(j=0; j<Elem_info[et].con_len; ++j) {
                /* ---------- uncomment to print out connectivity values ---------- */
/*              fprintf(stderr," %d",conns[i][j]);  */
                if (conns[i][j] <= 0 || conns[i][j] > Pbuild[p].nn ) {
                  fprintf(stderr,"\n****************************\n");
                  fprintf(stderr,"Connectivity value out of bounds: \n");
                  fprintf(stderr,"Either less than zero or greater than \n");
                  fprintf(stderr,"  number of nodes in part!!  \n");
                  fprintf(stderr,"i = %d   j = %d  conns[i][j] = %d \n",i,j,conns[i][j]);
                  fprintf(stderr,"****************************\n");
                }
              }
                /* ---------- uncomment to print out connectivity values ---------- */
/*            fprintf(stderr,"\n"); */
            }
#endif
            /* Last element of the type
             *-------------------------*/
            i = ne - 1;
            if(i > 0) {
              fprintf(stderr,"   %s Element %d of %d:\n",Elem_info[et].name,i+1,ne);
              if(Element_labels) {
                fprintf(stderr,"      id: %d\n",elemids[i]);
              }
              fprintf(stderr,"      connectivity:");
              for(j=0; j<Elem_info[et].con_len; ++j) {
                fprintf(stderr," %d",conns[i][j]);
              }
              fprintf(stderr,"\n");
            }
          }

          /* Free the allocated memory
           *--------------------------*/
          if(NULL != conns) {
            if(NULL != *conns) {
              free(*conns);
              *conns = NULL;
            }
            free(conns);
            conns = NULL;
          }

          if(Element_labels) {
            free(elemids);
          }
        }
      }

      /* Get the coords
       *---------------*/
      nn = Pbuild[p].nn;

      if(nn > 0) {

        coords = (float **) calloc(3,sizeof(float *));
        if(coords == (float **) NULL) {
          fprintf(stderr,"Error: allocating coords array\n");
          return(Z_ERR);
        }
        else {
          for(i=0; i<3; ++i) {
            coords[i] = (float *) calloc((nn+1),sizeof(float));
            if(coords[i] == (float *) NULL) {
              fprintf(stderr,"Error: allocating coords array\n");
              return(Z_ERR);
            }
          }
        }

        if(Node_labels) {
          nodeids = (int *) calloc((nn+1),sizeof(int));
          if(nodeids == (int *) NULL) {
            fprintf(stderr,"Error: allocating nodeids array\n");
            return(Z_ERR);
          }
        }


        err = USERD_get_part_coords(pn,coords);
        if(err == Z_ERR) {
          fprintf(stderr,"Error: getting unstructured coords\n");
          return(Z_ERR);
        }

        if(Node_labels) {
          err = USERD_get_part_node_ids(pn,nodeids);
          if(err == Z_ERR) {
            fprintf(stderr,"Error: getting nodeids\n");
            return(Z_ERR);
          }
        }

        /* Echo "some" info
         *-----------------*/

        /* First node
         *-----------*/
        i = 1;
        fprintf(stderr,"   Node %d of %d:\n",i,nn);
        if(Node_labels) {
          fprintf(stderr,"      id: %d\n",nodeids[i]);
        }
        fprintf(stderr,"      x y z coordinates: %g %g %g\n",
                coords[0][i], coords[1][i], coords[2][i]);
        mm[0] = mm[1] = coords[0][i];
        mm[2] = mm[3] = coords[1][i];
        mm[4] = mm[5] = coords[2][i];


        /* Last node
         *----------*/
        i = nn;
        if(i > 1) {
          fprintf(stderr,"   Node %d of %d:\n",i,nn);
          if(Node_labels) {
            fprintf(stderr,"      id: %d\n",nodeids[i]);
          }
          fprintf(stderr,"      x y z coordinates: %g %g %g\n",
                  coords[0][i], coords[1][i], coords[2][i]);
        }

        /* Min and max coordinate values
         *------------------------------*/
        for(i=2; i<=nn; ++i) {
          if(coords[0][i] < mm[0]) {
            mm[0] = coords[0][i];
          }
          if(coords[0][i] > mm[1]) {
            mm[1] = coords[0][i];
          }
          if(coords[1][i] < mm[2]) {
            mm[2] = coords[1][i];
          }
          if(coords[1][i] > mm[3]) {
            mm[3] = coords[1][i];
          }
          if(coords[2][i] < mm[4]) {
            mm[4] = coords[2][i];
          }
          if(coords[2][i] > mm[5]) {
            mm[5] = coords[2][i];
          }
        }

        fprintf(stderr,"   Coordinate ranges:\n");
        fprintf(stderr,"      min x: %g\n",mm[0]);
        fprintf(stderr,"      max x: %g\n",mm[1]);
        fprintf(stderr,"      min y: %g\n",mm[2]);
        fprintf(stderr,"      max y: %g\n",mm[3]);
        fprintf(stderr,"      min z: %g\n",mm[4]);
        fprintf(stderr,"      max z: %g\n",mm[5]);


        /* Free the allocated memory
         *--------------------------*/
        for(i=0; i<3; ++i) {
          free(coords[i]);
        }
        free(coords);
        if(Node_labels) {
          free(nodeids);
        }
      }
    }


    /*---------------------
     * For structured parts
     *---------------------*/
    else {

      empty_part = FALSE;
      if(Pbuild[p].ne[0] == 0 &&
         Pbuild[p].ne[1] == 0 &&
         Pbuild[p].ne[2] == 0) {
        empty_part = TRUE;
      }

      if(!empty_part) {

        /* Get the block coords
         *---------------------*/
        for(comp=0; comp<3; ++comp) {
          if(Pbuild[p].ne[comp] < 1) {
            bdim[comp] = 1;
          }
          else {
            bdim[comp] = Pbuild[p].ne[comp];
          }
        }
        nn = bdim[0] * bdim[1] * bdim[2];

        bd1 = bdim[0]-1;
        if(bd1 < 1) {
          bd1 = 1;
        }
        bd2 = bdim[1]-1;
        if(bd2 < 1) {
          bd2 = 1;
        }
        bd3 = bdim[2]-1;
        if(bd3 < 1) {
          bd3 = 1;
        }
        ne = bd1 * bd2 * bd3;

        /* Determine cell type
         *--------------------*/
        num_dims = 3;
        for(i=0; i<3; ++i) {
          if(bdim[i] == 1) {
            --num_dims;
          }
        }
        if(num_dims == 3) {
          cell_type = Z_HEX08;
        }
        else if(num_dims == 2) {
          cell_type = Z_QUA04;
        }
        else {
          cell_type = Z_BAR02;
        }

        coords = (float **) calloc(num_dims,sizeof(float *));
        if(coords == (float **) NULL) {
          fprintf(stderr,"Error: allocating coords array\n");
          return(Z_ERR);
        }
        else {
          for(i=0; i<num_dims; ++i) {
            coords[i] = (float *) calloc(nn,sizeof(float));
            if(coords[i] == (float *) NULL) {
              fprintf(stderr,"Error: allocating coords array\n");
              return(Z_ERR);
            }
          }
        }

        /* Get the coords
         *---------------*/
        for(comp=0; comp<num_dims; ++comp) {
          err = USERD_get_block_coords_by_component(pn,comp,coords[comp]);
          if(err == Z_ERR) {
            fprintf(stderr,"Error: getting block coords\n");
            return(Z_ERR);
          }
        }


#if (defined GT_USERD_API_200)

        if(Node_labels) {
          nodeids = (int *) calloc(nn,sizeof(int));
          if(nodeids == (int *) NULL) {
            fprintf(stderr,"Error: allocating nodeids array\n");
            return(Z_ERR);
          }
        }
        /* Get the node ids - if any
         *--------------------------*/
        if(Node_labels) {
          err = USERD_get_part_node_ids(pn,nodeids);
          if(err == Z_ERR) {
            fprintf(stderr,"Error: getting nodeids\n");
            return(Z_ERR);
          }
        }
#endif

        /* Echo "some" info
         *-----------------*/

        /* First node
         *-----------*/
        if(nn > 0) {
          i = 0;
          fprintf(stderr,"   Node %d of %d:\n",i+1,nn);

#if (defined GT_USERD_API_200)

          if(Node_labels) {
            fprintf(stderr,"      id: %d\n",nodeids[i]);
          }
#endif
          if(num_dims == 3) {
            fprintf(stderr,"      x y z coordinates: %g %g %g\n",
                    coords[0][i], coords[1][i], coords[2][i]);
            mm[0] = mm[1] = coords[0][i];
            mm[2] = mm[3] = coords[1][i];
            mm[4] = mm[5] = coords[2][i];
          }
          else if(num_dims == 2) {
            fprintf(stderr,"      x y coordinates: %g %g\n",
                    coords[0][i], coords[1][i]);
            mm[0] = mm[1] = coords[0][i];
            mm[2] = mm[3] = coords[1][i];
          }
          else {
            fprintf(stderr,"      x coordinates: %g\n",
                    coords[0][i]);
            mm[0] = mm[1] = coords[0][i];
          }


          /* Last node
           *----------*/
          i = nn-1;
          if(i > 1) {
            fprintf(stderr,"   Node %d of %d:\n",i+1,nn);

#if (defined GT_USERD_API_200)
            if(Node_labels) {
              fprintf(stderr,"      id: %d\n",nodeids[i]);
            }
#endif
            if(num_dims == 3) {
              fprintf(stderr,"      x y z coordinates: %g %g %g\n",
                      coords[0][i], coords[1][i], coords[2][i]);
            }
            else if(num_dims == 2) {
              fprintf(stderr,"      x y coordinates: %g %g\n",
                      coords[0][i], coords[1][i]);
            }
            else {
              fprintf(stderr,"      x coordinates: %g\n",
                      coords[0][i]);
            }
          }
        }

        /* Min and max coordinate values
         *------------------------------*/
        for(i=1; i<nn; ++i) {
          if(coords[0][i] < mm[0]) {
            mm[0] = coords[0][i];
          }
          if(coords[0][i] > mm[1]) {
            mm[1] = coords[0][i];
          }
          if(num_dims > 1) {
            if(coords[1][i] < mm[2]) {
              mm[2] = coords[1][i];
            }
            if(coords[1][i] > mm[3]) {
              mm[3] = coords[1][i];
            }
          }
          if(num_dims > 2) {
            if(coords[2][i] < mm[4]) {
              mm[4] = coords[2][i];
            }
            if(coords[2][i] > mm[5]) {
              mm[5] = coords[2][i];
            }
          }
        }

        fprintf(stderr,"   Coordinate ranges:\n");
        fprintf(stderr,"      min x: %g\n",mm[0]);
        fprintf(stderr,"      max x: %g\n",mm[1]);
        if(num_dims > 1) {
          fprintf(stderr,"      min y: %g\n",mm[2]);
          fprintf(stderr,"      max y: %g\n",mm[3]);
        }
        if(num_dims > 2) {
          fprintf(stderr,"      min z: %g\n",mm[4]);
          fprintf(stderr,"      max z: %g\n",mm[5]);
        }

        /* Free the allocated memory - so far
         *-----------------------------------*/
        for(i=0; i<num_dims; ++i) {
          free(coords[i]);
        }
        free(coords);

#if (defined GT_USERD_API_200)
        if(Node_labels) {
          free(nodeids);
        }
#endif

        /* Get the block iblanking - if any
         *---------------------------------*/
        if(Pbuild[p].type == Z_IBLANKED) {

          iblanking = (int *) calloc(nn,sizeof(int));
          if(iblanking == (int *) NULL) {
            fprintf(stderr,"Error: allocating iblanking array\n");
            return(Z_ERR);
          }

          err = USERD_get_block_iblanking(pn,iblanking);
          if(err == Z_ERR) {
            fprintf(stderr,"Error: getting block iblanking\n");
            return(Z_ERR);
          }

          /* Echo "some" info
           *-----------------*/
          ib[Z_EXT]    = 0;
          ib[Z_INT]    = 0;
          ib[Z_BND]    = 0;
          ib[Z_INTBND] = 0;
          ib[Z_SYM]    = 0;

          for(i=0; i<nn; ++i) {
            ++ib[iblanking[i]];
          }

          fprintf(stderr,"   Iblanking breakdown:\n");
          fprintf(stderr,"      Number of Z_EXT:    %d\n",ib[Z_EXT]);
          fprintf(stderr,"      Number of Z_INT:    %d\n",ib[Z_INT]);
          fprintf(stderr,"      Number of Z_BND:    %d\n",ib[Z_BND]);
          fprintf(stderr,"      Number of Z_INTBND: %d\n",ib[Z_INTBND]);
          fprintf(stderr,"      Number of Z_SYM:    %d\n",ib[Z_SYM]);

          free(iblanking);
        }

#if (defined GT_USERD_API_200)

        /* Get the ghost flags - if any
         *-----------------------------*/
        if(Pbuild[p].ghosts) {

          ghost_flag = (int *) calloc(ne,sizeof(int));
          if(ghost_flag == (int *) NULL) {
            fprintf(stderr,"Error: allocating ghost_flag array\n");
            return(Z_ERR);
          }

          err = USERD_get_block_ghost_flags(pn,ghost_flag);
          if(err == Z_ERR) {
            fprintf(stderr,"Error: getting block ghost flags\n");
            return(Z_ERR);
          }

          /* Echo "some" info
           *-----------------*/
          num_ghosts = 0;

          for(i=0; i<ne; ++i) {
            if(ghost_flag[i] > 0) {
              ++num_ghosts;
            }
          }

          fprintf(stderr,"   Block Ghost flag breakdown:\n");
          fprintf(stderr,"      %d ghost cells out of %d total cells\n",
                  num_ghosts,ne);

          free(ghost_flag);
        }

        /* Get the element ids - if any
         *-----------------------------*/
        if(Element_labels) {

          elemids = (int *) calloc(ne,sizeof(int));
          if(elemids == (int *) NULL) {
            fprintf(stderr,"Error: allocating elemids array\n");
            return(Z_ERR);
          }


          et = cell_type;
          err = USERD_get_part_element_ids_by_type(pn,et,elemids);
          if(err == Z_ERR) {
            fprintf(stderr,"Error: getting element ids\n");
            return(Z_ERR);
          }

          /* First element of the type
           *--------------------------*/
          i = 0;
          fprintf(stderr,"   %s Element %d of %d:\n",Elem_info[et].name,i+1,ne);
          fprintf(stderr,"      id: %d\n",elemids[i]);

          /* Last element of the type
           *-------------------------*/
          i = ne - 1;
          if(i > 0) {
            fprintf(stderr,"   %s Element %d of %d:\n",Elem_info[et].name,i+1,ne);
            fprintf(stderr,"      id: %d\n",elemids[i]);
          }

          free(elemids);
        }
#endif
      }
      else {
        fprintf(stderr,"   Empty structured part\n");
      }
    }

    /* Get border availability
     *------------------------*/
    err = USERD_get_border_availability(pn,num_elems);
    if(err == Z_OK) {

      /* Get border elements - if any
       *-----------------------------*/
      for(et=0; et<Z_MAXTYPE; ++et) {
        ne = num_elems[et];
        if(ne > 0) {

          conns = (int **) calloc(ne,sizeof(int *));
          if(conns == (int **) NULL) {
            fprintf(stderr,"Error: allocating border conns array\n");
            return(Z_ERR);
          }
          else {
            for(i=0; i<ne; ++i) {
              conns[i] = (int *) calloc(Elem_info[et].con_len,sizeof(int));
              if(conns[i] == (int *) NULL) {
                fprintf(stderr,"Error: allocating border conns array\n");
                return(Z_ERR);
              }
            }
          }

          parent_type = (short *) calloc(ne,sizeof(short));
          if(parent_type == (short *) NULL) {
            fprintf(stderr,"Error: allocating border parent_type array\n");
            return(Z_ERR);
          }

          parent_num = (int *) calloc(ne,sizeof(int));
          if(parent_num == (int *) NULL) {
            fprintf(stderr,"Error: allocating border parent_num array\n");
            return(Z_ERR);
          }


          err = USERD_get_border_elements_by_type(pn,
                                                  et,
                                                  conns,
                                                  parent_type,
                                                  parent_num);
          if(err == Z_ERR) {
            fprintf(stderr,"Error: getting border elements\n");
            return(Z_ERR);
          }


          /* Echo "some" info
           *-----------------*/

          /* First element of the type
           *--------------------------*/
          i = 0;
          fprintf(stderr,"   %s border element %d of %d:\n",
                  Elem_info[et].name,i+1,ne);
          fprintf(stderr,"      Parent type: %s\n",
                  Elem_info[parent_type[i]].name);
          fprintf(stderr,"      Parent num:  %d\n",parent_num[i]);
          fprintf(stderr,"      connectivity:");
          for(j=0; j<Elem_info[et].con_len; ++j) {
            fprintf(stderr," %d",conns[i][j]);
          }
          fprintf(stderr,"\n");

          /* Last element of the type
           *-------------------------*/
          i = ne - 1;
          if(i > 0) {
            fprintf(stderr,"   %s border element %d of %d:\n",
                    Elem_info[et].name,i+1,ne);
            fprintf(stderr,"      Parent type: %s\n",
                    Elem_info[parent_type[i]].name);
            fprintf(stderr,"      Parent num:  %d\n",parent_num[i]);
            fprintf(stderr,"      connectivity:");
            for(j=0; j<Elem_info[et].con_len; ++j) {
              fprintf(stderr," %d",conns[i][j]);
            }
            fprintf(stderr,"\n");
          }


          /* Free the allocated memory
           *--------------------------*/
          for(i=0; i<ne; ++i) {
            free(conns[i]);
          }
          free(conns);
          free(parent_type);
          free(parent_num);
        }
      }
    }
  } /* end for p = 0 to Num_parts */

  return(Z_OK);
}


/*----------------
 * gold_var_loader
 *----------------*/
static int
gold_var_loader(int var_time_step)
{
  int i, j;
  int err;
  int v, vn;
  int var_timeset;
  int p, pn;
  int et, e1, e2;
  int num_comps;
  int num_dims;
  int nsize;
  int comp;
  int bdim[3];
  int ne;
  int cell_type;
  float constant_val;
  char line1[Z_BUFL];
  char line2[Z_BUFL];
  float *values;
  float minv,maxv;
  int bd1,bd2,bd3;


  fprintf(stderr,"\n--------------- var_loader --------------\n");

  for(v=0; v<Num_vars; ++v) {
    vn = v + 1;

    if(v > 0) {
      fprintf(stderr,"\n");
    }
    if(Varinfo[v].classify == Z_PER_NODE) {
      fprintf(stderr," Z_PER_NODE Variable %d:\n",vn);
    }
    else {
      fprintf(stderr," Z_PER_ELEM Variable %d:\n",vn);
    }


    if(Num_time_sets > 0) {
      /* Get the timeset used for the variable
       *---------------------------------------*/
      var_timeset = Varinfo[v].timeset;

      /* Get the number of time steps for this timeset
       *----------------------------------------------*/
      Num_time_steps = USERD_get_num_of_time_steps(var_timeset);
      if(Num_time_steps < 1) {
        fprintf(stderr," Error: Number of time steps returned: %d\n",
                Num_time_steps);
        fprintf(stderr," (Must be >0 to be okay)\n");
        return(Z_ERR);
      }
      if(var_time_step > (Num_time_steps - 1)) {
        var_time_step = Num_time_steps - 1;
      }

      /* Set the timeset and step - to first step by default, but
       * can set it at others using -vts command argument
       *---------------------------------------------------------*/
      USERD_set_time_set_and_step(var_timeset,var_time_step);

      fprintf(stderr,"   Using timeset:   %d  (step range is %d through %d)\n",
              var_timeset,0,Num_time_steps-1);
      fprintf(stderr,"   Using time step: %d\n",var_time_step);
    }


    /* Constants
     *----------*/
    if(Varinfo[v].type == Z_CONSTANT) {

      constant_val = USERD_get_constant_val(vn,FALSE);
      fprintf(stderr,"   Constant (%s):\n",Varinfo[v].description);
      fprintf(stderr,"     value: %g\n",constant_val);

      if(Varinfo[v].complex) {
        constant_val = USERD_get_constant_val(vn,TRUE);
        fprintf(stderr,"   value (imag): %g\n",constant_val);
      }
    }

    /* Scalars, Vectors, Tensors
     *--------------------------*/
    else {

      /* Get the var description line
       *-----------------------------*/
      err = USERD_get_descrip_lines(Z_VARI,vn,FALSE,line1,line2);
      if(err == Z_OK) {
        fprintf(stderr,"   Desc line: %s\n",line1);
      }
      else {
        fprintf(stderr,"Error: getting var description line\n");
        return(Z_ERR);
      }

      if(Varinfo[v].complex) {
        err = USERD_get_descrip_lines(Z_VARI,vn,TRUE,line1,line2);
        if(err == Z_OK) {
          fprintf(stderr,"   Desc line (imag): %s\n",line1);
        }
        else {
          fprintf(stderr,"Error: getting var description line (imag)\n");
          return(Z_ERR);
        }
      }


      /* Get the values by component
       *-----------------------------*/
      if(Varinfo[v].type == Z_SCALAR) {
        num_comps = 1;
      }
      else if(Varinfo[v].type == Z_VECTOR) {
        num_comps = 3;
      }
      else if(Varinfo[v].type == Z_TENSOR) {
        num_comps = 6;
      }
      else if(Varinfo[v].type == Z_TENSOR9) {
        num_comps = 9;
      }


      /* Per_Node
       *---------*/
      if(Varinfo[v].classify == Z_PER_NODE) {

        for(p=0; p<Num_parts; ++p) {
          pn = p + 1;

          if(Pbuild[p].type == Z_UNSTRUCTURED) {
            nsize = Pbuild[p].nn;
          }
          else {
            for(comp=0; comp<3; ++comp) {
              if(Pbuild[p].ne[comp] < 1) {
                bdim[comp] = 1;
              }
              else {
                bdim[comp] = Pbuild[p].ne[comp];
              }
            }
            nsize = bdim[0] * bdim[1] * bdim[2];
          }


          fprintf(stderr,"   For part %d, with %d nodes:\n",pn,nsize);

          if(nsize > 0) {
            values = (float *) calloc((nsize+1),sizeof(float));
            if(values == (float *) NULL) {
              fprintf(stderr,"Error: allocating variable values\n");
              return(Z_ERR);
            }

            for(comp=0; comp<num_comps; ++comp) {

              err = USERD_get_var_by_component(vn,
                                               pn,
                                               Varinfo[v].type,
                                               0,
                                               FALSE,
                                               comp,
                                               values);
              if(err == Z_UNDEF) {
                fprintf(stderr,"  Variable not defined on this part\n");
              }


              /* For the component, show 1st node, last node, min, max values
               *-------------------------------------------------------------*/
              minv = maxv = values[1];
              for(i=2; i<=nsize; ++i) {
                if(values[i] < minv) {
                  minv = values[i];
                }
                if(values[i] > maxv) {
                  maxv = values[i];
                }
              }

              fprintf(stderr,"     For component %d: \n",comp);
              fprintf(stderr,"       node %10d value:   %g\n",1,values[1]);
              fprintf(stderr,"       node %10d value:   %g\n",nsize,values[nsize]);
              fprintf(stderr,"       min value:               %g\n",minv);
              fprintf(stderr,"       max value:               %g\n",maxv);

              if(Varinfo[v].complex) {
                err = USERD_get_var_by_component(vn,
                                                 pn,
                                                 Varinfo[v].type,
                                                 0,
                                                 FALSE,
                                                 comp,
                                                 values);
                if(err == Z_UNDEF) {
                  fprintf(stderr,"  Variable not defined on this part\n");
                }

                /* For the component, show 1st node, last node, min, max values
                 *-------------------------------------------------------------*/
                minv = maxv = values[1];
                for(i=2; i<=nsize; ++i) {
                  if(values[i] < minv) {
                    minv = values[i];
                  }
                  if(values[i] > maxv) {
                    maxv = values[i];
                  }
                }

                fprintf(stderr,"     For component %d (imag): \n",comp);
                fprintf(stderr,"       node %10d value:   %g\n",1,values[1]);
                fprintf(stderr,"       node %10d value:   %g\n",nsize,values[nsize]);
                fprintf(stderr,"       min value:               %g\n",minv);
                fprintf(stderr,"       max value:               %g\n",maxv);
              }
            }
            free(values);
          }
        }
      }

      /* Per_Elem
       *---------*/
      else {
        for(p=0; p<Num_parts; ++p) {
          pn = p + 1;

          if(Pbuild[p].type != Z_UNSTRUCTURED) {

            for(comp=0; comp<3; ++comp) {
              if(Pbuild[p].ne[comp] < 1) {
                bdim[comp] = 1;
              }
              else {
                bdim[comp] = Pbuild[p].ne[comp];
              }
            }

            bd1 = bdim[0]-1;
            if(bd1 < 1) {
              bd1 = 1;
            }
            bd2 = bdim[1]-1;
            if(bd2 < 1) {
              bd2 = 1;
            }
            bd3 = bdim[2]-1;
            if(bd3 < 1) {
              bd3 = 1;
            }
            nsize = bd1 * bd2 * bd3;


            /* Determine cell type
             *--------------------*/
            num_dims = 3;
            for(i=0; i<3; ++i) {
              if(bdim[i] == 1) {
                --num_dims;
              }
            }
            if(num_dims == 3) {
              cell_type = Z_HEX08;
            }
            else if(num_dims == 2) {
              cell_type = Z_QUA04;
            }
            else {
              cell_type = Z_BAR02;
            }
          }

          if(Pbuild[p].type == Z_UNSTRUCTURED) {
            e1 = 0;
            e2 = Z_MAXTYPE-1;
          }
          else {
            e1 = e2 = cell_type;
          }

          for(et=e1; et<=e2; ++et) {

            if(Pbuild[p].type == Z_UNSTRUCTURED) {
              nsize = Pbuild[p].ne[et];
            }

            if(nsize > 0) {

              fprintf(stderr,"   For part %d, with %d elems of type %s:\n",
                      pn,nsize,Elem_info[et].name);


              values = (float *) calloc((nsize+1),sizeof(float));
              if(values == (float *) NULL) {
                fprintf(stderr,"Error: allocating variable values\n");
                return(Z_ERR);
              }

              for(comp=0; comp<num_comps; ++comp) {

                err = USERD_get_var_by_component(vn,
                                                 pn,
                                                 Varinfo[v].type,
                                                 et,
                                                 FALSE,
                                                 comp,
                                                 values);
                if(err == Z_UNDEF) {
                  fprintf(stderr,"  Variable not defined on this part\n");
                }

                /* For the component, show 1st elem, last elem, min, max values
                 *-------------------------------------------------------------*/
                minv = maxv = values[1];
                for(i=2; i<=nsize; ++i) {
                  if(values[i] < minv) {
                    minv = values[i];
                  }
                  if(values[i] > maxv) {
                    maxv = values[i];
                  }
                }

                fprintf(stderr,"     For component %d: \n",comp);
                fprintf(stderr,"       elem %10d value:  %g\n",1,values[1]);
                fprintf(stderr,"       elem %10d value:  %g\n",nsize,values[nsize]);
                fprintf(stderr,"       min value:              %g\n",minv);
                fprintf(stderr,"       max value:              %g\n",maxv);

                if(Varinfo[v].complex) {
                  err = USERD_get_var_by_component(vn,
                                                   pn,
                                                   Varinfo[v].type,
                                                   et,
                                                   FALSE,
                                                   comp,
                                                   values);
                  if(err == Z_UNDEF) {
                    fprintf(stderr,"  Variable not defined on this part\n");
                  }

                  /* For the component, show 1st elem, last elem, min, max values
                   *-------------------------------------------------------------*/
                  minv = maxv = values[1];
                  for(i=2; i<=nsize; ++i) {
                    if(values[i] < minv) {
                      minv = values[i];
                    }
                    if(values[i] > maxv) {
                      maxv = values[i];
                    }
                  }

                  fprintf(stderr,"     For component %d (imag): \n",comp);
                  fprintf(stderr,"       elem %10d value:  %g\n",1,values[1]);
                  fprintf(stderr,"       elem %10d value:  %g\n",nsize,values[nsize]);
                  fprintf(stderr,"       min value:              %g\n",minv);
                  fprintf(stderr,"       max value:              %g\n",maxv);

                }
              }
              free(values);
            }
          }
        }
      }
    }
  }

  return(Z_OK);
}

#else

/*-------------
 * part_builder
 *-------------*/
static int
part_builder(int geom_time_step)
{
  int i, j;
  int err;
  int p, pn;
  int et, ne;
  int *elemids[Z_MAXTYPE];
  int **conns[Z_MAXTYPE];
  int nn;
  int comp;
  int bdim[3];
  int ib[5];
  int num_dims;
  int cell_type;
  float mm[6];
  float **coords;
  int *nodeids;
  int *iblanking;
  CRD *crds;
  int bd1,bd2,bd3;


  fprintf(stderr,"\n------------- part_builder --------------\n");


  if(Num_time_steps > 1) {
    if(geom_time_step > (Num_time_steps - 1)) {
      geom_time_step = Num_time_steps - 1;
    }

    /* Set the time step - to first step by default, but
     * can set it at others using -gts command argument
     *---------------------------------------------------*/
    USERD_set_time_step(geom_time_step);

    fprintf(stderr," Using time step: %d  (where range is %d through %d\n",
            geom_time_step,0,Num_time_steps-1);
  }


  /* Get the global coords
   *----------------------*/
  nn = USERD_get_number_of_global_nodes();

  if(nn > 0) {

    crds = (CRD *) calloc(nn,sizeof(CRD));
    if(crds == (CRD *) NULL) {
      fprintf(stderr,"Error: allocating crds array\n");
      return(Z_ERR);
    }

    if(Node_labels) {
      nodeids = (int *) calloc(nn,sizeof(int));
      if(nodeids == (int *) NULL) {
        fprintf(stderr,"Error: allocating nodeids array\n");
        return(Z_ERR);
      }
    }


    err = USERD_get_global_coords(crds);
    if(err == Z_ERR) {
      fprintf(stderr,"Error: getting unstructured coords\n");
      return(Z_ERR);
    }

    if(Node_labels) {
      err = USERD_get_global_node_ids(nodeids);
      if(err == Z_ERR) {
        fprintf(stderr,"Error: getting nodeids\n");
        return(Z_ERR);
      }
    }

    /* Echo "some" info
     *-----------------*/

    /* First node
     *-----------*/
    i = 0;
    fprintf(stderr,"   Node %d of %d:\n",i+1,nn);
    if(Node_labels) {
      fprintf(stderr,"      id: %d\n",nodeids[i]);
    }
    fprintf(stderr,"      x y z coordinates: %g %g %g\n",
            crds[i].xyz[0], crds[i].xyz[1], crds[i].xyz[2]);
    mm[0] = mm[1] = crds[i].xyz[0];
    mm[2] = mm[3] = crds[i].xyz[1];
    mm[4] = mm[5] = crds[i].xyz[2];


    /* Last node
     *----------*/
    i = nn-1;
    if(i > 0) {
      fprintf(stderr,"   Node %d of %d:\n",i+1,nn);
      if(Node_labels) {
        fprintf(stderr,"      id: %d\n",nodeids[i]);
      }
      fprintf(stderr,"      x y z coordinates: %g %g %g\n",
              crds[i].xyz[0], crds[i].xyz[1], crds[i].xyz[2]);
    }

    /* Min and max coordinate values
     *------------------------------*/
    for(i=1; i<nn; ++i) {
      if(crds[i].xyz[0] < mm[0]) {
        mm[0] = crds[i].xyz[0];
      }
      if(crds[i].xyz[0] > mm[1]) {
        mm[1] = crds[i].xyz[0];
      }
      if(crds[i].xyz[1] < mm[2]) {
        mm[2] = crds[i].xyz[1];
      }
      if(crds[i].xyz[1] > mm[3]) {
        mm[3] = crds[i].xyz[1];
      }
      if(crds[i].xyz[2] < mm[4]) {
        mm[4] = crds[i].xyz[2];
      }
      if(crds[i].xyz[2] > mm[5]) {
        mm[5] = crds[i].xyz[2];
      }
    }

    fprintf(stderr,"   Global coordinate ranges:\n");
    fprintf(stderr,"      min x: %g\n",mm[0]);
    fprintf(stderr,"      max x: %g\n",mm[1]);
    fprintf(stderr,"      min y: %g\n",mm[2]);
    fprintf(stderr,"      max y: %g\n",mm[3]);
    fprintf(stderr,"      min z: %g\n",mm[4]);
    fprintf(stderr,"      max z: %g\n",mm[5]);


    /* Free the allocated memory
     *--------------------------*/
    free(crds);
    if(Node_labels) {
      free(nodeids);
    }
  }



  for(p=0; p<Num_parts; ++p) {
    pn = p+1;

    fprintf(stderr,"\n");
    fprintf(stderr," Part %d:\n",pn);

    /*-----------------------
     * For unstructured parts
     *-----------------------*/
    if(Pbuild[p].type == Z_UNSTRUCTURED) {

      for(et=0; et<Z_MAXTYPE; ++et) {
        ne = Pbuild[p].ne[et];

        if(ne > 0) {

          conns[et] = (int **) calloc(ne,sizeof(int *));
          if(conns[et] == (int **) NULL) {
            fprintf(stderr,"Error: allocating conns array\n");
            return(Z_ERR);
          }
          else {
            for(i=0; i<ne; ++i) {
              conns[et][i] = (int *) calloc(Elem_info[et].con_len,sizeof(int));
              if(conns[et][i] == (int *) NULL) {
                fprintf(stderr,"Error: allocating conns array\n");
                return(Z_ERR);
              }
            }
          }

          if(Element_labels) {
            elemids[et] = (int *) calloc(ne,sizeof(int));
            if(elemids[et] == (int *) NULL) {
              fprintf(stderr,"Error: allocating elemids array\n");
              return(Z_ERR);
            }
          }
        }
      }

      /* Get the elements
       *-----------------*/
      err = USERD_get_element_connectivities_for_part(pn,conns);
      if(err == Z_ERR) {
        fprintf(stderr,"Error: getting element connectivities\n");
        return(Z_ERR);
      }

      /* Get the element ids - if any
       *-----------------------------*/
      if(Element_labels) {
        err = USERD_get_element_ids_for_part(pn,elemids);
        if(err == Z_ERR) {
          fprintf(stderr,"Error: getting element ids\n");
          return(Z_ERR);
        }
      }

      /* Echo "some" info
       *-----------------*/
      for(et=0; et<Z_MAXTYPE; ++et) {
        ne = Pbuild[p].ne[et];

        if(ne > 0) {

          /* First element of the type
           *--------------------------*/
          i = 0;
          fprintf(stderr,"   %s Element %d of %d:\n",Elem_info[et].name,i+1,ne);
          if(Element_labels) {
            fprintf(stderr,"      id: %d\n",elemids[et][i]);
          }
          fprintf(stderr,"      connectivity:");
          for(j=0; j<Elem_info[et].con_len; ++j) {
            fprintf(stderr," %d",conns[et][i][j]);
          }
          fprintf(stderr,"\n");

          /* Last element of the type
           *-------------------------*/
          i = ne - 1;
          if(i > 0) {
            fprintf(stderr,"   %s Element %d of %d:\n",Elem_info[et].name,i+1,ne);
            if(Element_labels) {
              fprintf(stderr,"      id: %d\n",elemids[et][i]);
            }
            fprintf(stderr,"      connectivity:");
            for(j=0; j<Elem_info[et].con_len; ++j) {
              fprintf(stderr," %d",conns[et][i][j]);
            }
            fprintf(stderr,"\n");
          }
        }
      }

      /* Free the allocated memory
       *--------------------------*/
      for(et=0; et<Z_MAXTYPE; ++et) {
        ne = Pbuild[p].ne[et];

        if(ne > 0) {
          for(i=0; i<ne; ++i) {
            free(conns[et][i]);
          }
          free(conns[et]);

          if(Element_labels) {
            free(elemids[et]);
          }
        }
      }
    }


    /*---------------------
     * For structured parts
     *---------------------*/
    else {

      /* Get the block coords
       *---------------------*/
      for(comp=0; comp<3; ++comp) {
        if(Pbuild[p].ne[comp] < 1) {
          bdim[comp] = 1;
        }
        else {
          bdim[comp] = Pbuild[p].ne[comp];
        }
      }
      nn = bdim[0] * bdim[1] * bdim[2];

      bd1 = bdim[0]-1;
      if(bd1 < 1) {
        bd1 = 1;
      }
      bd2 = bdim[1]-1;
      if(bd2 < 1) {
        bd2 = 1;
      }
      bd3 = bdim[2]-1;
      if(bd3 < 1) {
        bd3 = 1;
      }
      ne = bd1 * bd2 * bd3;


      /* Determine cell type
       *--------------------*/
      num_dims = 3;
      for(i=0; i<3; ++i) {
        if(bdim[i] == 1) {
          --num_dims;
        }
      }
      if(num_dims == 3) {
        cell_type = Z_HEX08;
      }
      else if(num_dims == 2) {
        cell_type = Z_QUA04;
      }
      else {
        cell_type = Z_BAR02;
      }

      coords = (float **) calloc(num_dims,sizeof(float *));
      if(coords == (float **) NULL) {
        fprintf(stderr,"Error: allocating coords array\n");
        return(Z_ERR);
      }
      else {
        for(i=0; i<num_dims; ++i) {
          coords[i] = (float *) calloc(nn,sizeof(float));
          if(coords[i] == (float *) NULL) {
            fprintf(stderr,"Error: allocating coords array\n");
            return(Z_ERR);
          }
        }
      }

      /* Get the coords
       *---------------*/
      for(comp=0; comp<num_dims; ++comp) {
        err = USERD_get_block_coords_by_component(pn,comp,coords[comp]);
        if(err == Z_ERR) {
          fprintf(stderr,"Error: getting block coords\n");
          return(Z_ERR);
        }
      }


      /* Echo "some" info
       *-----------------*/

      /* First node
       *-----------*/
      if(nn > 0) {
        i = 0;
        fprintf(stderr,"   Node %d of %d:\n",i+1,nn);

        if(num_dims == 3) {
          fprintf(stderr,"      x y z coordinates: %g %g %g\n",
                  coords[0][i], coords[1][i], coords[2][i]);
          mm[0] = mm[1] = coords[0][i];
          mm[2] = mm[3] = coords[1][i];
          mm[4] = mm[5] = coords[2][i];
        }
        else if(num_dims == 2) {
          fprintf(stderr,"      x y coordinates: %g %g\n",
                  coords[0][i], coords[1][i]);
          mm[0] = mm[1] = coords[0][i];
          mm[2] = mm[3] = coords[1][i];
        }
        else {
          fprintf(stderr,"      x coordinates: %g\n",
                  coords[0][i]);
          mm[0] = mm[1] = coords[0][i];
        }


        /* Last node
         *----------*/
        i = nn-1;
        if(i > 1) {
          fprintf(stderr,"   Node %d of %d:\n",i+1,nn);

          if(num_dims == 3) {
            fprintf(stderr,"      x y z coordinates: %g %g %g\n",
                    coords[0][i], coords[1][i], coords[2][i]);
          }
          else if(num_dims == 2) {
            fprintf(stderr,"      x y coordinates: %g %g\n",
                    coords[0][i], coords[1][i]);
          }
          else {
            fprintf(stderr,"      x coordinates: %g\n",
                    coords[0][i]);
          }
        }
      }

      /* Min and max coordinate values
       *------------------------------*/
      for(i=2; i<=nn; ++i) {
        if(coords[0][i] < mm[0]) {
          mm[0] = coords[0][i];
        }
        if(coords[0][i] > mm[1]) {
          mm[1] = coords[0][i];
        }
        if(num_dims > 1) {
          if(coords[1][i] < mm[2]) {
            mm[2] = coords[1][i];
          }
          if(coords[1][i] > mm[3]) {
            mm[3] = coords[1][i];
          }
        }
        if(num_dims > 2) {
          if(coords[2][i] < mm[4]) {
            mm[4] = coords[2][i];
          }
          if(coords[2][i] > mm[5]) {
            mm[5] = coords[2][i];
          }
        }
      }

      fprintf(stderr,"   Coordinate ranges:\n");
      fprintf(stderr,"      min x: %g\n",mm[0]);
      fprintf(stderr,"      max x: %g\n",mm[1]);
      if(num_dims > 1) {
        fprintf(stderr,"      min y: %g\n",mm[2]);
        fprintf(stderr,"      max y: %g\n",mm[3]);
      }
      if(num_dims > 2) {
        fprintf(stderr,"      min z: %g\n",mm[4]);
        fprintf(stderr,"      max z: %g\n",mm[5]);
      }

      /* Free the allocated memory - so far
       *-----------------------------------*/
      for(i=0; i<num_dims; ++i) {
        free(coords[i]);
      }
      free(coords);


      /* Get the block iblanking - if any
       *---------------------------------*/
      if(Pbuild[p].type == Z_IBLANKED) {

        iblanking = (int *) calloc(nn,sizeof(int));
        if(iblanking == (int *) NULL) {
          fprintf(stderr,"Error: allocating iblanking array\n");
          return(Z_ERR);
        }

        err = USERD_get_block_iblanking(pn,iblanking);
        if(err == Z_ERR) {
          fprintf(stderr,"Error: getting block iblanking\n");
          return(Z_ERR);
        }

        /* Echo "some" info
         *-----------------*/
        ib[Z_EXT]    = 0;
        ib[Z_INT]    = 0;
        ib[Z_BND]    = 0;
        ib[Z_INTBND] = 0;
        ib[Z_SYM]    = 0;

        for(i=0; i<nn; ++i) {
          ++ib[iblanking[i]];
        }

        fprintf(stderr,"   Iblanking breakdown:\n");
        fprintf(stderr,"      Number of Z_EXT:    %d\n",ib[Z_EXT]);
        fprintf(stderr,"      Number of Z_INT:    %d\n",ib[Z_INT]);
        fprintf(stderr,"      Number of Z_BND:    %d\n",ib[Z_BND]);
        fprintf(stderr,"      Number of Z_INTBND: %d\n",ib[Z_INTBND]);
        fprintf(stderr,"      Number of Z_SYM:    %d\n",ib[Z_SYM]);

        free(iblanking);
      }
    }
  }

  return(Z_OK);
}


/*-----------
 * var_loader
 *-----------*/
static int
var_loader(int var_time_step)
{
  int i, j, k;
  int err;
  int v, vn;
  int var_timeset;
  int p, pn;
  int et, e1, e2;
  int num_comps;
  int num_dims;
  int nsize;
  int comp;
  int bdim[3];
  int ne;
  int cell_type;
  float constant_val;
  char line1[Z_BUFL];
  char line2[Z_BUFL];
  float *values;
  float *tvalues;
  float minv[3],maxv[3];
  int bd1,bd2,bd3;


  fprintf(stderr,"\n--------------- var_loader --------------\n");

  if(Num_time_steps > 1 && v == 0) {
    if(var_time_step > (Num_time_steps - 1)) {
      var_time_step = Num_time_steps - 1;
    }

    /* Set the time step - to first step by default, but
     * can set it at others using -vts command argument
     *---------------------------------------------------------*/
    USERD_set_time_step(var_time_step);

    fprintf(stderr," Using time step: %d  (where range is %d through %d)\n\n",
            var_time_step,0,Num_time_steps-1);
  }

  for(v=0; v<Num_vars; ++v) {
    vn = v + 1;

    if(v > 0) {
      fprintf(stderr,"\n");
    }
    if(Varinfo[v].classify == Z_PER_NODE) {
      fprintf(stderr," Z_PER_NODE Variable %d:\n",vn);
    }
    else {
      fprintf(stderr," Z_PER_ELEM Variable %d:\n",vn);
    }

    /* Constants
     *----------*/
    if(Varinfo[v].type == Z_CONSTANT) {

      constant_val = USERD_get_constant_value(vn);
      fprintf(stderr,"   Constant (%s):\n",Varinfo[v].description);
      fprintf(stderr,"     value: %g\n",constant_val);
    }


    /* Scalars, Vectors, Tensors
     *--------------------------*/
    else {

      /* Get the var description line
       *-----------------------------*/
      err = USERD_get_description_lines(Z_VARI,vn,line1,line2);
      if(err == Z_OK) {
        fprintf(stderr,"   Desc line: %s\n",line1);
      }
      else {
        fprintf(stderr,"Error: getting var description line\n");
        return(Z_ERR);
      }


      /* Get the values by component
       *-----------------------------*/
      if(Varinfo[v].type == Z_SCALAR) {
        num_comps = 1;
      }
      else if(Varinfo[v].type == Z_VECTOR) {
        num_comps = 3;
      }
      else if(Varinfo[v].type == Z_TENSOR) {
        num_comps = 6;
      }
      else if(Varinfo[v].type == Z_TENSOR9) {
        num_comps = 9;
      }


      /* Per_Node
       *---------*/
      if(Varinfo[v].classify == Z_PER_NODE) {

        for(p=0; p<Num_parts; ++p) {
          pn = p + 1;

          if(Pbuild[p].type == Z_UNSTRUCTURED) {
            nsize = Pbuild[p].nn;
          }
          else {
            for(comp=0; comp<3; ++comp) {
              if(Pbuild[p].ne[comp] < 1) {
                bdim[comp] = 1;
              }
              else {
                bdim[comp] = Pbuild[p].ne[comp];
              }
            }
            nsize = bdim[0] * bdim[1] * bdim[2];
          }


          fprintf(stderr,"   For part %d, with %d nodes:\n",pn,nsize);

          if(nsize > 0) {
            values = (float *) calloc((num_comps * nsize),sizeof(float));
            if(values == (float *) NULL) {
              fprintf(stderr,"Error: allocating variable values\n");
              return(Z_ERR);
            }

            if(num_comps == 1) {

              if(Pbuild[p].type == Z_UNSTRUCTURED) {
                err = USERD_get_scalar_values(vn,
                                              pn,
                                              0,
                                              values);
                if(err == Z_ERR) {
                  fprintf(stderr,"Error: getting scalar values\n");
                  return(Z_ERR);
                }
              }
              else {
                err = USERD_get_block_scalar_values(pn,
                                                    vn,
                                                    values);
                if(err == Z_ERR) {
                  fprintf(stderr,"Error: getting block scalar values\n");
                  return(Z_ERR);
                }
              }

              /* For the component, show 1st node, last node, min, max values
               *-------------------------------------------------------------*/
              minv[0] = maxv[0] = values[0];
              for(i=0; i<nsize; ++i) {
                if(values[i] < minv[0]) {
                  minv[0] = values[i];
                }
                if(values[i] > maxv[0]) {
                  maxv[0] = values[i];
                }
              }

              fprintf(stderr,"       node %10d value: %g\n",1,values[0]);
              fprintf(stderr,"       node %10d value: %g\n",nsize,values[nsize-1]);
              fprintf(stderr,"       min value:             %g\n",minv[0]);
              fprintf(stderr,"       max value:             %g\n",maxv[0]);

            }

            else if(num_comps == 3) {

              if(Pbuild[p].type == Z_UNSTRUCTURED) {
                err = USERD_get_vector_values(vn,
                                              pn,
                                              0,
                                              values);
                if(err == Z_ERR) {
                  fprintf(stderr,"Error: getting vector values\n");
                  return(Z_ERR);
                }
              }
              else {

                tvalues = (float *) calloc(nsize,sizeof(float));
                if(tvalues == (float *) NULL) {
                  fprintf(stderr,"Error: allocating tvalues array\n");
                  return(Z_ERR);
                }

                for(i=0; i<3; ++i) {
                  err = USERD_get_block_vector_values_by_component(pn,
                                                                   vn,
                                                                   i,
                                                                   tvalues);
                  if(err == Z_ERR) {
                    fprintf(stderr,"Error: getting vector values\n");
                    return(Z_ERR);
                  }
                  for(j=0; j<nsize; ++j) {
                    k = j*3 + i;
                    values[k] = tvalues[j];
                  }
                }
                free(tvalues);
              }

              /* For the component, show 1st node, last node, min, max values
               *-------------------------------------------------------------*/
              minv[0] = maxv[0] = values[0];
              minv[1] = maxv[1] = values[1];
              minv[2] = maxv[2] = values[2];
              for(i=0; i<nsize; ++i) {
                j = i*3;
                for(k=0; k<3; ++k) {
                  if(values[j+k] < minv[k]) {
                    minv[k] = values[j+k];
                  }
                  if(values[j+k] > maxv[k]) {
                    maxv[k] = values[j+k];
                  }
                }
              }

              fprintf(stderr,"       node %10d values: %g %g %g\n",1,
                      values[0],values[1],values[2]);
              fprintf(stderr,"       node %10d values: %g %g %g\n",nsize,
                      values[3*nsize-3],values[3*nsize-2],values[3*nsize-1]);
              fprintf(stderr,"       min values:             %g %g %g\n",
                      minv[0],minv[1],minv[2]);
              fprintf(stderr,"       max values:             %g %g %g\n",
                      maxv[0],maxv[1],maxv[2]);

            }
            free(values);
          }
        }
      }

      /* Per_Elem
       *---------*/
      else {
        for(p=0; p<Num_parts; ++p) {
          pn = p + 1;

          if(Pbuild[p].type != Z_UNSTRUCTURED) {

            for(comp=0; comp<3; ++comp) {
              if(Pbuild[p].ne[comp] < 1) {
                bdim[comp] = 1;
              }
              else {
                bdim[comp] = Pbuild[p].ne[comp];
              }
            }
            bd1 = bdim[0]-1;
            if(bd1 < 1) {
              bd1 = 1;
            }
            bd2 = bdim[1]-1;
            if(bd2 < 1) {
              bd2 = 1;
            }
            bd3 = bdim[2]-1;
            if(bd3 < 1) {
              bd3 = 1;
            }
            nsize = bd1 * bd2 * bd3;


            /* Determine cell type
             *--------------------*/
            num_dims = 3;
            for(i=0; i<3; ++i) {
              if(bdim[i] == 1) {
                --num_dims;
              }
            }
            if(num_dims == 3) {
              cell_type = Z_HEX08;
            }
            else if(num_dims == 2) {
              cell_type = Z_QUA04;
            }
            else {
              cell_type = Z_BAR02;
            }
          }

          if(Pbuild[p].type == Z_UNSTRUCTURED) {
            e1 = 0;
            e2 = Z_MAXTYPE-1;
          }
          else {
            e1 = e2 = cell_type;
          }

          for(et=e1; et<=e2; ++et) {

            if(Pbuild[p].type == Z_UNSTRUCTURED) {
              nsize = Pbuild[p].ne[et];
            }

            if(nsize > 0) {

              fprintf(stderr,"   For part %d, with %d elems of type %s:\n",
                      pn,nsize,Elem_info[et].name);

              values = (float *) calloc((num_comps * nsize),sizeof(float));
              if(values == (float *) NULL) {
                fprintf(stderr,"Error: allocating variable values\n");
                return(Z_ERR);
              }

              if(num_comps == 1) {
                if(Pbuild[p].type == Z_UNSTRUCTURED) {
                  err = USERD_get_scalar_values(vn,
                                                pn,
                                                et,
                                                values);
                  if(err == Z_ERR) {
                    fprintf(stderr,"Error: getting scalar values\n");
                    return(Z_ERR);
                  }
                }
                else {
                  err = USERD_get_block_scalar_values(pn,
                                                      vn,
                                                      values);
                  if(err == Z_ERR) {
                    fprintf(stderr,"Error: getting block scalar values\n");
                    return(Z_ERR);
                  }
                }

                /* For the component, show 1st node, last node, min, max values
                 *-------------------------------------------------------------*/
                minv[0] = maxv[0] = values[0];
                for(i=1; i<nsize; ++i) {
                  if(values[i] < minv[0]) {
                    minv[0] = values[i];
                  }
                  if(values[i] > maxv[0]) {
                    maxv[0] = values[i];
                  }
                }

                fprintf(stderr,"       elem %10d value: %g\n",1,values[0]);
                fprintf(stderr,"       elem %10d value: %g\n",nsize,values[nsize-1]);
                fprintf(stderr,"       min value:             %g\n",minv[0]);
                fprintf(stderr,"       max value:             %g\n",maxv[0]);

              }

              else if(num_comps == 3) {

                if(Pbuild[p].type == Z_UNSTRUCTURED) {
                  err = USERD_get_vector_values(vn,
                                                pn,
                                                et,
                                                values);
                  if(err == Z_ERR) {
                    fprintf(stderr,"Error: getting vector values\n");
                    return(Z_ERR);
                  }
                }
                else {

                  tvalues = (float *) calloc(nsize,sizeof(float));
                  if(tvalues == (float *) NULL) {
                    fprintf(stderr,"Error: allocating tvalues array\n");
                    return(Z_ERR);
                  }

                  for(i=0; i<3; ++i) {
                    err = USERD_get_block_vector_values_by_component(pn,
                                                                     vn,
                                                                     i,
                                                                     tvalues);
                    if(err == Z_ERR) {
                      fprintf(stderr,"Error: getting vector values\n");
                      return(Z_ERR);
                    }
                    for(j=0; j<nsize; ++j) {
                      k = j*3 + i;
                      values[k] = tvalues[j];
                    }
                  }
                  free(tvalues);
                }

                /* For the component, show 1st node, last node, min, max values
                 *-------------------------------------------------------------*/
                minv[0] = maxv[0] = values[0];
                minv[1] = maxv[1] = values[1];
                minv[2] = maxv[2] = values[2];
                for(i=1; i<=nsize; ++i) {
                  j = i*3;
                  for(k=0; k<3; ++k) {
                    if(values[j+k] < minv[k]) {
                      minv[k] = values[j+k];
                    }
                    if(values[j+k] > maxv[k]) {
                      maxv[k] = values[j+k];
                    }
                  }
                }

                fprintf(stderr,"       elem %10d values: %g %g %g\n",1,
                        values[0],values[1],values[2]);
                fprintf(stderr,"       elem %10d values: %g %g %g\n",nsize,
                        values[3*nsize-3],values[3*nsize-2],values[3*nsize-1]);
                fprintf(stderr,"       min values:             %g %g %g\n",
                        minv[0],minv[1],minv[2]);
                fprintf(stderr,"       max values:             %g %g %g\n",
                        maxv[0],maxv[1],maxv[2]);

              }
              free(values);
            }
          }
        }
      }
    }
  }

  return(Z_OK);
}

#endif


#if (defined GT_USERD_API_202)


/*---------------
 * materials_info
 *---------------*/
static int
materials_info( void )
{
  int i,j;
  int err;
  int   *num_materials;
  int   *msids;
  char **msname;
  int   *mids;
  char **mdesc;


  fprintf(stderr,"\n------------ materials_info --------------\n");

  /* Get the number of variables
   *----------------------------*/
  Num_materials_sets = USERD_get_number_of_material_sets();
  if(Num_materials_sets < 0) {
    fprintf(stderr,"Error: getting the number of material sets\n");
    return(Z_ERR);
  }
  else {
    if(Num_materials_sets == 0) {
      fprintf(stderr," No materials sets in the model\n");
      return (Z_OK);
    }
    else if(Num_materials_sets > 1) {
      fprintf(stderr," Number of materials sets: %d\n",Num_materials_sets);
      fprintf(stderr," Currently, EnSight 7.6 only supports 1 material set\n");
      return(Z_ERR);
    }
    else {
      fprintf(stderr," Number of materials sets: %d\n",Num_materials_sets);
    }
  }

  /* Get the material set index list and names
   *------------------------------------------*/
  msids = (int *) calloc(Num_materials_sets,sizeof(int));
  if(msids == (int *)NULL) {
    fprintf(stderr," Problems allocating for material set ids\n");
    return(Z_ERR);
  }

  num_materials = (int *) calloc(Num_materials_sets,sizeof(int));
  if(num_materials == (int *)NULL) {
    fprintf(stderr," Problems allocating for material set num materials\n");
    return(Z_ERR);
  }

  msname = (char **) calloc(Num_materials_sets,sizeof(char *));
  if(msname == (char **)NULL) {
    fprintf(stderr," Problems allocating for material set names\n");
    return(Z_ERR);
  }
  else {
    for(i=0; i<Num_materials_sets; ++i) {
      msname[i] = (char *) calloc(Z_BUFL,sizeof(char));
      if(msname[i] == (char *)NULL) {
        fprintf(stderr," Problems allocating for material set names\n");
        return(Z_ERR);
      }
    }
  }

  err = USERD_get_matf_set_info(msids,msname);
  if(err == Z_ERR) {
    fprintf(stderr,"Error: getting material set info\n");
  }
  else {
    for(i=0; i<Num_materials_sets; ++i) {

      /* Echo some feedback
       *-------------------*/
      fprintf(stderr," For Material set %d:\n",i+1);

      fprintf(stderr,"   id:   %d\n",msids[i]);
      fprintf(stderr,"   name: %s\n",msname[i]);

      num_materials[i] = USERD_get_number_of_materials(i);
      if(num_materials[i] < 0) {
        fprintf(stderr,"Error: getting the number of materials in set %d\n",i);
        return (Z_ERR);
      }
      else if(num_materials[i] == 0) {
        fprintf(stderr," No materials in Materials set %d\n",i);
        return (Z_OK);
      }
      else {
        mids = (int *) calloc(num_materials[i],sizeof(int));
        if(mids == (int *)NULL) {
          fprintf(stderr," Problems allocating for material ids\n");
          return(Z_ERR);
        }

        mdesc = (char **) calloc(num_materials[i],sizeof(char *));
        if(mdesc == (char **)NULL) {
          fprintf(stderr," Problems allocating for material desc\n");
          return(Z_ERR);
        }
        else {
          for(j=0; j<num_materials[i]; ++j) {
            mdesc[j] = (char *) calloc(Z_BUFL,sizeof(char));
            if(mdesc[j] == (char *)NULL) {
              fprintf(stderr," Problems allocating for material desc\n");
              return(Z_ERR);
            }
          }
        }

        err = USERD_get_matf_var_info(i,mids,mdesc);
        if(err == Z_ERR) {
          fprintf(stderr,"Error: getting material info\n");
        }
        else {

          for(j=0; j<num_materials[i]; ++j) {
            /* Echo some feedback
             *-------------------*/
            fprintf(stderr,"   For Material %d:\n",j+1);

            fprintf(stderr,"     index:       %d\n",mids[j]);
            fprintf(stderr,"     description: %s\n",mdesc[j]);
          }
        }
      }
    }
  }

  /* Free the allocated memory
   *--------------------------*/
  for(i=0; i<Num_materials_sets; ++i) {
    free(msname[i]);

    for(j=0; j<num_materials[i]; ++j) {
      free(mdesc[j]);
    }
  }
  free(msname);
  free(msids);
  free(num_materials);
  free(mdesc);
  free(mids);

  return(Z_OK);
}


/*----------------------
 * gold_materials_loader
 *----------------------*/
static int
gold_materials_loader(int geom_time_step)
{
  int i, j, k, ms, nn;
  int err, err1, err2;
  int geom_timeset;
  int p, pn;
  int et, e1, e2;
  int num_dims;
  int comp;
  int bdim[3];
  int ne;
  int cell_type;
  int bd1,bd2,bd3;
  int *ivals;
  float *fvals;
  int do_num;
  int mixed_present;
  int matf_size, matfv_size;


  fprintf(stderr,"\n-------------- materials_loader --------------\n");

  if(Num_time_sets > 0) {
    /* Get the timeset used for the geometry
     *--------------------------------------*/
    geom_timeset = USERD_get_geom_timeset_number();

    /* Get the number of time steps for this timeset
     *----------------------------------------------*/
    Num_time_steps = USERD_get_num_of_time_steps(geom_timeset);
    if(Num_time_steps < 1) {
      fprintf(stderr," Error: Num time steps returned: %d\n",Num_time_steps);
      fprintf(stderr," (Must be >0 to be okay)\n");
      return(Z_ERR);
    }
    if(geom_time_step > (Num_time_steps - 1)) {
      geom_time_step = Num_time_steps - 1;
    }

    /* Set the timeset and step - to first step by default, but
     * can set it at others using -gts command argument
     *---------------------------------------------------------*/
    USERD_set_time_set_and_step(geom_timeset,geom_time_step);

    fprintf(stderr," Using timeset:   %d  (step range is %d through %d)\n",
            geom_timeset,0,Num_time_steps-1);
    fprintf(stderr," Using time step: %d\n",geom_time_step);
  }

  for(ms=0; ms<Num_materials_sets; ++ms) {
    fprintf(stderr,"\n");
    fprintf(stderr," Materials Set %d:\n",ms+1);

    for(p=0; p<Num_parts; ++p) {
      pn = p+1;

      fprintf(stderr,"\n");
      fprintf(stderr,"   Part %d:\n",pn);

      /*-----------------------
       * For unstructured parts
       *-----------------------*/
      if(Pbuild[p].type == Z_UNSTRUCTURED) {

        e1 = 0;
        e2 = Z_MAXTYPE;
      }
      else {
        for(comp=0; comp<3; ++comp) {
          if(Pbuild[p].ne[comp] < 1) {
            bdim[comp] = 1;
          }
          else {
            bdim[comp] = Pbuild[p].ne[comp];
          }
        }
        nn = bdim[0] * bdim[1] * bdim[2];

        bd1 = bdim[0]-1;
        if(bd1 < 1) {
          bd1 = 1;
        }
        bd2 = bdim[1]-1;
        if(bd2 < 1) {
          bd2 = 1;
        }
        bd3 = bdim[2]-1;
        if(bd3 < 1) {
          bd3 = 1;
        }
        ne = bd1 * bd2 * bd3;

        /* Determine cell type
         *--------------------*/
        num_dims = 3;
        for(i=0; i<3; ++i) {
          if(bdim[i] == 1) {
            --num_dims;
          }
        }
        if(num_dims == 3) {
          cell_type = Z_HEX08;
        }
        else if(num_dims == 2) {
          cell_type = Z_QUA04;
        }
        else {
          cell_type = Z_BAR02;
        }

        e1 = cell_type;
        e2 = cell_type + 1;
      }


      for(et=e1; et<e2; ++et) {

        if(Pbuild[p].type == Z_UNSTRUCTURED) {
          ne = Pbuild[p].ne[et];
        }

        if(ne > 0) {

          /* Get the material ids, if any
           *-----------------------------*/
          err = USERD_size_matf_data(ms,
                                     pn,
                                     et,
                                     Z_MAT_INDEX,
                                     &matf_size);
          if(err == Z_OK && matf_size > 0) {


            /* Go get the material ids
             *------------------------*/
            ivals = (int *) calloc(matf_size,sizeof(int));
            if(ivals == (int *)NULL) {
              fprintf(stderr," Problems allocating for material ids\n");
              return(Z_ERR);
            }
            err = USERD_load_matf_data(ms,
                                       pn,
                                       et,
                                       Z_MAT_INDEX,
                                       ivals,
                                       fvals);
            if(err == Z_OK) {
              if(matf_size < 20) {
                fprintf(stderr,"     Printing all mat ids for %s elements\n",
                        Elem_info[et].name);
                do_num = matf_size;
              }
              else {
                fprintf(stderr,"     Printing first 20 mat ids for %s elements\n",
                        Elem_info[et].name);
                do_num = 20;
              }

              /* See if any mixed materials
               *---------------------------*/
              mixed_present = FALSE;
              for(k=0; k<matf_size; ++k) {
                if(ivals[k] < 0) {
                  mixed_present = TRUE;
                  break;
                }
              }

              /* Feedback
               *---------*/
              for(k=0; k<do_num; ++k) {
                fprintf(stderr,"       mat id[%d] = %d\n",k,ivals[k]);
              }
              free(ivals);
            }
            else {
              fprintf(stderr,"     Trouble getting mat ids for %s elements\n",
                      Elem_info[et].name);
              free(ivals);
              return(Z_ERR);
            }
          }
          else {
            fprintf(stderr,"     %s elements have no material ids\n",
                    Elem_info[et].name);
          }


          /* Get the mixed material ids, if any
           *-----------------------------------*/
          if(mixed_present) {
            err1 = USERD_size_matf_data(ms,
                                        pn,
                                        et,
                                        Z_MIX_INDEX,
                                        &matf_size);
            err2 = USERD_size_matf_data(ms,
                                        pn,
                                        et,
                                        Z_MIX_VALUE,
                                        &matfv_size);

            if(err1 == Z_OK &&
               err2 == Z_OK &&
               matf_size > 0 &&
               matfv_size > 0) {

              /* Go get the material ids
               *------------------------*/
              ivals = (int *) calloc(matf_size,sizeof(int));
              if(ivals == (int *)NULL) {
                fprintf(stderr," Problems allocating for mixed material ids\n");
                return(Z_ERR);
              }
              fvals = (float *) calloc(matfv_size,sizeof(float));
              if(fvals == (float *)NULL) {
                fprintf(stderr," Problems allocating for mixed material values\n");
                return(Z_ERR);
              }

              err1 = USERD_load_matf_data(ms,
                                          pn,
                                          et,
                                          Z_MIX_INDEX,
                                          ivals,
                                          fvals);

              err2 = USERD_load_matf_data(ms,
                                          pn,
                                          et,
                                          Z_MIX_VALUE,
                                          ivals,
                                          fvals);
              if(err1 == Z_OK &&
                 err2 == Z_OK) {
                if(matf_size < 20) {
                  fprintf(stderr,"     Printing all mixed mat ids for %s elements\n",
                          Elem_info[et].name);
                  do_num = matf_size;
                }
                else {
                  fprintf(stderr,"     Printing first 20 mixed mat ids for %s elements\n",
                          Elem_info[et].name);
                  do_num = 20;
                }
                for(k=0; k<do_num; ++k) {
                  fprintf(stderr,"       mixed mat id[%d] = %d\n",k,ivals[k]);
                }
                free(ivals);

                if(matfv_size < 20) {
                  fprintf(stderr,"     Printing all mixed mat values for %s elements\n",
                          Elem_info[et].name);
                  do_num = matfv_size;
                }
                else {
                  fprintf(stderr,"     Printing first 20 mixed mat values for %s elements\n",
                          Elem_info[et].name);
                  do_num = 20;
                }
                for(k=0; k<do_num; ++k) {
                  fprintf(stderr,"       mixed mat val[%d] = %f\n",k,fvals[k]);
                }
                free(fvals);
              }
              else {
                fprintf(stderr,"     Trouble getting mixed mat ids or vals for %s elements\n",
                        Elem_info[et].name);
                free(ivals);
                free(fvals);
                return(Z_ERR);
              }
            }
            else {
              fprintf(stderr,"     Trouble getting mixed mat sizes for %s elements\n",
                      Elem_info[et].name);
              return(Z_ERR);
            }
          }
          else {
            fprintf(stderr,"       (%s elements have no mixed material ids)\n",
                    Elem_info[et].name);
          }
        }
      }
    }
  }
  return(Z_OK);
}

#endif

/*--------------
 * entity_querys
 *--------------*/
static int
entity_querys(int var_time_step)
{
  int i, j;
  int err;
  int v, vn;
  int var_timeset;
  int p, pn;
  int et, e1, e2;
  int num_comps;
  int num_dims;
  int nsize;
  int comp;
  int bdim[3];
  int ne;
  int cell_type;
  char line1[Z_BUFL];
  char line2[Z_BUFL];
  float qvals[3];
  int bd1,bd2,bd3;


  fprintf(stderr,"\n-------------- entity_querys ------------\n");
  fprintf(stderr,"       (scalar & vector variables only)    \n");
  fprintf(stderr,"\n");

#if (defined USERD_API_100)

  if(Num_time_steps > 1) {
    /* Get the number of time steps for this timeset
     *----------------------------------------------*/
    if(var_time_step > (Num_time_steps - 1)) {
      var_time_step = Num_time_steps - 1;
    }

    /* Set the time step - to first step by default, but
     * can set it at others using -vts command argument
     *---------------------------------------------------------*/
    USERD_set_time_step(var_time_step);

    fprintf(stderr," Using time step: %d  (where range is %d through %d)\n\n",
            var_time_step,0,Num_time_steps-1);
  }
#endif

  for(v=0; v<Num_vars; ++v) {
    vn = v + 1;

    /* Scalar or vectors only
     *-----------------------*/
    if(Varinfo[v].type == Z_SCALAR || Varinfo[v].type == Z_VECTOR) {


      if(Varinfo[v].classify == Z_PER_NODE) {
        fprintf(stderr," Z_PER_NODE Variable %d:\n",vn);
      }
      else {
        fprintf(stderr," Z_PER_ELEM Variable %d:\n",vn);
      }

#if (defined GT_USERD_API_100)

      if(Num_time_sets > 0) {
        /* Get the timeset used for the variable
         *---------------------------------------*/
        var_timeset = Varinfo[v].timeset;

        /* Get the number of time steps for this timeset
         *----------------------------------------------*/
        Num_time_steps = USERD_get_num_of_time_steps(var_timeset);
        if(Num_time_steps < 1) {
          fprintf(stderr," Error: Number of time steps returned: %d\n",
                  Num_time_steps);
          fprintf(stderr," (Must be >0 to be okay)\n");
          return(Z_ERR);
        }
        if(var_time_step > (Num_time_steps - 1)) {
          var_time_step = Num_time_steps - 1;
        }

        /* Set the timeset and step - to first step by default, but
         * can set it at others using -vts command argument
         *---------------------------------------------------------*/
        USERD_set_time_set_and_step(var_timeset,var_time_step);

        fprintf(stderr,"   Using timeset:   %d  (step range is %d through %d)\n",
                var_timeset,0,Num_time_steps-1);
        fprintf(stderr,"   Using time step: %d\n",var_time_step);
      }
#endif


      /* Get the var description line
       *-----------------------------*/
#if (defined GT_USERD_API_100)
      err = USERD_get_descrip_lines(Z_VARI,vn,FALSE,line1,line2);
      if(err == Z_OK) {
        fprintf(stderr,"   Desc line: %s\n",line1);
      }
      else {
        fprintf(stderr,"Error: getting var description line\n");
        return(Z_ERR);
      }

      if(Varinfo[v].complex) {
        err = USERD_get_descrip_lines(Z_VARI,vn,TRUE,line1,line2);
        if(err == Z_OK) {
          fprintf(stderr,"   Desc line (imag): %s\n",line1);
        }
        else {
          fprintf(stderr,"Error: getting var description line (imag)\n");
          return(Z_ERR);
        }
      }
#else

      err = USERD_get_description_lines(Z_VARI,vn,line1,line2);
      if(err == Z_OK) {
        fprintf(stderr,"   Desc line: %s\n",line1);
      }
      else {
        fprintf(stderr,"Error: getting var description line\n");
        return(Z_ERR);
      }

#endif

      /* Get the values by component
       *-----------------------------*/
      if(Varinfo[v].type == Z_SCALAR) {
        num_comps = 1;
      }
      else if(Varinfo[v].type == Z_VECTOR) {
        num_comps = 3;
      }

      /* Per_Node
       *---------*/
      if(Varinfo[v].classify == Z_PER_NODE) {

        for(p=0; p<Num_parts; ++p) {
          pn = p + 1;

          if(Pbuild[p].type == Z_UNSTRUCTURED) {
            nsize = Pbuild[p].nn;
          }
          else {
            for(comp=0; comp<3; ++comp) {
              if(Pbuild[p].ne[comp] < 1) {
                bdim[comp] = 1;
              }
              else {
                bdim[comp] = Pbuild[p].ne[comp];
              }
            }
            nsize = bdim[0] * bdim[1] * bdim[2];
          }



          if(nsize > 0) {

            fprintf(stderr,"   For part %d, using node %d:\n",pn,nsize);

#if (defined GT_USERD_API_100)
            err = USERD_get_var_value_at_specific(vn,
                                                  nsize,
                                                  pn,
                                                  0,
                                                  var_time_step,
                                                  qvals,
                                                  FALSE);
#else
            err = USERD_get_variable_value_at_specific(vn,
                                                       nsize,
                                                       pn,
                                                       0,
                                                       var_time_step,
                                                       qvals);
#endif
            if(err == Z_NOT_IMPLEMENTED) {
              fprintf(stderr,"  Node and element queries not implemented\n");
              return(Z_OK);
            }
            else if(err == Z_ERR) {
              fprintf(stderr,"     Could not get value\n");
            }
            else {

              /* For the component, show 1st node, last node, min, max values
               *-------------------------------------------------------------*/
              if(Varinfo[v].type == Z_SCALAR) {
                fprintf(stderr,"     Scalar value is: %g\n",qvals[0]);
              }
              else {
                fprintf(stderr,"     Vector values are: %g %g %g\n",
                        qvals[0],qvals[1],qvals[2]);
              }

#if (defined GT_USERD_API_100)
              if(Varinfo[v].complex) {

                err = USERD_get_var_value_at_specific(vn,
                                                      nsize,
                                                      pn,
                                                      0,
                                                      var_time_step,
                                                      qvals,
                                                      TRUE);

                if(err == Z_ERR) {
                  fprintf(stderr,"     Could not get imag value\n");
                }
                else {

                  /* For the component, show 1st node, last node, min, max values
                   *-------------------------------------------------------------*/
                  if(Varinfo[v].type == Z_SCALAR) {
                    fprintf(stderr,"     Scalar value (imag) is: %g\n",qvals[0]);
                  }
                  else {
                    fprintf(stderr,"     Vector values (imag) are: %g %g %g\n",
                            qvals[0],qvals[1],qvals[2]);
                  }
                }
              }
#endif

            }
          }
        }
      }

      /* Per_Elem
       *---------*/
      else {
        for(p=0; p<Num_parts; ++p) {
          pn = p + 1;

          if(Pbuild[p].type != Z_UNSTRUCTURED) {

            for(comp=0; comp<3; ++comp) {
              if(Pbuild[p].ne[comp] < 1) {
                bdim[comp] = 1;
              }
              else {
                bdim[comp] = Pbuild[p].ne[comp];
              }
            }
            bd1 = bdim[0]-1;
            if(bd1 < 1) {
              bd1 = 1;
            }
            bd2 = bdim[1]-1;
            if(bd2 < 1) {
              bd2 = 1;
            }
            bd3 = bdim[2]-1;
            if(bd3 < 1) {
              bd3 = 1;
            }
            nsize = bd1 * bd2 * bd3;


            /* Determine cell type
             *--------------------*/
            num_dims = 3;
            for(i=0; i<3; ++i) {
              if(bdim[i] == 1) {
                --num_dims;
              }
            }
            if(num_dims == 3) {
              cell_type = Z_HEX08;
            }
            else if(num_dims == 2) {
              cell_type = Z_QUA04;
            }
            else {
              cell_type = Z_BAR02;
            }
          }

          if(Pbuild[p].type == Z_UNSTRUCTURED) {
            e1 = 0;
            e2 = Z_MAXTYPE-1;
          }
          else {
            e1 = e2 = cell_type;
          }

          for(et=e1; et<=e2; ++et) {

            if(Pbuild[p].type == Z_UNSTRUCTURED) {
              nsize = Pbuild[p].ne[et];
            }

            if(nsize > 0) {


              fprintf(stderr,"   For part %d, using elem %d of type %s:\n",
                      pn,nsize,Elem_info[et].name);


#if (defined GT_USERD_API_100)
              err = USERD_get_var_value_at_specific(vn,
                                                    nsize,
                                                    pn,
                                                    et,
                                                    var_time_step,
                                                    qvals,
                                                    FALSE);
#else
            err = USERD_get_variable_value_at_specific(vn,
                                                       nsize,
                                                       pn,
                                                       et,
                                                       var_time_step,
                                                       qvals);
#endif

              if(err == Z_NOT_IMPLEMENTED) {
                fprintf(stderr,"  Node and element queries not implemented\n");
                return(Z_OK);
              }
              else if(err == Z_ERR) {
                fprintf(stderr,"     Could not get value\n");
              }
              else {
                if(Varinfo[v].type == Z_SCALAR) {
                  fprintf(stderr,"     Scalar value is: %g\n",qvals[0]);
                }
                else {
                  fprintf(stderr,"     Vector values are: %g %g %g\n",
                          qvals[0],qvals[1],qvals[2]);
                }

#if (defined GT_USERD_API_100)
                if(Varinfo[v].complex) {

                  err = USERD_get_var_value_at_specific(vn,
                                                        nsize,
                                                        pn,
                                                        et,
                                                        var_time_step,
                                                        qvals,
                                                        TRUE);
                  if(err == Z_ERR) {
                    fprintf(stderr,"     Could not get imag value\n");
                  }
                  else {
                    if(Varinfo[v].type == Z_SCALAR) {
                      fprintf(stderr,"     Scalar value (imag) is: %g\n",qvals[0]);
                    }
                    else {
                      fprintf(stderr,"     Vector values (imag) are: %g %g %g\n",
                              qvals[0],qvals[1],qvals[2]);
                    }
                  }
                }
#endif
              }
            }
          }
        }
      }
    }
  }

  return(Z_OK);
}


/*--------------
 * exercise_bkup
 *--------------*/
static int
exercise_bkup( void )
{
  int err;
  FILE *arcfile;

  fprintf(stderr,"\n------------ exercise_archive -----------\n");

  arcfile = fopen("test.arc","wb");
  if(arcfile == (FILE *)NULL) {
    fprintf(stderr,"Error: opening test archive file\n");
    return(Z_ERR);
  }
  err = USERD_bkup(arcfile,Z_SAVE_ARCHIVE);
  if(err == Z_ERR) {
    fprintf(stderr,"Error: saving to test archive file\n");
    return(Z_ERR);
  }
  fclose(arcfile);

  arcfile = fopen("test.arc","rb");
  err = USERD_bkup(arcfile,Z_REST_ARCHIVE);
  if(err == Z_ERR) {
    fprintf(stderr,"Error: restoring from test archive file\n");
    return(Z_ERR);
  }

  fprintf(stderr," Archive test completed\n");

  fclose(arcfile);

  return(Z_OK);
}

/* -------------------------------------------------------
 *  threshold_operator1 & 2 can be one of the following
 *    Z_ELE_FAILED_NONE,           - disables checking
 *     Z_ELE_FAILED_GREATER,        - greater than
 *     Z_ELE_FAILED_LESS,           - less than
 *     Z_ELE_FAILED_EQUAL,          - equal
 *     Z_ELE_FAILED_NOT_EQUAL,      - not equal
 *     Z_ELE_FAILED_MANY            - not used
 *
 * logic_criteria2
 *      Z_ELE_FAILED_LOGIC_NONE,
 *      Z_ELE_FAILED_LOGIC_AND,
 *      Z_ELE_FAILED_LOGIC_OR,
 *      Z_ELE_FAILED_LOGIC_MANY
 *
 * ------------------------------------------------------ */
int load_fail_defaults(void)
{
  int check_for_failed = FALSE;
  int cri1 = 0;                 /* Criteria1 ELE_FAILED_GREATER, etc */
  int cri2 = 0;
  int  logic_cri2 = 0;        /* Logic for criteria 2  ELE_FAILED_LOGIC_NONE, AND, etc */
  float val1 = 0.0;           /* failure threshold 1 */
  float  val2= 0.0;           /* failure threshold 2 */
  char failed_var_name[Z_MXVARIABLEDESC]={EOS};

  check_for_failed =  USERD_get_uns_failed_params( failed_var_name,
                                                   &val1, &val2, &cri1, &cri2,
                                                   &logic_cri2 );
  if (check_for_failed == TRUE) {
    fprintf(stderr,"Failed element criteria info \n");
    fprintf(stderr,"Variable name = %s\n",failed_var_name);
    fprintf(stderr,"Criteria 1 = %d\n",cri1);
    fprintf(stderr,"Criteria 2 = %d\n",cri1);
    fprintf(stderr,"Logic criteria = %d\n",logic_cri2);
    fprintf(stderr,"Value 1 = %f\n",val1);
    fprintf(stderr,"Value 2 = %f\n",val2);
  } else {
    fprintf(stderr,"No Failed elements\n");
  }
  return(Z_OK);
}


/* End of File */

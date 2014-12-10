/*--------------------------------------------------------------*/
/* Prototype Header file for EnSight External Reader            */
/* DSO Library Routines                                         */
/*                                                              */
/* intended to be included from global_extern.h only            */
/*--------------------------------------------------------------*/
/*  *************************************************************
 *   Copyright 1998 Computational Engineering International, Inc.
 *   All Rights Reserved.
 *
 *        Restricted Rights Legend
 *
 *   Use, duplication, or disclosure of this
 *   software and its documentation by the
 *   Government is subject to restrictions as
 *   set forth in subdivision [(b)(3)(ii)] of
 *   the Rights in Technical Data and Computer
 *   Software clause at 52.227-7013.
 *  *************************************************************
 */
#ifndef GLOBAL_EXTERN_PROTO_H
#define GLOBAL_EXTERN_PROTO_H

#include <stdio.h>

#ifdef WIN32
#define W32IMPORT __declspec( dllimport )
#define W32EXPORT __declspec( dllexport )
#else
#define W32IMPORT extern
#define W32EXPORT extern
#endif

/*----------------------
 * Same in All Versions
 *----------------------*/
W32EXPORT int
USERD_get_number_of_model_parts( void );

W32EXPORT int
USERD_get_block_coords_by_component(int block_number,
                                    int which_component,
                                    float *coord_array);

W32EXPORT int
USERD_get_block_iblanking(int block_number,
                          int *iblank_array);

W32EXPORT int
USERD_get_block_scalar_values(int block_number,
                              int which_scalar,
                              float *scalar_array);

W32EXPORT int
USERD_get_block_vector_values_by_component(int block_number,
                                           int which_vector,
                                           int which_component,
                                           float *vector_array);

W32EXPORT int
USERD_get_name_of_reader(char reader_name[Z_MAX_USERD_NAME],
                         int *two_fields);

/*
 * This mechanism is used to mark the fact that a given
 * reader cannot be unloaded.  We set this by default for
 * C++ based readers as there are known issues with unloading
 * a C++ DLL on certain platforms (Linux).
 */
W32EXPORT int
USERD_reader_unloadable(void);

#ifdef __cplusplus
/*
 * Define a macro that defines the cpp function as part of the
 * USERD_get_name_of_reader declaration
 */
#ifndef NO_AUTO_UNLOADABLE_CODE

#if defined(LINUX) || defined(SGI)

#define USERD_get_name_of_reader \
  USERD_reader_unloadable(void) { return(0); } \
int USERD_get_name_of_reader

#endif

#endif

#endif

W32EXPORT int
USERD_get_reader_descrip(char descrip[Z_MAXFILENP]);

W32EXPORT int
USERD_set_filenames(char filename_1[],
                    char filename_2[],
                    char the_path[],
                    int swapbytes);

W32EXPORT int
USERD_get_number_of_files_in_dataset( void );

W32EXPORT int
USERD_get_dataset_query_file_info(Z_QFILES *qfiles);

W32EXPORT int
USERD_get_changing_geometry_status( void );

W32EXPORT int
USERD_get_node_label_status( void );

W32EXPORT int
USERD_get_element_label_status( void );

W32EXPORT int
USERD_get_number_of_variables( void );

W32EXPORT void
USERD_stop_part_building( void );

W32EXPORT int
USERD_bkup(FILE *archive_file,
           int backup_type);

/* -----------------------------------
 *   Optional routine allows getting data
 *      from the reader to modify server/client behavior
 * ------------------------------------ */
W32EXPORT int
USERD_get_extra_data(int *target,
      int *nints, int *nflts, int *nchrs,
      int *pints, float *pflts, char *pchrs);

/* ----------------------------
 *  Extra "Before" GUI stuff available for all versions of API
 *  Note: this API suite is entirely optional...
 * --------------------------- */
W32EXPORT void USERD_get_extra_gui_numbers(
         int *num_Toggles,
         int *num_pulldowns,
         int *num_fields
);

W32EXPORT int USERD_get_extra_gui_defaults(
         char **toggle_Title,             /* [num_toggles][Z_LEN_GUI_TITLE_STR] */
         int *toggle_default_status,      /* [num_toggles] */
         char **pulldown_Title,           /* [num_pulldowns][Z_LEN_GUI_TITLE_STR] */
         int *pulldown_number_in_list,    /* [num_pulldowns] */
         int *pulldown_default_selection, /* [num_pulldowns] */
         char ***pulldown_item_strings,   /* [num_pulldowns][Z_MAX_NUM_GUI_PULL_ITEMS][Z_LEN_GUI_PULL_STR] */
         char **field_Title,              /* [num_fields][Z_LEN_GUI_TITLE_STR] */
         char **field_user_string         /* [num_fields][Z_LEN_GUI_FIELD_STR] */
);

W32EXPORT void USERD_set_extra_gui_data(
                  int *toggle,            /* [num_toggle] */
                  int *pulldown,          /* [num_pulldown] */
                  char **field_text       /* [num_fields][Z_LEN_GUI_FIELD_STR] */
);

/* ----------------------------
 *  Extra "After" GUI stuff available for all versions of API
 *  Note: this API suite is entirely optional...
 * --------------------------- */
W32EXPORT void USERD_get_var_extract_gui_numbers(
         int *num_Toggles,
         int *num_pulldowns,
         int *num_fields
);

W32EXPORT int USERD_get_var_extract_gui_defaults(
         char **toggle_Title,             /* [num_toggles][Z_LEN_GUI_TITLE_STR] */
         int *toggle_default_status,      /* [num_toggles] */
         char **pulldown_Title,           /* [num_pulldowns][Z_LEN_GUI_TITLE_STR] */
         int *pulldown_number_in_list,    /* [num_pulldowns] */
         int *pulldown_default_selection, /* [num_pulldowns] */
         char ***pulldown_item_strings,   /* [num_pulldowns][Z_MAX_NUM_GUI_PULL_ITEMS][Z_LEN_GUI_PULL_STR] */
         char **field_Title,              /* [num_fields][Z_LEN_GUI_TITLE_STR] */
         char **field_user_string         /* [num_fields][Z_LEN_GUI_FIELD_STR] */
);

W32EXPORT void USERD_set_var_extract_gui_data(
                  int *toggle,            /* [num_toggle] */
                  int *pulldown,          /* [num_pulldown] */
                  char **field_text       /* [num_fields][Z_LEN_GUI_FIELD_STR] */ );

/* --------------------
 * xy-query data routines
 * -------------------- */
W32EXPORT int USERD_get_num_xy_queries(void);

W32EXPORT int USERD_get_xy_query_info(
        int query_num,
        char *query_name,
        char *query_xtitle,
        char *query_ytitle,
        int *query_num_pairs);

W32EXPORT int USERD_get_xy_query_data(
        int query_num,
        int num_vals,
        float *x_vals,
        float *y_vals);


/* This routine added so the reader can know if we are at the "right" side of
 * an interval - namely, interpolation between steps is being done in EnSight
 * (It can be in any version of EnSight)
 *----------------------------------------------------------------------------*/
W32EXPORT void
USERD_set_right_side( void );

/*---------------------------------------------
 * Routines that get the geometry in buffers,
 * used for Unstructured Auto Distribute
 * (Optional)
 *---------------------------------------------*/
W32EXPORT int
USERD_get_part_coords_in_buffers(int part_number,
                                 float **coord_array,
                                 int first,
                                 int n_beg,
                                 int n_end,
                                 int buffer_size,
                                 int *num_returned);

W32EXPORT int
USERD_get_part_node_ids_in_buffers(int part_number,
                                   int *nodeid_array,
                                   int first,
                                   int n_beg,
                                   int n_end,
                                   int buffer_size,
                                   int *num_returned);

W32EXPORT int
USERD_get_part_elements_by_type_in_buffers(int part_number,
                                           int element_type,
                                           int **conn_array,
                                           int first,
                                           int e_beg,
                                           int e_end,
                                           int buffer_size,
                                           int *num_returned);

W32EXPORT int
USERD_get_part_element_ids_by_type_in_buffers(int part_number,
                                              int element_type,
                                              int *elemid_array,
                                              int first,
                                              int e_beg,
                                              int e_end,
                                              int buffer_size,
                                              int *num_returned);
W32EXPORT int
USERD_get_var_by_component_in_buffers(int which_variable,
                                      int which_part,
                                      int var_type,
                                      int which_type,
                                      int imag_data,
                                      int component,
                                      float *var_array,
                                      int first,
                                      int ne_beg,
                                      int ne_end,
                                      int buffer_size,
                                      int leftside,
                                      int *num_returned);

W32EXPORT int
USERD_get_nsided_conn_in_buffers(int part_number,
                                 int *num_nodes_per_elem_array,
                                 int *nsided_conn_array,
                                 int first,
                                 int e_beg,
                                 int e_end,
                                 int buffer_size,
                                 int *num_returned);

W32EXPORT int
USERD_get_nfaced_conn_in_buffers(int part_number,
                                 int *nfaced_fpe_arrray,
                                 int *nfaced_npf_arrray,
                                 int *nfaced_conn_array,
                                 int first,
                                 int e_beg,
                                 int e_end,
                                 int buffer_size,
                                 int *num_returned);


/*-----------------------
 * For Version 1.000 Only
 *-----------------------*/
#if defined USERD_API_100

W32EXPORT int
USERD_get_number_of_global_nodes( void );

W32EXPORT int
USERD_get_global_coords(CRD *coord_array);

W32EXPORT int
USERD_get_global_node_ids(int *nodeid_array);

W32EXPORT int
USERD_get_element_connectivities_for_part(int part_number,
                                          int **conn_array[Z_MAXTYPE]);

W32EXPORT int
USERD_get_element_ids_for_part(int part_number,
                               int *elemid_array[Z_MAXTYPE]);

W32EXPORT int
USERD_get_vector_values(int which_vector,
                        int which_part,
                        int which_type,
                        float *vector_array);

W32EXPORT int
USERD_get_part_build_info(int *part_id,
                          int *part_types,
                          char *part_descriptions[Z_BUFL],
                          int *number_of_elements[Z_MAXTYPE],
                          int *ijk_dimensions[3],
                          int *iblanking_options[6]);

W32EXPORT int
USERD_get_scalar_values(int which_scalar,
                        int which_part,
                        int which_type,
                        float *scalar_array);

W32EXPORT int
USERD_get_variable_info(char **var_description,
                        char **var_filename,
                        int *var_type,
                        int *var_classify);

W32EXPORT int
USERD_get_description_lines(int which_type,
                            int which_var,
                            char line1[Z_BUFL],
                            char line2[Z_BUFL]);

W32EXPORT int
USERD_get_variable_value_at_specific(int which_var,
                                     int which_node_or_elem,
                                     int which_part,
                                     int which_elem_type,
                                     int time_step,
                                     float values[3]);

W32EXPORT float
USERD_get_constant_value(int which_var);

W32EXPORT int
USERD_get_solution_times(float *solution_times);
W32EXPORT void
USERD_set_time_step(int time_step);

W32EXPORT int
USERD_get_number_of_time_steps(void);

#endif


/*----------------------
 * New For Version 2.000
 *----------------------*/
#if !defined USERD_API_100

W32EXPORT int
USERD_get_part_coords(int part_number,
                      float **coord_array);

W32EXPORT int
USERD_get_part_node_ids(int part_number,
                        int *nodeid_array);

W32EXPORT int
USERD_get_part_elements_by_type(int part_number,
                                int element_type,
                                int **conn_array);
W32EXPORT int
USERD_get_part_element_ids_by_type(int part_number,
                                   int element_type,
                                   int *elemid_array);

W32EXPORT int
USERD_get_reader_version(char version_number[Z_MAX_USERD_NAME]);

W32EXPORT int
USERD_get_reader_release(char version_number[Z_MAX_USERD_NAME]);

W32EXPORT int
USERD_get_var_by_component(int which_variable,
                           int which_part,
                           int var_type,
                           int which_type,
                           int complex,
                           int component,
                           float *var_array);

W32EXPORT int
USERD_get_maxsize_info(int *max_number_of_nodes,
                       int *max_number_of_elements[Z_MAXTYPE],
                       int *max_ijk_dimensions[3]);

W32EXPORT void
USERD_exit_routine( void );

W32EXPORT int
USERD_get_gold_variable_info(char **var_description,
                             char **var_filename,
                             int *var_type,
                             int *var_classify,
                             int *var_complex,
                             char **var_ifilename,
                             float *var_freq,
                             int *var_contran,
                             int *var_timeset);
W32EXPORT int
USERD_get_model_extents( float extents[6] );

W32EXPORT int
USERD_get_descrip_lines(int which_type,
                        int which_var,
                        int imag_data,
                        char line1[Z_BUFL],
                        char line2[Z_BUFL]);

W32EXPORT int
USERD_get_var_value_at_specific(int which_var,
                                int which_node_or_elem,
                                int which_part,
                                int which_elem_type,
                                int time_step,
                                float values[3],
                                int imag_data);

W32EXPORT float
USERD_get_constant_val(int which_var, int imag_data);

W32EXPORT int
USERD_get_geom_timeset_number(void);

W32EXPORT int
USERD_get_number_of_timesets(void);

W32EXPORT int
USERD_get_timeset_description(int timeset_number,
                              char timeset_description[Z_BUFL]);

W32EXPORT int
USERD_get_sol_times(int timeset_number,
                    float *solution_times);
W32EXPORT void
USERD_set_time_set_and_step(int timeset_number,
                            int time_step);
W32EXPORT int
USERD_get_num_of_time_steps(int timeset_number);

W32EXPORT int
USERD_get_border_availability(int part_number,
                              int number_of_elements[Z_MAXTYPE]);

W32EXPORT int
USERD_get_border_elements_by_type(int part_number,
                                  int element_type,
                                  int **conn_array,
                                  short *parent_element_type,
                                  int *parent_element_num);

W32EXPORT void
USERD_set_server_number(int serv_num,
                        int tot_servs);

#endif


/*----------------------
 * New For Version 2.010
 *----------------------*/
#if defined USERD_API_201 || defined USERD_API_202 || defined USERD_API_203 || defined USERD_API_204 || defined USERD_API_205 || defined USERD_API_206 || defined USERD_API_207 || defined USERD_API_208 || defined USERD_API_209 || defined USERD_API_210
W32EXPORT int
USERD_get_ghosts_in_model_flag( void );

W32EXPORT int
USERD_get_ghosts_in_block_flag(int block_number);

W32EXPORT int
USERD_get_block_ghost_flags(int block_number,
                            int *ghost_flags);
#endif

/*--------------------------
 * Modified at Version 2.030
 *--------------------------*/
#if defined USERD_API_200 || defined USERD_API_201 || defined USERD_API_202

W32EXPORT int
USERD_get_gold_part_build_info(int *part_id,
                               int *part_types,
                               char *part_descriptions[Z_BUFL],
                               int *number_of_nodes,
                               int *number_of_elements[Z_MAXTYPE],
                               int *ijk_dimensions[3],
                               int *iblanking_options[6]);
#endif

#if defined USERD_API_203 || defined USERD_API_204 || defined USERD_API_205 || defined USERD_API_206 || defined USERD_API_207 || defined USERD_API_208 || defined USERD_API_209 || defined USERD_API_210
W32EXPORT int
USERD_get_gold_part_build_info(int *part_id,
                               int *part_types,
                               char *part_descriptions[Z_BUFL],
                               int *number_of_nodes,
                               int *number_of_elements[Z_MAXTYPE],
                               int *ijk_dimensions[9],
                               int *iblanking_options[6]);
#endif


/*----------------------
 * New For Version 2.030
 *----------------------*/
#if defined USERD_API_203 || defined USERD_API_204 || defined USERD_API_205 || defined USERD_API_206 || defined USERD_API_207 || defined USERD_API_208 || defined USERD_API_209 || defined USERD_API_210
W32EXPORT int
USERD_get_number_of_material_sets( void );

W32EXPORT int
USERD_get_matf_set_info(int *mat_set_ids,
                        char **mat_set_name);

W32EXPORT int
USERD_get_number_of_materials( int set_index );

W32EXPORT int
USERD_get_matf_var_info(int set_index,
                        int *mat_ids,
                        char **mat_desc);

W32EXPORT int
USERD_size_matf_data(int set_index,
                     int part_id,
                     int wtyp,
                     int mat_type,
                     int *matf_size );

W32EXPORT int
USERD_load_matf_data( int set_index,
                      int part_id,
                      int wtyp,
                      int mat_type,
                      int *ids_list,
                      float *val_list );

W32EXPORT int
USERD_get_nsided_conn( int part_number,
                       int *nsided_conn_array );

W32EXPORT int
USERD_get_nfaced_nodes_per_face( int part_number,
                                 int *nfaced_npf_array );

W32EXPORT int
USERD_get_nfaced_conn( int part_number,
                       int *nfaced_conn_array );

#endif

/*----------------------
 * New For Version 2.040
 *----------------------*/
#if defined USERD_API_204 || defined USERD_API_205 || defined USERD_API_206 || defined USERD_API_207 || defined USERD_API_208 || defined USERD_API_209 || defined USERD_API_210

W32EXPORT int
USERD_get_uns_failed_params(
                char *fail_var_name,           /* variable name to be used in failure
                                            must be scalar, per elem      */
                float *threshold_val1,     /* number to compare for failure */
                float *threshold_val2,     /* number to compare for failure */
                int *threshold_operator1,   /* Z_GREATER_THAN, Z_LESS_THAN,
                                            Z_EQUAL_TO */
                int *threshold_operator2,   /* Z_GREATER_THAN, Z_LESS_THAN,
                                            Z_EQUAL_TO */
                int *logic_criteria2

                );

#endif

/*----------------------
** New For Version 2.050
**----------------------*/
#if defined USERD_API_205 || defined USERD_API_206 || defined USERD_API_207 || defined USERD_API_208 || defined USERD_API_209 || defined USERD_API_210

W32EXPORT int
USERD_get_number_of_species( int set_index );

W32EXPORT int
USERD_get_matsp_info(int    set_index,
                      int   *sp_ids,
                      char **sp_desc,
                      int   *sppermatcnt,
                      int   *sppermatlis);

W32EXPORT int
USERD_rigidbody_existence( void );

#endif

/*--------------------------------------------
 * New at 2.05, but modified for Version 2.080
 *-------------------------------------------- */
#if defined USERD_API_205 || defined USERD_API_206 || defined USERD_API_207
W32EXPORT int
USERD_rigidbody_values(int part_number,
                       float values[10]);
#endif

#if defined USERD_API_208 || defined USERD_API_209 || defined USERD_API_210
W32EXPORT int
USERD_rigidbody_values(int part_number,
                       float values[14]);
#endif




/*----------------------
** New For Version 2.060
**----------------------*/
#if defined USERD_API_206 || defined USERD_API_207 || defined USERD_API_208 || defined USERD_API_209 || defined USERD_API_210

W32EXPORT int
USERD_get_structured_reader_cinching( void );

W32EXPORT int
USERD_set_block_range_and_stride(int file_pn,
                                 int mini, int maxi, int stepi,
                                 int minj, int maxj, int stepj,
                                 int mink, int maxk, int stepk);
#endif

/*----------------------
** New For Version 2.070
**----------------------*/
#if defined USERD_API_207 || defined USERD_API_208 || defined USERD_API_209 || defined USERD_API_210

/* non-optional functions go here */

#endif

/* This is optional; defaults to 'Set file' and 'Set results' if not
 * defined.  If 'two_fields' is true, then both labels must have a
 * non-NULL string otherwise the defaults will be used.
 */
W32EXPORT void
USERD_set_filename_button_labels(char filename_label_1[Z_MAX_USERD_NAME],
                                 char filename_label_2[Z_MAX_USERD_NAME]);

/* This is optional; defaults to TRUE if not defined. */
W32EXPORT int
USERD_prefer_auto_distribute(void);



/*----------------------
** New For Version 2.090
**----------------------*/
#if defined USERD_API_209 || defined USERD_API_210

/* non-optional functions go here */

#endif

/* These are optional */
W32EXPORT int
USERD_get_vglyph_counts(int *num_vglyph_vectors,
                        int *num_vglyph_timelines);

W32EXPORT int
USERD_get_vglyph_timeline_info(int vtl,
                               int *id,
                               int *numtimes,
                               int *before,
                               int *amidst,
                               int *after);

W32EXPORT int
USERD_get_vglyph_timeline_times(int vtl,
                                float *times);

W32EXPORT int
USERD_get_vglyph_vector_info(int vg,
                             int *id,
                             char *description,
                             int *type,
                             int *time_condition,
                             int *time_line,
                             int *part,
                             int *nidloc,
                             int *eidloc);

W32EXPORT int
USERD_get_vglyph_vector_values(int vg,
                               float **values);

W32EXPORT int
USERD_get_vglyph_vector_xyzloc(int vg,
                               float **xyzloc);

/*----------------------
** New For Version 2.100
**----------------------*/
#if defined USERD_API_210

W32EXPORT int
USERD_get_mat_scalars_desc(int set_index,
                           char **mesv_desc);
#endif

/* These are optional */
W32EXPORT int
USERD_get_matf_set_type(int set_index);

/* special, optional functions */
W32EXPORT void
USERD_reset_routine(void);

/*--------------------------------------------------------------------*/
#endif /*GLOBAL_EXTERN_PROTO_H*/

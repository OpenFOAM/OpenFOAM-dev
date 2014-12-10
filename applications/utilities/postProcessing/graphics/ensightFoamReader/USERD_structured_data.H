int USERD_get_block_vector_values_by_component
(
    int block_number,
    int which_vector,
    int which_component,
    float *vector_array
)
{
  return(Z_OK);
}

int USERD_get_block_coords_by_component
(
    int block_number,
    int which_component,
    float *coord_array
)
{
  return(Z_OK);
}


int USERD_get_block_iblanking
(
    int block_number,
    int *iblank_array
)
{
  return(Z_OK);
}

int USERD_get_block_scalar_values
(
    int block_number,
    int which_scalar,
    float *scalar_array
)
{
  return(Z_OK);
}
int USERD_get_ghosts_in_model_flag( void )
{
  return(Z_OK);
}

int USERD_get_ghosts_in_block_flag
(
    int block_number
)
{
  return(Z_OK);
}

int USERD_get_block_ghost_flags
(
    int block_number,
    int *ghost_flags
)
{
  return(Z_OK);
}


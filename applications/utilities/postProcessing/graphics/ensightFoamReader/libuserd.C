/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    libuserd-foam

Description
    EnSight library module to read OpenFOAM data directly without translation

    It can currently handle most cell types.

    See also: README_USERD_2.0
    24 Sep 2001: NN - Added support for Ensight API 2.0
    02 Sep 2002: NN - Added support for ghost cells
    14 Mar 2004: NN - Added patches to the parts

\*---------------------------------------------------------------------------*/

#include <stdio.h>

#include "fvCFD.H"
#include "IOobjectList.H"
#include "Cloud.H"
#include "passiveParticle.H"
#include "fvMesh.H"
#include "cellModeller.H"
#include "globalFoam.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

extern "C"
{

#include "USERD_API.H"
#include "global_extern.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// same API as in 1.0
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "USERD_bkup.H"
#include "USERD_get_name_of_reader.H"
#include "USERD_set_filenames.H"
#include "USERD_get_number_of_model_parts.H"
#include "USERD_get_changing_geometry_status.H"
#include "USERD_get_dataset_query_file_info.H"
#include "USERD_get_element_label_status.H"
#include "USERD_get_node_label_status.H"
#include "USERD_get_number_of_files_in_dataset.H"
#include "USERD_get_number_of_variables.H"
#include "USERD_stop_part_building.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// slightly changed with 2.0 from 1.0
// (to handle complex variables - not used by OpenFOAM anyway)
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "USERD_get_constant_val.H"
#include "USERD_get_descrip_lines.H"
#include "USERD_get_var_value_at_specific.H"
#include "USERD_get_gold_variable_info.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// critical changes with 2.0 from 1.0
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "USERD_get_gold_part_build_info.H"
#include "USERD_get_num_of_time_steps.H"
#include "USERD_get_sol_times.H"
#include "USERD_set_time_set_and_step.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// new additions with 2.0 from 1.0
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "USERD_get_var_by_component.H"
#include "USERD_get_part_coords.H"
#include "USERD_get_part_node_ids.H"
#include "USERD_get_part_elements_by_type.H"
#include "USERD_get_part_element_ids_by_type.H"

#include "USERD_exit_routine.H"
#include "USERD_get_model_extents.H"
#include "USERD_get_reader_version.H"
#include "USERD_get_reader_release.H"
#include "USERD_get_number_timesets.H"
#include "USERD_get_timeset_description.H"
#include "USERD_get_geom_timeset_number.H"

#include "USERD_get_border_availability.H"
#include "USERD_get_border_elements_by_type.H"

#include "USERD_get_maxsize_info.H"
#include "USERD_set_server_number.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// new additions with 2.03 from 2.02
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "USERD_get_number_of_material_sets.H"
#include "USERD_get_matf_set_info.H"
#include "USERD_get_number_of_materials.H"
#include "USERD_get_matf_var_info.H"
#include "USERD_size_matf_data.H"
#include "USERD_load_matf_data.H"
#include "USERD_get_nsided_conn.H"
#include "USERD_get_nfaced_nodes_per_face.H"
#include "USERD_get_nfaced_conn.H"

//**********************************************************************
//======================================================================
// STRUCTURED DATA STUFF - not used in OpenFOAM
//======================================================================
//**********************************************************************

#include "USERD_structured_data.H"

}

// ************************************************************************ //

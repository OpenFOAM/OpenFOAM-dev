/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "LagrangianSubFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

#define DEFINE_LAGRANGIAN_SUB_CAPTALIZED_TYPE_FIELDS_(Type)                    \
    defineTemplate2TypeNameAndDebug(LagrangianSub##Type##Field, 0);            \
    defineTemplate2TypeNameAndDebugWithName                                    \
    (                                                                          \
        LagrangianSub##Type##SubField,                                         \
        STR(LagrangianSub##Type##Field),                                       \
        0                                                                      \
    );

#define DEFINE_LAGRANGIAN_SUB_CAPTALIZED_TYPE_FIELDS(Type)                     \
    DEFINE_LAGRANGIAN_SUB_CAPTALIZED_TYPE_FIELDS_(Type)

#define DEFINE_LAGRANGIAN_SUB_TYPE_FIELDS(Type, nullArg)                       \
    DEFINE_LAGRANGIAN_SUB_CAPTALIZED_TYPE_FIELDS(CAPITALIZE(Type))

DEFINE_LAGRANGIAN_SUB_TYPE_FIELDS(label, );
FOR_ALL_FIELD_TYPES(DEFINE_LAGRANGIAN_SUB_TYPE_FIELDS);

#undef DEFINE_LAGRANGIAN_SUB_CAPTALIZED_TYPE_FIELDS_
#undef DEFINE_LAGRANGIAN_SUB_CAPTALIZED_TYPE_FIELDS
#undef DEFINE_LAGRANGIAN_SUB_TYPE_FIELDS

}

// ************************************************************************* //

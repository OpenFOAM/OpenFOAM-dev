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

#include "LagrangianFieldSources.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

#define makeLagrangianFieldSource(LagrangianTypeFieldSource)                   \
    defineNamedTemplateTypeNameAndDebug(LagrangianTypeFieldSource, 0);         \
    template<>                                                                 \
    int LagrangianTypeFieldSource::disallowGenericLagrangianFieldSource        \
    (                                                                          \
        debug::debugSwitch("disallowGenericLagrangianFieldSource", 0)          \
    );                                                                         \
    defineTemplateRunTimeSelectionTable(LagrangianTypeFieldSource, null);      \
    defineTemplateRunTimeSelectionTable(LagrangianTypeFieldSource, dictionary)


#define makeLagrangianTemplateFieldSource(fieldType, nullArg)                  \
    makeLagrangianFieldSource                                                  \
    (                                                                          \
        CAT3(Lagrangian, CAPITALIZE(fieldType), FieldSource)                   \
    )


makeLagrangianTemplateFieldSource(label, );
FOR_ALL_FIELD_TYPES(makeLagrangianTemplateFieldSource);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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

Description
    Macros for instantiating solid chemistry models

\*---------------------------------------------------------------------------*/

#ifndef makeSolidChemistryModel_H
#define makeSolidChemistryModel_H

#include "addToRunTimeSelectionTable.H"
#include "solidChemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeSolidChemistryModel(sChemistry, SS, Comp, SThermo, GThermo)       \
                                                                              \
    typedef SS<Comp, SThermo, GThermo> SS##Comp##SThermo##GThermo;            \
                                                                              \
    typedef sChemistry<Comp, SThermo> sChemistryl##Comp##SThermo;             \
                                                                              \
    defineTemplateTypeNameAndDebugWithName                                    \
    (                                                                         \
        sChemistryl##Comp##SThermo,                                           \
        (word(sChemistry::typeName_()) + "<"#Comp"," + SThermo::typeName() +  \
        + ">").c_str(),                                                       \
        0                                                                     \
    );                                                                        \
                                                                              \
    defineTemplateTypeNameAndDebugWithName                                    \
    (                                                                         \
        SS##Comp##SThermo##GThermo,                                           \
        (word(SS::typeName_()) + "<"#Comp"," + SThermo::typeName() + "," +    \
        GThermo::typeName() + ">").c_str(),                                   \
        0                                                                     \
    );


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

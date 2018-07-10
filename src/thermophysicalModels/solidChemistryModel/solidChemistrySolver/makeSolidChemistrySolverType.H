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

Description
    Macros for instantiating solid chemistry models based on compressibility
    and transport types

\*---------------------------------------------------------------------------*/

#ifndef makeSolidChemistrySolverType_H
#define makeSolidChemistrySolverType_H

#include "addToRunTimeSelectionTable.H"

#include "noChemistrySolver.H"
#include "EulerImplicit.H"
#include "ode.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeSolidChemistrySolverType(SS, Schem, Comp, SThermo, GThermo)        \
                                                                               \
    typedef SS<Schem<Comp, SThermo, GThermo>>                                  \
        SS##Schem##Comp##SThermo##GThermo;                                     \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        SS##Schem##Comp##SThermo##GThermo,                                     \
        (#SS"<" + word(Schem<Comp, SThermo, GThermo>::typeName_())             \
      + "<"#Comp"," + SThermo::typeName()                                      \
      + ","  + GThermo::typeName() + ">>").c_str(),                            \
        0                                                                      \
    );                                                                         \
                                                                               \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        Comp,                                                                  \
        SS##Schem##Comp##SThermo##GThermo,                                     \
        thermo                                                                 \
    );


#define makeSolidChemistrySolverTypes(SolidChem, Comp, SThermo, GThermo)       \
                                                                               \
    makeSolidChemistrySolverType                                               \
    (                                                                          \
        noChemistrySolver,                                                     \
        SolidChem,                                                             \
        Comp,                                                                  \
        SThermo,                                                               \
        GThermo                                                                \
    );                                                                         \
                                                                               \
    makeSolidChemistrySolverType                                               \
    (                                                                          \
        ode,                                                                   \
        SolidChem,                                                             \
        Comp,                                                                  \
        SThermo,                                                               \
        GThermo                                                                \
    );


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

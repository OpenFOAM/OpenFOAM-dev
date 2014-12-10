/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

InClass
    Foam::fluidThermo

Description
    Macros for creating basic fluid thermo packages

\*---------------------------------------------------------------------------*/

#ifndef makeThermo_H
#define makeThermo_H

#include "fluidThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeThermoTypedefs(BaseThermo,Cthermo,Mixture,Transport,Type,Thermo,EqnOfState,Specie)\
                                                                              \
typedef                                                                       \
    Transport                                                                 \
    <                                                                         \
        species::thermo                                                       \
        <                                                                     \
            Thermo                                                            \
            <                                                                 \
                EqnOfState                                                    \
                <                                                             \
                    Specie                                                    \
                >                                                             \
            >,                                                                \
            Type                                                              \
        >                                                                     \
    > Transport##Type##Thermo##EqnOfState##Specie;                            \
                                                                              \
typedef                                                                       \
    Cthermo                                                                   \
    <                                                                         \
        BaseThermo,                                                           \
        Mixture<Transport##Type##Thermo##EqnOfState##Specie>                  \
    > Cthermo##Mixture##Transport##Type##Thermo##EqnOfState##Specie;          \
                                                                              \
defineTemplateTypeNameAndDebugWithName                                        \
(                                                                             \
    Cthermo##Mixture##Transport##Type##Thermo##EqnOfState##Specie,            \
    (                                                                         \
        #Cthermo"<"#Mixture"<"                                                \
      + Transport##Type##Thermo##EqnOfState##Specie::typeName()               \
      + ">>"                                                                  \
    ).c_str(),                                                                \
    0                                                                         \
);


#define makeThermo(BaseThermo,Cthermo,Mixture,Transport,Type,Thermo,EqnOfState,Specie)\
                                                                              \
makeThermoTypedefs                                                            \
(                                                                             \
    BaseThermo,                                                               \
    Cthermo,                                                                  \
    Mixture,                                                                  \
    Transport,                                                                \
    Type,                                                                     \
    Thermo,                                                                   \
    EqnOfState,                                                               \
    Specie                                                                    \
)                                                                             \
                                                                              \
addToRunTimeSelectionTable                                                    \
(                                                                             \
    basicThermo,                                                              \
    Cthermo##Mixture##Transport##Type##Thermo##EqnOfState##Specie,            \
    fvMesh                                                                    \
);                                                                            \
                                                                              \
addToRunTimeSelectionTable                                                    \
(                                                                             \
    fluidThermo,                                                              \
    Cthermo##Mixture##Transport##Type##Thermo##EqnOfState##Specie,            \
    fvMesh                                                                    \
);                                                                            \
                                                                              \
addToRunTimeSelectionTable                                                    \
(                                                                             \
    BaseThermo,                                                               \
    Cthermo##Mixture##Transport##Type##Thermo##EqnOfState##Specie,            \
    fvMesh                                                                    \
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

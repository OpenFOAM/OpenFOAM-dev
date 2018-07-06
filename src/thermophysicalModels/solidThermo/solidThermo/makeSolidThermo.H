/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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
    Foam::solidThermo

Description
    Macros for creating solid thermo packages

\*---------------------------------------------------------------------------*/

#ifndef makeSolidThermo_H
#define makeSolidThermo_H

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#define makeSolidThermo(BaseThermo,Cthermo,Mixture,Transport,Type,Thermo,EqnOfState,Specie)\
                                                                               \
                                                                               \
typedef                                                                        \
    Transport                                                                  \
    <                                                                          \
        species::thermo                                                        \
        <                                                                      \
            Thermo                                                             \
            <                                                                  \
                EqnOfState                                                     \
                <                                                              \
                    Specie                                                     \
                >                                                              \
            >,                                                                 \
            Type                                                               \
        >                                                                      \
    > Transport##Type##Thermo##EqnOfState##Specie;                             \
                                                                               \
typedef                                                                        \
    heThermo                                                                   \
    <                                                                          \
        BaseThermo,                                                            \
        Mixture<Transport##Type##Thermo##EqnOfState##Specie>                   \
    > heThermo##Mixture##Transport##Type##Thermo##EqnOfState##Specie;          \
                                                                               \
typedef                                                                        \
    Cthermo                                                                    \
    <                                                                          \
        BaseThermo,                                                            \
        Mixture<Transport##Type##Thermo##EqnOfState##Specie>                   \
    > Cthermo##Mixture##Transport##Type##Thermo##EqnOfState##Specie;           \
                                                                               \
                                                                               \
defineTemplateTypeNameAndDebugWithName                                         \
(                                                                              \
    Cthermo##Mixture##Transport##Type##Thermo##EqnOfState##Specie,             \
    (                                                                          \
        #Cthermo"<"#Mixture"<"                                                 \
      + Transport##Type##Thermo##EqnOfState##Specie::typeName()                \
      + ">>"                                                                   \
    ).c_str(),                                                                 \
    0                                                                          \
);                                                                             \
                                                                               \
                                                                               \
addToRunTimeSelectionTable                                                     \
(                                                                              \
    basicThermo,                                                               \
    Cthermo##Mixture##Transport##Type##Thermo##EqnOfState##Specie,             \
    fvMesh                                                                     \
);                                                                             \
                                                                               \
addToRunTimeSelectionTable                                                     \
(                                                                              \
    BaseThermo,                                                                \
    Cthermo##Mixture##Transport##Type##Thermo##EqnOfState##Specie,             \
    fvMesh                                                                     \
);                                                                             \
                                                                               \
addToRunTimeSelectionTable                                                     \
(                                                                              \
    BaseThermo,                                                                \
    Cthermo##Mixture##Transport##Type##Thermo##EqnOfState##Specie,             \
    dictionary                                                                 \
);



#define makeSolidThermoPhysicsType(BaseThermo,Cthermo,Mixture,SolidPhysicsType)\
                                                                               \
                                                                               \
                                                                               \
typedef                                                                        \
    heThermo                                                                   \
    <                                                                          \
        BaseThermo,                                                            \
        Mixture<SolidPhysicsType>                                              \
    > heThermo##Mixture##SolidPhysicsType;                                     \
                                                                               \
typedef                                                                        \
    Cthermo                                                                    \
    <                                                                          \
        BaseThermo,                                                            \
        Mixture<SolidPhysicsType>                                              \
    > Cthermo##Mixture##SolidPhysicsType;                                      \
                                                                               \
                                                                               \
defineTemplateTypeNameAndDebugWithName                                         \
(                                                                              \
    Cthermo##Mixture##SolidPhysicsType,                                        \
    (                                                                          \
        #Cthermo"<"#Mixture"<"                                                 \
      + SolidPhysicsType::typeName()                                           \
      + ">>"                                                                   \
    ).c_str(),                                                                 \
    0                                                                          \
);                                                                             \
                                                                               \
                                                                               \
addToRunTimeSelectionTable                                                     \
(                                                                              \
    basicThermo,                                                               \
    Cthermo##Mixture##SolidPhysicsType,                                        \
    fvMesh                                                                     \
);                                                                             \
                                                                               \
addToRunTimeSelectionTable                                                     \
(                                                                              \
    BaseThermo,                                                                \
    Cthermo##Mixture##SolidPhysicsType,                                        \
    fvMesh                                                                     \
);                                                                             \
                                                                               \
addToRunTimeSelectionTable                                                     \
(                                                                              \
    BaseThermo,                                                                \
    Cthermo##Mixture##SolidPhysicsType,                                        \
    dictionary                                                                 \
);
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

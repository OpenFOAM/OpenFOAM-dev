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
    Foam::solidThermo

Description
    Macros for creating reacting solid thermo packages

\*---------------------------------------------------------------------------*/

#ifndef makeReactingSolidThermo_H
#define makeReactingSolidThermo_H

#include "addToRunTimeSelectionTable.H"
#include "SpecieMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#define makeSolidReactionMixtureThermo(BaseThermo,Cthermo,Mixture,ThermoPhys) \
                                                                              \
typedef SpecieMixture                                                         \
    <                                                                         \
        Mixture                                                               \
        <                                                                     \
            ThermoPhys                                                        \
        >                                                                     \
    > SpecieMixture##Mixture##ThermoPhys;                                     \
                                                                              \
typedef                                                                       \
    heThermo                                                                  \
    <                                                                         \
        BaseThermo,                                                           \
        SpecieMixture##Mixture##ThermoPhys                                    \
    > heThermo##Mixture##ThermoPhys;                                          \
                                                                              \
typedef                                                                       \
    Cthermo                                                                   \
    <                                                                         \
        BaseThermo,                                                           \
        SpecieMixture##Mixture##ThermoPhys                                    \
    > Cthermo##Mixture##ThermoPhys;                                           \
                                                                              \
                                                                              \
defineTemplateTypeNameAndDebugWithName                                        \
(                                                                             \
    Cthermo##Mixture##ThermoPhys,                                             \
    (                                                                         \
        #Cthermo"<"#Mixture"<"                                                \
      + ThermoPhys::typeName()                                                \
      + ">>"                                                                  \
    ).c_str(),                                                                \
    0                                                                         \
);                                                                            \
                                                                              \
addToRunTimeSelectionTable                                                    \
(                                                                             \
    basicThermo,                                                              \
    Cthermo##Mixture##ThermoPhys,                                             \
    fvMesh                                                                    \
);                                                                            \
                                                                              \
addToRunTimeSelectionTable                                                    \
(                                                                             \
    BaseThermo,                                                               \
    Cthermo##Mixture##ThermoPhys,                                             \
    fvMesh                                                                    \
);                                                                            \
                                                                              \
addToRunTimeSelectionTable                                                    \
(                                                                             \
    BaseThermo,                                                               \
    Cthermo##Mixture##ThermoPhys,                                             \
    dictionary                                                                \
);


#define makeReactingSolidThermo(BaseThermo,CThermo,Mixture,Transport,Type,Thermo,EqnOfState,Specie) \
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
makeSolidReactionMixtureThermo                                                 \
(                                                                              \
    BaseThermo,                                                                \
    CThermo,                                                                   \
    Mixture,                                                                   \
    Transport##Type##Thermo##EqnOfState##Specie                                \
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

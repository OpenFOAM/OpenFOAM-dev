/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#ifndef makeReactionThermo_H
#define makeReactionThermo_H

#include "addToRunTimeSelectionTable.H"
#include "SpecieMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeReactionMixtureThermo(BaseThermo,CThermo,MixtureThermo,Mixture,ThermoPhys)   \
                                                                               \
typedef MixtureThermo                                                          \
<                                                                              \
    CThermo,                                                                   \
    SpecieMixture                                                              \
    <                                                                          \
        Mixture                                                                \
        <                                                                      \
            ThermoPhys                                                         \
        >                                                                      \
    >                                                                          \
> MixtureThermo##Mixture##ThermoPhys;                                          \
                                                                               \
defineTemplateTypeNameAndDebugWithName                                         \
(                                                                              \
    MixtureThermo##Mixture##ThermoPhys,                                        \
    (#MixtureThermo"<"#Mixture"<" + ThermoPhys::typeName() + ">>").c_str(),    \
    0                                                                          \
);                                                                             \
                                                                               \
addToRunTimeSelectionTable                                                     \
(                                                                              \
    basicThermo,                                                               \
    MixtureThermo##Mixture##ThermoPhys,                                        \
    fvMesh                                                                     \
);                                                                             \
                                                                               \
addToRunTimeSelectionTable                                                     \
(                                                                              \
    fluidThermo,                                                               \
    MixtureThermo##Mixture##ThermoPhys,                                        \
    fvMesh                                                                     \
);                                                                             \
                                                                               \
addToRunTimeSelectionTable                                                     \
(                                                                              \
    BaseThermo,                                                                \
    MixtureThermo##Mixture##ThermoPhys,                                        \
    fvMesh                                                                     \
);                                                                             \
                                                                               \
addToRunTimeSelectionTable                                                     \
(                                                                              \
    CThermo,                                                                   \
    MixtureThermo##Mixture##ThermoPhys,                                        \
    fvMesh                                                                     \
);


#define makeReactionThermo(BaseThermo,CThermo,MixtureThermo,Mixture,Transport,Type,Thermo,EqnOfState,Specie) \
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
makeReactionMixtureThermo                                                      \
(                                                                              \
    BaseThermo,                                                                \
    CThermo,                                                                   \
    MixtureThermo,                                                             \
    Mixture,                                                                   \
    Transport##Type##Thermo##EqnOfState##Specie                                \
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

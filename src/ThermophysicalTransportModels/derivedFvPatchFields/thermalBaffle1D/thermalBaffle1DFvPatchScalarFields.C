/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "thermalBaffle1DFvPatchScalarField.H"
#include "forSolids.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

typedef
    constIsoSolidTransport
    <
        species::thermo
        <
            eConstThermo
            <
                rhoConst<specie>
            >,
            sensibleInternalEnergy
        >
    > eConstSolidThermoPhysics;

typedef
    compressible::thermalBaffle1DFvPatchScalarField<eConstSolidThermoPhysics>
    thermalBaffle1DHConstSolidThermoPhysicsFvPatchScalarField;

defineTemplateTypeNameAndDebugWithName
(
    thermalBaffle1DHConstSolidThermoPhysicsFvPatchScalarField,
    "compressible::thermalBaffle1D<eConstSolidThermoPhysics>",
    0
);

addToPatchFieldRunTimeSelection
(
    fvPatchScalarField,
    thermalBaffle1DHConstSolidThermoPhysicsFvPatchScalarField
);

typedef
    exponentialSolidTransport
    <
        species::thermo
        <
            ePowerThermo
            <
                rhoConst<specie>
            >,
            sensibleInternalEnergy
        >
    > ePowerSolidThermoPhysics;

typedef
    compressible::thermalBaffle1DFvPatchScalarField<ePowerSolidThermoPhysics>
    thermalBaffle1DHPowerSolidThermoPhysicsFvPatchScalarField;

defineTemplateTypeNameAndDebugWithName
(
    thermalBaffle1DHPowerSolidThermoPhysicsFvPatchScalarField,
    "compressible::thermalBaffle1D<ePowerSolidThermoPhysics>",
    0
);

addToPatchFieldRunTimeSelection
(
    fvPatchScalarField,
    thermalBaffle1DHPowerSolidThermoPhysicsFvPatchScalarField
);

}

// ************************************************************************* //

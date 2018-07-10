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

\*---------------------------------------------------------------------------*/

#include "noThermo.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace thermalBaffleModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(noThermo, 0);

addToRunTimeSelectionTable(thermalBaffleModel, noThermo, mesh);
addToRunTimeSelectionTable(thermalBaffleModel, noThermo, dictionary);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool noThermo::read()
{
    return regionModel1D::read();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

noThermo::noThermo(const word& modelType, const fvMesh& mesh)
:
    thermalBaffleModel(mesh)
{}


noThermo::noThermo
(
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    thermalBaffleModel(modelType, mesh, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

noThermo::~noThermo()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void noThermo::preEvolveRegion()
{}


void noThermo::evolveRegion()
{}


const tmp<volScalarField> noThermo::Cp() const
{
    FatalErrorInFunction
        << "Cp field not available for " << type()
        << abort(FatalError);

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "noThermo::Cp",
                time().timeName(),
                primaryMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            primaryMesh(),
            dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0.0)
        )
    );
}

const volScalarField& noThermo::kappaRad() const
{
    FatalErrorInFunction
        << "kappa field not available for " << type()
        << abort(FatalError);
    return volScalarField::null();
}


const volScalarField& noThermo::rho() const
{
    FatalErrorInFunction
        << "rho field not available for " << type()
        << abort(FatalError);
    return volScalarField::null();
}


const volScalarField& noThermo::kappa() const
{
   FatalErrorInFunction
        << "K field not available for " << type()
        << abort(FatalError);
    return volScalarField::null();
}


const volScalarField& noThermo::T() const
{
    FatalErrorInFunction
        << "T field not available for " << type()
        << abort(FatalError);
    return volScalarField::null();
}


const solidThermo& noThermo::thermo() const
{
    FatalErrorInFunction
        << "T field not available for " << type()
        << abort(FatalError);
    return NullObjectRef<solidThermo>();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace thermalBaffleModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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

#include "isotropic.H"
#include "fvmLaplacian.H"
#include "fvcLaplacian.H"
#include "fvcSnGrad.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidThermophysicalTransportModels
{
    defineTypeNameAndDebug(isotropic, 0);
    addToRunTimeSelectionTable
    (
        solidThermophysicalTransportModel,
        isotropic,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidThermophysicalTransportModels::isotropic::isotropic
(
    const solidThermo& thermo
)
:
    solidThermophysicalTransportModel(typeName, thermo)
{
    if (!thermo.isotropic())
    {
        FatalIOErrorInFunction(*this)
            << "Cannot instantiate an isotropic transport model "
               "with anisotropic solid properties"
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::dictionary&
Foam::solidThermophysicalTransportModels::isotropic::coeffDict() const
{
    return dictionary::null;
}


bool Foam::solidThermophysicalTransportModels::isotropic::read()
{
    return true;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::solidThermophysicalTransportModels::isotropic::q() const
{
    return surfaceScalarField::New
    (
        "q",
       -fvc::interpolate(kappa())*fvc::snGrad(thermo().T())
    );
}


Foam::tmp<Foam::scalarField>
Foam::solidThermophysicalTransportModels::isotropic::qCorr
(
    const label patchi
) const
{
    return tmp<scalarField>(nullptr);
}


Foam::tmp<Foam::fvScalarMatrix>
Foam::solidThermophysicalTransportModels::isotropic::divq
(
    volScalarField& e
) const
{
    const solidThermo& thermo = this->thermo();

    // Return heat flux source as an implicit energy correction
    // to the temperature gradient flux
    return
       -fvc::laplacian(kappa(), thermo.T())
       -fvm::laplacianCorrection(kappa()/thermo.Cv(), e);
}


void Foam::solidThermophysicalTransportModels::isotropic::correct()
{
    solidThermophysicalTransportModel::correct();
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "spectEddyVisc.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(spectEddyVisc, 0);
addToRunTimeSelectionTable(LESModel, spectEddyVisc, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void spectEddyVisc::updateSubGridScaleFields(const volTensorField& gradU)
{
    const volScalarField Re(sqr(delta())*mag(symm(gradU))/nu());
    for (label i=0; i<5; i++)
    {
        nuSgs_ =
            nu()
           /(
                 scalar(1)
               - exp(-cB_*pow(nu()/(nuSgs_ + nu()), 1.0/3.0)*pow(Re, -2.0/3.0))
            );
    }

    nuSgs_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

spectEddyVisc::spectEddyVisc
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    LESModel(modelName, U, phi, transport, turbulenceModelName),
    GenEddyVisc(U, phi, transport),

    cB_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cB",
            coeffDict_,
            8.22
        )
    ),
    cK1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cK1",
            coeffDict_,
            0.83
        )
    ),
    cK2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cK2",
            coeffDict_,
            1.03
        )
    ),
    cK3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cK3",
            coeffDict_,
            4.75
        )
    ),
    cK4_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cK4",
            coeffDict_,
            2.55
        )
    )
{
    printCoeffs();

    updateSubGridScaleFields(fvc::grad(U));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volScalarField> spectEddyVisc::k() const
{
    const volScalarField eps(2*nuEff()*magSqr(symm(fvc::grad(U()))));

    return
        cK1_*pow(delta()*eps, 2.0/3.0)
       *exp(-cK2_*pow(delta(), -4.0/3.0)*nu()/pow(eps, 1.0/3.0))
      - cK3_*sqrt(eps*nu())
       *erfc(cK4_*pow(delta(), -2.0/3.0)*sqrt(nu())*pow(eps, -1.0/6.0));
}


void spectEddyVisc::correct(const tmp<volTensorField>& gradU)
{
    GenEddyVisc::correct(gradU);
    updateSubGridScaleFields(gradU());
}


bool spectEddyVisc::read()
{
    if (GenEddyVisc::read())
    {
        cB_.readIfPresent(coeffDict());
        cK1_.readIfPresent(coeffDict());
        cK2_.readIfPresent(coeffDict());
        cK3_.readIfPresent(coeffDict());
        cK4_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //

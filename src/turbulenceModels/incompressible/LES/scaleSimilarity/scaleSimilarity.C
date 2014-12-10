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

#include "scaleSimilarity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(scaleSimilarity, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

scaleSimilarity::scaleSimilarity
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    LESModel(modelName, U, phi, transport, turbulenceModelName),
    filterPtr_(LESfilter::New(U.mesh(), coeffDict())),
    filter_(filterPtr_())
{
    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volScalarField> scaleSimilarity::k() const
{
    return(0.5*(filter_(magSqr(U())) - magSqr(filter_(U()))));
}


tmp<volScalarField> scaleSimilarity::epsilon() const
{
    tmp<volSymmTensorField> D = symm(fvc::grad(U()));

    return((filter_(sqr(U())) - sqr(filter_(U()))) && D);
}


tmp<volSymmTensorField> scaleSimilarity::B() const
{
    return(filter_(sqr(U())) - sqr(filter_(U())));
}


tmp<volSymmTensorField> scaleSimilarity::devReff() const
{
    return dev(B());
}


tmp<fvVectorMatrix> scaleSimilarity::divDevReff(volVectorField& U) const
{
    return fvm::Su(fvc::div(devReff()), U);
}


tmp<fvVectorMatrix> scaleSimilarity::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    return fvm::Su(fvc::div(rho*devReff()), U);
}


void scaleSimilarity::correct(const tmp<volTensorField>&)
{}


bool scaleSimilarity::read()
{
    if (LESModel::read())
    {
        filter_.read(coeffDict());

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

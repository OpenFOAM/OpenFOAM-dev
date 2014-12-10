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

#include "GenSGSStress.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameWithName(GenSGSStress, "GenSGSStress");


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

GenSGSStress::GenSGSStress
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    LESModel(modelName, U, phi, transport, turbulenceModelName),

    ce_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ce",
            coeffDict_,
            1.048
        )
    ),

    couplingFactor_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "couplingFactor",
            coeffDict_,
            0.0
        )
    ),

    B_
    (
        IOobject
        (
            "B",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    nuSgs_
    (
        IOobject
        (
            "nuSgs",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        nu(),
        B_.boundaryField().types()
    )
{
    if (couplingFactor_.value() < 0.0 || couplingFactor_.value() > 1.0)
    {
        FatalErrorIn
        (
            "GenSGSStress::GenSGSStress"
            "(const volVectorField& U, const surfaceScalarField& phi,"
            "transportModel& transport)"
        )   << "couplingFactor = " << couplingFactor_
            << " is not in range 0 - 1" << nl
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> GenSGSStress::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            B_ - nu()*dev(twoSymm(fvc::grad(U())))
        )
    );
}


tmp<fvVectorMatrix> GenSGSStress::divDevReff
(
    volVectorField& U
) const
{
    if (couplingFactor_.value() > 0.0)
    {
        return
        (
            fvc::div(B_ + couplingFactor_*nuSgs_*fvc::grad(U))
          + fvc::laplacian
            (
                (1.0 - couplingFactor_)*nuSgs_, U, "laplacian(nuEff,U)"
            )
          - fvm::laplacian(nuEff(), U)
        );
    }
    else
    {
        return
        (
            fvc::div(B_)
          + fvc::laplacian(nuSgs_, U, "laplacian(nuEff,U)")
          - fvm::laplacian(nuEff(), U)
        );
    }
}


tmp<fvVectorMatrix> GenSGSStress::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    volScalarField muEff("muEff", rho*nuEff());

    if (couplingFactor_.value() > 0.0)
    {
        return
        (
            fvc::div(rho*B_ + couplingFactor_*rho*nuSgs_*fvc::grad(U))
          + fvc::laplacian
            (
                (1.0 - couplingFactor_)*rho*nuSgs_, U, "laplacian(muEff,U)"
            )
          - fvm::laplacian(muEff, U)
        );
    }
    else
    {
        return
        (
            fvc::div(rho*B_)
          + fvc::laplacian(rho*nuSgs_, U, "laplacian(muEff,U)")
          - fvm::laplacian(muEff, U)
        );
    }
}


bool GenSGSStress::read()
{
    if (LESModel::read())
    {
        ce_.readIfPresent(coeffDict());

        couplingFactor_.readIfPresent(coeffDict());

        if (couplingFactor_.value() < 0.0 || couplingFactor_.value() > 1.0)
        {
            FatalErrorIn("GenSGSStress::read()")
                << "couplingFactor = " << couplingFactor_
                << " is not in range 0 - 1"
                << exit(FatalError);
        }

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

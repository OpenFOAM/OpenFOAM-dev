/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2024 OpenFOAM Foundation
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

#include "sixDoFAcceleration.H"
#include "fvMatrices.H"
#include "geometricOneField.H"
#include "uniformDimensionedFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(sixDoFAcceleration, 0);
    addToRunTimeSelectionTable(fvModel, sixDoFAcceleration, dictionary);
    addBackwardCompatibleToRunTimeSelectionTable
    (
        fvModel,
        sixDoFAcceleration,
        dictionary,
        sixDoFAccelerationSource,
        "sixDoFAccelerationSource"
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::sixDoFAcceleration::readCoeffs(const dictionary& dict)
{
    UName_ = dict.lookupOrDefault<word>("U", "U");

    acceleration_.reset
    (
        Function1<vector>::New
        (
            "acceleration",
            mesh().time().userUnits(),
            dimAcceleration,
            dict
        ).ptr()
    );

    angularVelocity_.reset
    (
        Function1<vector>::New
        (
            "angularVelocity",
            mesh().time().userUnits(),
            unitRadians/dimTime,
            dict
        ).ptr()
    );

    angularAcceleration_.reset
    (
        Function1<vector>::New
        (
            "angularAcceleration",
            mesh().time().userUnits(),
            unitRadians/sqr(dimTime),
            dict
        ).ptr()
    );
}


template<class AlphaFieldType, class RhoFieldType>
void Foam::fv::sixDoFAcceleration::addForce
(
    const AlphaFieldType& alpha,
    const RhoFieldType& rho,
    const volVectorField& U,
    fvMatrix<vector>& eqn
) const
{
    const dimensionedVector a
   (
       "a",
       dimAcceleration,
       acceleration_->value(mesh().time().value())
   );

    // If gravitational force is present combine with the linear acceleration
    if (mesh().foundObject<uniformDimensionedVectorField>("g"))
    {
        uniformDimensionedVectorField& g =
            mesh().lookupObjectRef<uniformDimensionedVectorField>("g");

        const uniformDimensionedScalarField& hRef =
            mesh().lookupObject<uniformDimensionedScalarField>("hRef");

        g = g_ - a;

        dimensionedScalar ghRef(- mag(g)*hRef);

        mesh().lookupObjectRef<volScalarField>("gh") = (g & mesh().C()) - ghRef;

        mesh().lookupObjectRef<surfaceScalarField>("ghf") =
            (g & mesh().Cf()) - ghRef;
    }
    // ... otherwise include explicitly in the momentum equation
    else
    {
        eqn -= alpha*rho*a;
    }

    const dimensionedVector Omega
    (
        "Omega",
        dimless/dimTime,
        angularVelocity_->value(mesh().time().value())
    );

    const dimensionedVector dOmegaDT
    (
        "dOmegaDT",
        dimless/sqr(dimTime),
        angularAcceleration_->value(mesh().time().value())
    );

    eqn -=
    (
        alpha*rho*(2*Omega ^ U)                  // Coriolis force
      + alpha*rho*(Omega ^ (Omega ^ mesh().C())) // Centrifugal force
      + alpha*rho*(dOmegaDT ^ mesh().C())        // Angular acceleration force
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::sixDoFAcceleration::sixDoFAcceleration
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(name, modelType, mesh, dict),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    acceleration_(nullptr),
    angularVelocity_(nullptr),
    angularAcceleration_(nullptr),
    g_
    (
        mesh.foundObject<uniformDimensionedVectorField>("g")
      ? dimensionedVector(mesh.lookupObject<uniformDimensionedVectorField>("g"))
      : dimensionedVector("g", dimAcceleration, Zero)
    )
{
    readCoeffs(coeffs(dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::sixDoFAcceleration::addSupFields() const
{
    return wordList(1, UName_);
}


void Foam::fv::sixDoFAcceleration::addSup
(
    const volVectorField& U,
    fvMatrix<vector>& eqn
) const
{
    addForce(geometricOneField(), geometricOneField(), U, eqn);
}


void Foam::fv::sixDoFAcceleration::addSup
(
    const volScalarField& rho,
    const volVectorField& U,
    fvMatrix<vector>& eqn
) const
{
    addForce(geometricOneField(), rho, U, eqn);
}


void Foam::fv::sixDoFAcceleration::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    fvMatrix<vector>& eqn
) const
{
    addForce(alpha, rho, U, eqn);
}


bool Foam::fv::sixDoFAcceleration::movePoints()
{
    return true;
}


void Foam::fv::sixDoFAcceleration::topoChange(const polyTopoChangeMap&)
{}


void Foam::fv::sixDoFAcceleration::mapMesh(const polyMeshMap& map)
{}


void Foam::fv::sixDoFAcceleration::distribute(const polyDistributionMap&)
{}


bool Foam::fv::sixDoFAcceleration::read(const dictionary& dict)
{
    if (fvModel::read(dict))
    {
        readCoeffs(coeffs(dict));
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //

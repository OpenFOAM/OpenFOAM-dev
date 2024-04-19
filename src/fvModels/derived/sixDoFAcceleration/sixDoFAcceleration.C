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
#include "makeFunction1s.H"
#include "makeTableReaders.H"
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef Foam::fv::sixDoFAcceleration::accelerationVectors avType;

template<>
const char* const avType::vsType::typeName = "vectorVector";

template<>
const char* const avType::vsType::componentNames[] = {"x", "y", "z"};

template<>
const avType avType::vsType::vsType::zero(avType::uniform(vector::uniform(0)));

template<>
const avType avType::vsType::one(avType::uniform(vector::uniform(1)));

template<>
const avType avType::vsType::max(avType::uniform(vector::uniform(vGreat)));

template<>
const avType avType::vsType::min(avType::uniform(vector::uniform(-vGreat)));

template<>
const avType avType::vsType::rootMax
(
    avType::uniform(vector::uniform(rootVGreat))
);

template<>
const avType avType::vsType::rootMin
(
    avType::uniform(vector::uniform(-rootVGreat))
);

template<>
const avType avType::vsType::nan(avType::uniform(vector::uniform(NaN)));

namespace Foam
{
    makeFunction1s(avType, nullArg);
    makeFoamTableReaders(avType, nullArg);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::sixDoFAcceleration::readCoeffs()
{
    UName_ = coeffs().lookupOrDefault<word>("U", "U");

    accelerations_.reset
    (
        Function1<accelerationVectors>::New
        (
            "accelerations",
            mesh().time().userUnits(),
            unitNone,
            coeffs()
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
    const Vector<vector> accelerations
    (
        accelerations_->value(mesh().time().value())
    );

    // If gravitational force is present combine with the linear acceleration
    if (mesh().foundObject<uniformDimensionedVectorField>("g"))
    {
        uniformDimensionedVectorField& g =
            mesh().lookupObjectRef<uniformDimensionedVectorField>("g");

        const uniformDimensionedScalarField& hRef =
            mesh().lookupObject<uniformDimensionedScalarField>("hRef");

        g = g_ - dimensionedVector("a", dimAcceleration, accelerations.x());

        dimensionedScalar ghRef(- mag(g)*hRef);

        mesh().lookupObjectRef<volScalarField>("gh") = (g & mesh().C()) - ghRef;

        mesh().lookupObjectRef<surfaceScalarField>("ghf") =
            (g & mesh().Cf()) - ghRef;
    }
    // ... otherwise include explicitly in the momentum equation
    else
    {
        const dimensionedVector a("a", dimAcceleration, accelerations.x());
        eqn -= alpha*rho*a;
    }

    dimensionedVector Omega
    (
        "Omega",
        dimensionSet(0, 0, -1, 0, 0),
        accelerations.y()
    );

    dimensionedVector dOmegaDT
    (
        "dOmegaDT",
        dimensionSet(0, 0, -2, 0, 0),
        accelerations.z()
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
    UName_(coeffs().lookupOrDefault<word>("U", "U")),
    accelerations_(nullptr),
    g_
    (
        mesh.foundObject<uniformDimensionedVectorField>("g")
      ? dimensionedVector(mesh.lookupObject<uniformDimensionedVectorField>("g"))
      : dimensionedVector("g", dimAcceleration, Zero)
    )
{
    readCoeffs();
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
        readCoeffs();
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2024 OpenFOAM Foundation
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

#include "zeroDimensionalFixedPressureConstraint.H"
#include "zeroDimensionalFixedPressureModel.H"
#include "fvModels.H"
#include "fvMatrix.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(zeroDimensionalFixedPressureConstraint, 0);
    addToRunTimeSelectionTable
    (
        fvConstraint,
        zeroDimensionalFixedPressureConstraint,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::fv::zeroDimensionalFixedPressureModel&
Foam::fv::zeroDimensionalFixedPressureConstraint::model() const
{
    const fvModels& models = fvModels::New(mesh());

    forAll(models, i)
    {
        if (isA<zeroDimensionalFixedPressureModel>(models[i]))
        {
            return refCast<const zeroDimensionalFixedPressureModel>
            (
                models[i]
            );
        }
    }

    FatalErrorInFunction
        << "The " << typeName << " fvConstraint requires a corresponding "
        << zeroDimensionalFixedPressureModel::typeName << " fvModel"
        << exit(FatalError);

    return NullObjectRef<zeroDimensionalFixedPressureModel>();
}


template<class AlphaFieldType>
Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::zeroDimensionalFixedPressureConstraint::massSource
(
    const AlphaFieldType& alpha,
    const volScalarField::Internal& rho
) const
{
    // Source does not exist yet. Return zero.
    if (!sourcePtr_.valid())
    {
        return
            volScalarField::Internal::New
            (
                typedName("source"),
                mesh(),
                dimensionedScalar(dimMass/dimVolume/dimTime, 0)
            );
    }

    // Source for mass-based pressure equations
    if (sourcePtr_->dimensions() == dimMass/dimVolume/dimTime)
    {
        return alpha*sourcePtr_();
    }

    // Source for volume-based pressure equations
    if (sourcePtr_->dimensions() == dimless/dimTime)
    {
        return alpha*rho*sourcePtr_();
    }

    FatalErrorInFunction
        << "Pressure equation dimensions not recognised"
        << exit(FatalError);

    return tmp<volScalarField::Internal>(nullptr);
}


void Foam::fv::zeroDimensionalFixedPressureConstraint::readCoeffs()
{
    pName_ = coeffs().lookupOrDefault<word>("p", "p");

    rhoName_ = coeffs().lookupOrDefault<word>("rho", "rho");

    p_.reset
    (
        Function1<scalar>::New
        (
            "pressure",
            mesh().time().userUnits(),
            dimPressure,
            coeffs()
        ).ptr()
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::zeroDimensionalFixedPressureConstraint::
zeroDimensionalFixedPressureConstraint
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvConstraint(name, modelType, mesh, dict),
    pName_(word::null),
    rhoName_(word::null),
    p_(nullptr),
    sourcePtr_(nullptr)
{
    if (mesh.nGeometricD() != 0)
    {
        FatalIOErrorInFunction(dict)
            << "Zero-dimensional fvConstraint applied to a "
            << mesh.nGeometricD() << "-dimensional mesh"
            << exit(FatalIOError);
    }

    readCoeffs();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::zeroDimensionalFixedPressureConstraint::
~zeroDimensionalFixedPressureConstraint()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList
Foam::fv::zeroDimensionalFixedPressureConstraint::constrainedFields() const
{
    return wordList(1, pName_);
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::zeroDimensionalFixedPressureConstraint::pEqnSource
(
    const volScalarField& rho,
    fvMatrix<scalar>& pEqn
) const
{
    // Ensure the corresponding fvModel exits
    model();

    // Return zero if the source does not yet exist
    if (!sourcePtr_.valid())
    {
        return
            volScalarField::Internal::New
            (
                typedName("source"),
                mesh(),
                dimensionedScalar(pEqn.dimensions()/dimVolume, 0)
            );
    }

    // Return the source, multiplying by density if needed
    if (sourcePtr_->dimensions() == pEqn.dimensions()/dimVolume)
    {
        return sourcePtr_();
    }
    else if (sourcePtr_->dimensions() == pEqn.dimensions()/dimMass)
    {
        return rho()*sourcePtr_();
    }
    else
    {
        FatalErrorInFunction
            << "Dimensions of equation for pressure "
            << pEqn.psi().name() << " not recognised"
            << exit(FatalError);

        return tmp<volScalarField::Internal>(nullptr);
    }
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::zeroDimensionalFixedPressureConstraint::massSource
(
    const volScalarField::Internal& rho
) const
{
    return massSource<geometricOneField>(geometricOneField(), rho);
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::zeroDimensionalFixedPressureConstraint::massSource
(
    const volScalarField::Internal& alpha,
    const volScalarField::Internal& rho
) const
{
    return massSource<volScalarField::Internal>(alpha, rho);
}


bool Foam::fv::zeroDimensionalFixedPressureConstraint::constrain
(
    fvMatrix<scalar>& pEqn,
    const word& fieldName
) const
{
    // Create the source field if it does not already exist
    if (!sourcePtr_.valid())
    {
        sourcePtr_.set
        (
            new volScalarField::Internal
            (
                IOobject
                (
                    typedName("source"),
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh(),
                dimensionedScalar(pEqn.dimensions()/dimVolume, 0)
            )
        );
    }

    // Remove the previous iteration's source from the pressure equation
    pEqn += sourcePtr_();

    // Set the source as the residual of the pressure equation when evaluated
    // at the desired pressure
    sourcePtr_() =
        pEqn
      & volScalarField::Internal::New
        (
            "p",
            mesh(),
            dimensionedScalar
            (
                dimPressure,
                p_->value(mesh().time().value())
            )
        );

    // Add the source to the pressure equation to force the pressure towards
    // the desired value
    pEqn -= sourcePtr_();

    return true;
}


bool Foam::fv::zeroDimensionalFixedPressureConstraint::movePoints()
{
    return true;
}


void Foam::fv::zeroDimensionalFixedPressureConstraint::topoChange
(
    const polyTopoChangeMap& map
)
{}


void Foam::fv::zeroDimensionalFixedPressureConstraint::mapMesh
(
    const polyMeshMap& map
)
{}


void Foam::fv::zeroDimensionalFixedPressureConstraint::distribute
(
    const polyDistributionMap& map
)
{}


bool Foam::fv::zeroDimensionalFixedPressureConstraint::read
(
    const dictionary& dict
)
{
    if (fvConstraint::read(dict))
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

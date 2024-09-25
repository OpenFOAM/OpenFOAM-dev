/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2024 OpenFOAM Foundation
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

#include "relativeVelocityModel.H"
#include "fixedValueFvPatchFields.H"
#include "slipFvPatchFields.H"
#include "partialSlipFvPatchFields.H"
#include "fvcGrad.H"
#include "fvcDiv.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(relativeVelocityModel, 0);
    defineRunTimeSelectionTable(relativeVelocityModel, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions   * * * * * * * * * * * //

Foam::wordList Foam::relativeVelocityModel::UdmPatchFieldTypes() const
{
    const volVectorField& U = mixture_.U();

    wordList UdmTypes
    (
        U.boundaryField().size(),
        calculatedFvPatchScalarField::typeName
    );

    forAll(U.boundaryField(), i)
    {
        if
        (
            isA<fixedValueFvPatchVectorField>(U.boundaryField()[i])
         || isA<slipFvPatchVectorField>(U.boundaryField()[i])
         || isA<partialSlipFvPatchVectorField>(U.boundaryField()[i])
        )
        {
            UdmTypes[i] = fixedValueFvPatchVectorField::typeName;
        }
    }

    return UdmTypes;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativeVelocityModel::relativeVelocityModel
(
    const dictionary& dict,
    const incompressibleDriftFluxMixture& mixture,
    const uniformDimensionedVectorField& g
)
:
    mixture_(mixture),
    g_(g),
    Udm_
    (
        IOobject
        (
            "Udm",
            mixture_.time().name(),
            mixture_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mixture_.mesh(),
        dimensionedVector(dimVelocity, Zero),
        UdmPatchFieldTypes()
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::relativeVelocityModel> Foam::relativeVelocityModel::New
(
    const dictionary& dict,
    const incompressibleDriftFluxMixture& mixture,
    const uniformDimensionedVectorField& g
)
{
    const word modelType(dict.lookup(typeName));

    Info<< "Selecting relative velocity model " << modelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorInFunction(dict)
            << "Unknown time scale model type " << modelType
            << ", constructor not in hash table" << nl << nl
            << "    Valid time scale model types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << abort(FatalIOError);
    }

    return
        autoPtr<relativeVelocityModel>
        (
            cstrIter()
            (
                dict.optionalSubDict(modelType + "Coeffs"),
                mixture,
                g
            )
        );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::relativeVelocityModel::~relativeVelocityModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField>
Foam::relativeVelocityModel::acceleration() const
{
    // Dispersed phase velocity
    // const volVectorField Ud(mixture_.U() + Udm_);

    // Use the mixture rather than the dispersed-phase velocity to approximate
    // the dispersed-phase acceleration to improve stability as only the mixture
    // momentum equation is coupled to continuity and pressure
    //
    // This approximation is valid only in the limit of small drift-velocity.
    // For large drift-velocity an Euler-Euler approach should be used in
    // which both the continuous and dispersed-phase momentum equations are
    // solved and coupled to the pressure.
    const volVectorField& Ud = mixture_.U();

    return g_ - (Ud & fvc::grad(Ud));
}


Foam::tmp<Foam::volSymmTensorField> Foam::relativeVelocityModel::tauDm() const
{
    const volScalarField betac(mixture_.alphac()*mixture_.rhoc());
    const volScalarField betad(mixture_.alphad()*mixture_.rhod());

    // Calculate the relative velocity of the continuous phase w.r.t the mean
    const volVectorField Ucm(betad*Udm_/betac);

    return volSymmTensorField::New
    (
        "tauDm",
        betad*sqr(Udm_) + betac*sqr(Ucm)
    );
}


Foam::tmp<Foam::volVectorField> Foam::relativeVelocityModel::divDevTau() const
{
    return fvc::div(tauDm());
}


void Foam::relativeVelocityModel::correct()
{
    Udm_ = UdmCoeff()*acceleration();
}

// ************************************************************************* //

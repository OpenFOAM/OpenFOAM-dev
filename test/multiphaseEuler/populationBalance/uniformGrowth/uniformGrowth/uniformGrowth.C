/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "uniformGrowth.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(uniformGrowth, 0);
    addToRunTimeSelectionTable(fvModel, uniformGrowth, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::wordList Foam::fv::uniformGrowth::phaseNames() const
{
    wordList names(popBal_.uniquePhases().size());
    forAll(popBal_.uniquePhases(), uniquePhasei)
    {
        names[uniquePhasei] =
            popBal_.uniquePhases()[uniquePhasei].name();
    }
    return names;
}


Foam::wordList Foam::fv::uniformGrowth::alphaNames() const
{
    wordList names(popBal_.uniquePhases().size());
    forAll(popBal_.uniquePhases(), uniquePhasei)
    {
        names[uniquePhasei] =
            popBal_.uniquePhases()[uniquePhasei].volScalarField::name();
    }
    return names;
}


Foam::wordList Foam::fv::uniformGrowth::rhoNames() const
{
    wordList names(popBal_.uniquePhases().size());
    forAll(popBal_.uniquePhases(), uniquePhasei)
    {
        names[uniquePhasei] =
            popBal_.uniquePhases()[uniquePhasei].rho().name();
    }
    return names;
}


void Foam::fv::uniformGrowth::readCoeffs(const dictionary& dict)
{
    if (dict.lookup<word>("populationBalance") != popBal_.name())
    {
        FatalIOErrorInFunction(dict)
            << "Cannot change the population balance model of a " << type()
            << " model at run time" << exit(FatalIOError);
    }

    massFlowRate_.reset
    (
        Function1<scalar>::New
        (
            "massFlowRate",
            mesh().time().userUnits(),
            dimMass/dimTime,
            dict
        ).ptr()
    );
}


template<class Type>
void Foam::fv::uniformGrowth::addSupType
(
    const VolField<Type>& field,
    fvMatrix<Type>& eqn
) const
{
    DebugInFunction
        << "field=" << field.name()
        << ", eqnField=" << eqn.psi().name() << endl;

    FatalErrorInFunction
        << "Cannot add a mass source for field " << field.name()
        << " to equation for " << eqn.psi().name() << " because this field's "
        << "equation was not recognised as being in phase-conservative form"
        << exit(FatalError);
}


void Foam::fv::uniformGrowth::addSupType
(
    const volScalarField& rhoOrField,
    fvMatrix<scalar>& eqn
) const
{
    DebugInFunction
        << "rhoOrField=" << rhoOrField.name()
        << ", eqnField=" << eqn.psi().name() << endl;

    // Not recognised. Fail.
    addSupType<scalar>(rhoOrField, eqn);
}


template<class Type>
void Foam::fv::uniformGrowth::addSupType
(
    const volScalarField& rho,
    const VolField<Type>& field,
    fvMatrix<Type>& eqn
) const
{
    DebugInFunction
        << "rho=" << rho.name()
        << ", field=" << field.name()
        << ", eqnField=" << eqn.psi().name() << endl;

    // Multiphase mass-weighted mixture property equation
    if
    (
        rho.group() == word::null
     && rho.dimensions() == dimDensity
     && field.group() == word::null
    )
    {
        fvSpecificSource::addSupType(rho, field, eqn);
    }
    // Not recognised. Fail.
    else
    {
        addSupType<Type>(field, eqn);
    }
}


void Foam::fv::uniformGrowth::addSupType
(
    const volScalarField& alphaOrRho,
    const volScalarField& rhoOrField,
    fvMatrix<scalar>& eqn
) const
{
    DebugInFunction
        << "alphaOrRho=" << alphaOrRho.name()
        << ", rhoOrField=" << rhoOrField.name()
        << ", eqnField=" << eqn.psi().name() << endl;

    // Multiphase continuity equation
    if (findIndex(rhoNames_, rhoOrField.name()) != -1)
    {
        eqn += S(alphaOrRho.name());
    }
    // Try the general type method
    else
    {
        addSupType<scalar>(alphaOrRho, rhoOrField, eqn);
    }
}


template<class Type>
void Foam::fv::uniformGrowth::addSupType
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const VolField<Type>& field,
    fvMatrix<Type>& eqn
) const
{
    DebugInFunction
        << "alpha=" << alpha.name()
        << ", rho=" << rho.name()
        << ", field=" << field.name()
        << ", eqnField=" << eqn.psi().name() << endl;

    // Multiphase property equation
    fvSpecificSource::addSupType(alpha, rho, field, eqn);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::uniformGrowth::uniformGrowth
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvSpecificSource(name, modelType, mesh, dict),
    popBal_
    (
        mesh().lookupObject<populationBalanceModel>
        (
            coeffs(dict).lookup<word>("populationBalance")
        )
    ),
    phaseNames_(phaseNames()),
    alphaNames_(alphaNames()),
    rhoNames_(rhoNames()),
    massFlowRate_(nullptr)
{
    readCoeffs(coeffs(dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::uniformGrowth::addsSupToField(const word& fieldName) const
{
    const word fieldPhaseName = IOobject::group(fieldName);

    return
        fieldPhaseName == word::null
     || findIndex(phaseNames_, fieldPhaseName) != -1;
}


Foam::tmp<Foam::volScalarField::Internal> Foam::fv::uniformGrowth::S
(
    const word& fieldName
) const
{
    const label i = findIndex(phaseNames_, IOobject::group(fieldName));
    const phaseModel& phase = popBal_.fluid().phases()[phaseNames_[i]];

    tmp<volScalarField::Internal> tSumNPhase =
        volScalarField::Internal::New
        (
            IOobject::groupName(name() + ":sumN", phase.name()),
            mesh(),
            dimensionedScalar(inv(dimVolume), scalar(0))
        );
    tmp<volScalarField::Internal> tSumN =
        volScalarField::Internal::New
        (
            name() + ":sumN",
            mesh(),
            dimensionedScalar(inv(dimVolume), scalar(0))
        );

    forAll(popBal_.fs(), i)
    {
        tmp<volScalarField::Internal> tN =
            popBal_.phases()[i]()*popBal_.f(i)/popBal_.v(i);

        if (&popBal_.phases()[i] == &phase)
        {
            tSumNPhase.ref() += tN();
        }

        tSumN.ref() += tN;
    }

    return
        tSumNPhase
       /max(tSumN, dimensionedScalar(inv(dimVolume), rootVSmall))
       *dimensionedScalar
        (
            dimMass/dimTime,
            massFlowRate_->value(mesh().time().value())
        )
       /sum(mesh().V());
}


void Foam::fv::uniformGrowth::addSup(fvMatrix<scalar>& eqn) const
{
    DebugInFunction
        << "eqnField=" << eqn.psi().name() << endl;

    FatalErrorInFunction
        << "Field-less mass sources are not possible"
        << exit(FatalError);
}


FOR_ALL_FIELD_TYPES
(
    IMPLEMENT_FV_MODEL_ADD_FIELD_SUP,
    fv::uniformGrowth
);


FOR_ALL_FIELD_TYPES
(
    IMPLEMENT_FV_MODEL_ADD_RHO_FIELD_SUP,
    fv::uniformGrowth
);


FOR_ALL_FIELD_TYPES
(
    IMPLEMENT_FV_MODEL_ADD_ALPHA_RHO_FIELD_SUP,
    fv::uniformGrowth
);


bool Foam::fv::uniformGrowth::read(const dictionary& dict)
{
    if (fvSpecificSource::read(dict))
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

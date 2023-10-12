/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
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

#include "massTransferBase.H"
#include "fvMatrices.H"
#include "physicalProperties.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(massTransferBase, 0);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::massTransferBase::readCoeffs()
{
    phaseNames_ = coeffs().lookup<Pair<word>>("phases");

    alphaNames_ =
        coeffs().lookupOrDefault<Pair<word>>
        (
            "alphas",
            Pair<word>
            (
                IOobject::groupName("alpha", phaseNames_.first()),
                IOobject::groupName("alpha", phaseNames_.second())
            )
        );

    rhoNames_ =
        coeffs().lookupOrDefault<Pair<word>>
        (
            "rhos",
            Pair<word>
            (
                IOobject::groupName("rho", phaseNames_.first()),
                IOobject::groupName("rho", phaseNames_.second())
            )
        );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField::Internal> Foam::fv::massTransferBase::rho
(
    const label i
) const
{
    // Compressible case. Lookup the field.
    if (mesh().foundObject<volScalarField::Internal>(rhoNames()[i]))
    {
        return mesh().lookupObject<volScalarField::Internal>(rhoNames()[i]);
    }

    // Incompressible case. Read from the physical properties dictionary.
    const word physicalPropertiesName =
        IOobject::groupName(physicalProperties::typeName, phaseNames()[i]);

    if (mesh().foundObject<IOdictionary>(physicalPropertiesName))
    {
        const IOdictionary& physicalProperties =
            mesh().lookupObject<IOdictionary>(physicalPropertiesName);

        if (physicalProperties.found("rho"))
        {
            return
                volScalarField::Internal::New
                (
                    rhoNames()[i],
                    mesh(),
                    dimensionedScalar("rho", dimDensity, physicalProperties)
                );
        }
    }

    // Fail.
    FatalErrorInFunction
        << "Could not determine the density " << rhoNames()[i]
        << " for phase " << phaseNames()[i]
        << exit(FatalError);
    return tmp<volScalarField::Internal>(nullptr);
}


void Foam::fv::massTransferBase::addSupType
(
    const volScalarField& alphaOrField,
    fvMatrix<scalar>& eqn
) const
{
    DebugInFunction
        << "alphaOrField=" << alphaOrField.name()
        << ", eqnField=" << eqn.psi().name() << endl;

    const label i = index(alphaNames(), alphaOrField.name());

    // Incompressible continuity equation
    if (i != -1)
    {
        eqn += S(alphaOrField.name())/rho(i);
    }
    // Not recognised. Fail.
    else
    {
        addSupType<scalar>(alphaOrField, eqn);
    }
}


void Foam::fv::massTransferBase::addSupType
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

    // Compressible continuity equation
    if
    (
        index(alphaNames(), alphaOrRho.name()) != -1
     && index(rhoNames(), rhoOrField.name()) != -1
    )
    {
        eqn += S(alphaOrRho.name());
    }
    // Mixture property equation
    else
    {
        addSupType<scalar>(alphaOrRho, rhoOrField, eqn);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::massTransferBase::massTransferBase
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvSpecificSource(name, modelType, mesh, dict),
    phaseNames_(word::null, word::null),
    alphaNames_(word::null, word::null),
    rhoNames_(word::null, word::null)
{
    readCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::massTransferBase::addsSupToField(const word& fieldName) const
{
    const word group = IOobject::group(fieldName);

    return
        group == word::null
     || group == phaseNames_.first()
     || group == phaseNames_.second();
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::massTransferBase::S(const word& fieldName) const
{
    return sign(phaseNames(), IOobject::group(fieldName))*mDotByV();
}


void Foam::fv::massTransferBase::addSup(fvMatrix<scalar>& eqn) const
{
    DebugInFunction
        << "eqnField=" << eqn.psi().name() << endl;

    FatalErrorInFunction
        << "Field-less mass transfers are not possible"
        << exit(FatalError);
}


FOR_ALL_FIELD_TYPES(IMPLEMENT_FV_MODEL_ADD_FIELD_SUP, fv::massTransferBase);


FOR_ALL_FIELD_TYPES
(
    IMPLEMENT_FV_MODEL_ADD_RHO_FIELD_SUP,
    fv::massTransferBase
);


FOR_ALL_FIELD_TYPES
(
    IMPLEMENT_FV_MODEL_ADD_ALPHA_RHO_FIELD_SUP,
    fv::massTransferBase
);


bool Foam::fv::massTransferBase::read(const dictionary& dict)
{
    if (fvSpecificSource::read(dict))
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

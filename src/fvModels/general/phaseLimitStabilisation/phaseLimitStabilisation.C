/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2025 OpenFOAM Foundation
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

#include "phaseLimitStabilisation.H"
#include "fvMatrices.H"
#include "fvmSup.H"
#include "uniformDimensionedFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(phaseLimitStabilisation, 0);

    addToRunTimeSelectionTable
    (
        fvModel,
        phaseLimitStabilisation,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::phaseLimitStabilisation::readCoeffs(const dictionary& dict)
{
    fieldName_ = dict.lookup<word>("field");
    rateName_ = dict.lookup<word>("rate");
    residualAlpha_ = dict.lookup<scalar>("residualAlpha");
}


template<class Type, class RhoRateFieldType>
void Foam::fv::phaseLimitStabilisation::addStabilisation
(
    const volScalarField& alpha,
    const RhoRateFieldType& rhoRate,
    fvMatrix<Type>& eqn
) const
{
    eqn -= fvm::Sp(max(residualAlpha_ - alpha(), scalar(0))*rhoRate, eqn.psi());
}


template<class Type>
void Foam::fv::phaseLimitStabilisation::addSupType
(
    const volScalarField& alpha,
    const VolField<Type>& field,
    fvMatrix<Type>& eqn
) const
{
    addSupType(alpha, geometricOneField(), field, eqn);
}


template<class Type, class RhoFieldType>
void Foam::fv::phaseLimitStabilisation::addSupType
(
    const volScalarField& alpha,
    const RhoFieldType& rho,
    const VolField<Type>& field,
    fvMatrix<Type>& eqn
) const
{
    if (mesh().foundObject<volScalarField::Internal>(rateName_))
    {
        const volScalarField::Internal& rate =
            mesh().lookupObjectRef<volScalarField::Internal>(rateName_);

        addStabilisation(alpha, rho.internalField()*rate, eqn);
    }
    else if (mesh().foundObject<uniformDimensionedScalarField>(rateName_))
    {
        const uniformDimensionedScalarField& rate =
            mesh().lookupObjectRef<uniformDimensionedScalarField>(rateName_);

        addStabilisation(alpha, rho.internalField()*rate, eqn);
    }
    else
    {
        FatalErrorInFunction
            << nl
            << "    request for " << volScalarField::Internal::typeName
            << " or " << uniformDimensionedScalarField::typeName
            << " " << rateName_ << " from objectRegistry " << this->name()
            << " failed\n    available objects of type "
            << volScalarField::Internal::typeName << " are" << nl
            << mesh().toc<volScalarField::Internal>()
            << "\n    available objects of type "
            << uniformDimensionedScalarField::typeName << " are" << nl
            << mesh().toc<uniformDimensionedScalarField>()
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::phaseLimitStabilisation::phaseLimitStabilisation
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(name, modelType, mesh, dict),
    fieldName_(word::null),
    rateName_(word::null),
    residualAlpha_(NaN)
{
    readCoeffs(coeffs(dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::phaseLimitStabilisation::addSupFields() const
{
    return wordList(1, fieldName_);
}


FOR_ALL_FIELD_TYPES
(
    IMPLEMENT_FV_MODEL_ADD_RHO_FIELD_SUP,
    fv::phaseLimitStabilisation
)


FOR_ALL_FIELD_TYPES
(
    IMPLEMENT_FV_MODEL_ADD_ALPHA_RHO_FIELD_SUP,
    fv::phaseLimitStabilisation
)


bool Foam::fv::phaseLimitStabilisation::movePoints()
{
    return true;
}


void Foam::fv::phaseLimitStabilisation::topoChange(const polyTopoChangeMap&)
{}


void Foam::fv::phaseLimitStabilisation::mapMesh(const polyMeshMap& map)
{}


void Foam::fv::phaseLimitStabilisation::distribute
(
    const polyDistributionMap&
)
{}


bool Foam::fv::phaseLimitStabilisation::read(const dictionary& dict)
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

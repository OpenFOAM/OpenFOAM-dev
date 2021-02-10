/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2021 OpenFOAM Foundation
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

#include "phaseLimitStabilization.H"
#include "fvMatrices.H"
#include "fvmSup.H"
#include "uniformDimensionedFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(phaseLimitStabilization, 0);

    addToRunTimeSelectionTable
    (
        option,
        phaseLimitStabilization,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::phaseLimitStabilization::readCoeffs()
{
    fieldName_ = coeffs_.lookup<word>("field");
    rateName_ = coeffs_.lookup<word>("rate");
    residualAlpha_ = coeffs_.lookup<scalar>("residualAlpha");
}


template<class Type>
void Foam::fv::phaseLimitStabilization::addSupType
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<Type>& eqn,
    const word& fieldName
) const
{
    const GeometricField<Type, fvPatchField, volMesh>& psi = eqn.psi();

    uniformDimensionedScalarField& rate =
        mesh_.lookupObjectRef<uniformDimensionedScalarField>(rateName_);

    eqn -= fvm::Sp(max(residualAlpha_ - alpha, scalar(0))*rho*rate, psi);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::phaseLimitStabilization::phaseLimitStabilization
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    option(name, modelType, dict, mesh),
    fieldName_(word::null),
    rateName_(word::null),
    residualAlpha_(NaN)
{
    readCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::phaseLimitStabilization::addSupFields() const
{
    return wordList(1, fieldName_);
}


FOR_ALL_FIELD_TYPES
(
    IMPLEMENT_FV_OPTION_ADD_ALPHA_RHO_SUP,
    phaseLimitStabilization
);


bool Foam::fv::phaseLimitStabilization::read(const dictionary& dict)
{
    if (option::read(dict))
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

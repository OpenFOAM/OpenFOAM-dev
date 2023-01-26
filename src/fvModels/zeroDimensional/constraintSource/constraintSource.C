/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "constraintSource.H"
#include "fluidThermo.H"
#include "fvModels.H"
#include "fvMatrix.H"
#include "fvmSup.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
namespace zeroDimensional
{
    defineTypeNameAndDebug(constraintSource, 0);
}
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::fv::zeroDimensional::constraintSource::addSupType
(
    fvMatrix<Type>& eqn,
    const word& fieldName
) const
{
    FatalErrorInFunction
        << "Cannot add a constraint source to field " << fieldName
        << " because this field's equation is not in mass-conservative form"
        << exit(FatalError);
}


void Foam::fv::zeroDimensional::constraintSource::addSupType
(
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    if (fieldName == "rho")
    {
        eqn += dmdt();
    }
    else
    {
        addSupType<scalar>(eqn, fieldName);
    }
}


template<class Type>
void Foam::fv::zeroDimensional::constraintSource::addSupType
(
    const volScalarField& rho,
    fvMatrix<Type>& eqn,
    const word& fieldName
) const
{
    eqn -= fvm::SuSp(-dmdt(), eqn.psi());
}


void Foam::fv::zeroDimensional::constraintSource::addSupType
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    if (fieldName == "rho")
    {
        eqn += dmdt();
    }
    else
    {
        addSupType<scalar>(rho, eqn, fieldName);
    }
}


template<class Type>
void Foam::fv::zeroDimensional::constraintSource::addSupType
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<Type>& eqn,
    const word& fieldName
) const
{
    FatalErrorInFunction
        << "Constraint sources do not support phase equations"
        << exit(FatalError);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::zeroDimensional::constraintSource::constraintSource
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    zeroDimensionalFvModel(name, modelType, mesh, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::zeroDimensional::constraintSource::~constraintSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::zeroDimensional::constraintSource::addsSupToField
(
    const word& fieldName
) const
{
    return true;
}


Foam::wordList Foam::fv::zeroDimensional::constraintSource::addSupFields() const
{
    return wordList(1, "rho");
}


FOR_ALL_FIELD_TYPES
(
    IMPLEMENT_FV_MODEL_ADD_SUP,
    fv::zeroDimensional::constraintSource
);


FOR_ALL_FIELD_TYPES
(
    IMPLEMENT_FV_MODEL_ADD_RHO_SUP,
    fv::zeroDimensional::constraintSource
);


FOR_ALL_FIELD_TYPES
(
    IMPLEMENT_FV_MODEL_ADD_ALPHA_RHO_SUP,
    fv::zeroDimensional::constraintSource
);


// ************************************************************************* //

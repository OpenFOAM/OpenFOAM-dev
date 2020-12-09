/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2020 OpenFOAM Foundation
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

#include "sixDoFAccelerationSource.H"
#include "fvMatrices.H"
#include "geometricOneField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(sixDoFAccelerationSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        sixDoFAccelerationSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "None.H"
#include "Constant.H"
#include "Uniform.H"
#include "ZeroConstant.H"
#include "OneConstant.H"
#include "Polynomial1.H"
#include "Sine.H"
#include "Square.H"
#include "Table.H"
#include "EmbeddedTableReader.H"
#include "FoamTableReader.H"
#include "Scale.H"
#include "CodedFunction1.H"

typedef Foam::fv::sixDoFAccelerationSource::accelerationVectors avType;

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

namespace Foam
{

    makeFunction1s(avType);

    defineTableReader(avType);
    makeTableReader(Embedded, avType);
    makeTableReader(Foam, avType);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::sixDoFAccelerationSource::sixDoFAccelerationSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    option(name, modelType, dict, mesh),
    accelerations_
    (
        Function1<accelerationVectors>::New("accelerations", coeffs_)
    ),
    UName_(coeffs_.lookupOrDefault<word>("U", "U")),
    g0_("g0", dimAcceleration, Zero)
{
    fieldNames_.setSize(1, UName_);
    applied_.setSize(1, false);

    if (mesh.foundObject<uniformDimensionedVectorField>("g"))
    {
        g0_ = mesh.lookupObject<uniformDimensionedVectorField>("g");
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::sixDoFAccelerationSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
) const
{
    addSup<geometricOneField>(geometricOneField(), eqn, fieldi);
}


void Foam::fv::sixDoFAccelerationSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
) const
{
    addSup<volScalarField>(rho, eqn, fieldi);
}


bool Foam::fv::sixDoFAccelerationSource::read(const dictionary& dict)
{
    if (option::read(dict))
    {
        accelerations_ = Function1<accelerationVectors>::New
        (
            "accelerations",
            dict
        );

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //

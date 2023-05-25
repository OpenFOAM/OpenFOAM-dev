/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2023 OpenFOAM Foundation
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

#include "liftModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(liftModel, 0);
    defineBlendedInterfacialModelTypeNameAndDebug(liftModel, 0);
    defineRunTimeSelectionTable(liftModel, dictionary);
}

const Foam::dimensionSet Foam::liftModel::dimF(1, -2, -2, 0, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::liftModel::liftModel
(
    const dictionary& dict,
    const phaseInterface& interface
)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::liftModel::~liftModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField> Foam::blendedLiftModel::F() const
{
    return evaluate(&liftModel::F, "F", liftModel::dimF, true);
}


Foam::tmp<Foam::surfaceScalarField> Foam::blendedLiftModel::Ff() const
{
    return evaluate(&liftModel::Ff, "Ff", liftModel::dimF*dimArea, true);
}


// ************************************************************************* //

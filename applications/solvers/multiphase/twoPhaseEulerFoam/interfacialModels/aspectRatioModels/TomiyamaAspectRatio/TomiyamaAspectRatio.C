/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 OpenFOAM Foundation
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

#include "TomiyamaAspectRatio.H"
#include "orderedPhasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace aspectRatioModels
{
    defineTypeNameAndDebug(TomiyamaAspectRatio, 0);
    addToRunTimeSelectionTable
    (
        aspectRatioModel,
        TomiyamaAspectRatio,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::aspectRatioModels::TomiyamaAspectRatio::TomiyamaAspectRatio
(
    const dictionary& dict,
    const orderedPhasePair& pair
)
:
    VakhrushevEfremov(dict, pair),
    wallDependentModel(pair.phase1().mesh())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::aspectRatioModels::TomiyamaAspectRatio::~TomiyamaAspectRatio()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::aspectRatioModels::TomiyamaAspectRatio::E() const
{
    return
        VakhrushevEfremov::E()
       *max
       (
           scalar(1) - 0.35*yWall()/pair_.dispersed().d(),
           scalar(0.65)
       );
}


// ************************************************************************* //

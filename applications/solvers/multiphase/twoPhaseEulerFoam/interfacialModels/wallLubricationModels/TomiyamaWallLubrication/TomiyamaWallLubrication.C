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

#include "TomiyamaWallLubrication.H"
#include "phasePair.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace wallLubricationModels
{
    defineTypeNameAndDebug(TomiyamaWallLubrication, 0);
    addToRunTimeSelectionTable
    (
        wallLubricationModel,
        TomiyamaWallLubrication,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallLubricationModels::TomiyamaWallLubrication::TomiyamaWallLubrication
(
    const dictionary& dict,
    const phasePair& pair
)
:
    wallLubricationModel(dict, pair),
    D_("Cwd", dimLength, dict.lookup("D"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallLubricationModels::TomiyamaWallLubrication::~TomiyamaWallLubrication()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField>
Foam::wallLubricationModels::TomiyamaWallLubrication::F() const
{
    volVectorField Ur(pair_.Ur());

    const volVectorField& n(nWall());
    const volScalarField& y(yWall());

    volScalarField Eo(pair_.Eo());

    return
        (
            pos(Eo - 1.0)*neg(Eo - 5.0)*exp(-0.933*Eo + 0.179)
          + pos(Eo - 5.0)*neg(Eo - 33.0)*(0.00599*Eo - 0.0187)
          + pos(Eo - 33.0)*0.179
        )
       *0.5
       *pair_.dispersed().d()
       *(
            1/sqr(y)
          - 1/sqr(D_ - y)
        )
       *pair_.dispersed()
       *pair_.continuous().rho()
       *magSqr(Ur - (Ur & n)*n)
       *n;
}


// ************************************************************************* //

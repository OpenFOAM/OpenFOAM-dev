/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2019 OpenFOAM Foundation
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

#include "cosine.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace partitioningModels
{
    defineTypeNameAndDebug(cosine, 0);
    addToRunTimeSelectionTable
    (
        partitioningModel,
        cosine,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallBoilingModels::partitioningModels::
cosine::cosine(const dictionary& dict)
:
    partitioningModel(),
    alphaLiquid1_(readScalar(dict.lookup("alphaLiquid1"))),
    alphaLiquid0_(readScalar(dict.lookup("alphaLiquid0")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallBoilingModels::partitioningModels::
cosine::~cosine()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::wallBoilingModels::partitioningModels::
cosine::fLiquid
(
    const scalarField& alphaLiquid
) const
{
    return
        pos0(alphaLiquid1_ - alphaLiquid)
       *(
            neg(alphaLiquid0_ - alphaLiquid)
           *(
                0.5
               *(
                    1 - cos
                    (
                        constant::mathematical::pi
                       *(alphaLiquid - alphaLiquid0_)
                       /(alphaLiquid1_ - alphaLiquid0_)
                    )
                )
            )
        )
      + neg(alphaLiquid1_ - alphaLiquid);
}


void Foam::wallBoilingModels::partitioningModels::
cosine::write(Ostream& os) const
{
    partitioningModel::write(os);
    writeEntry(os, "alphaLiquid1", alphaLiquid1_);
    writeEntry(os, "alphaLiquid0", alphaLiquid0_);
}


// ************************************************************************* //

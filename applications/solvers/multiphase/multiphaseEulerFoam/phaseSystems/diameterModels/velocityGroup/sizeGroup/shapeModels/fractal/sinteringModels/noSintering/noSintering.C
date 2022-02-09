/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2022 OpenFOAM Foundation
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

#include "noSintering.H"
#include "addToRunTimeSelectionTable.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace shapeModels
{
namespace sinteringModels
{
    defineTypeNameAndDebug(noSintering, 0);
    addToRunTimeSelectionTable(sinteringModel, noSintering, dictionary);
}
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::shapeModels::sinteringModels::noSintering::noSintering
(
    const dictionary& dict,
    const fractal& fractalShape
)
:
    sinteringModel(dict, fractalShape)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModels::shapeModels::sinteringModels::noSintering::~noSintering()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvScalarMatrix>
Foam::diameterModels::shapeModels::sinteringModels::noSintering::R() const
{
    const sizeGroup& fi = fractal_.SizeGroup();

    volScalarField::Internal R
    (
        IOobject
        (
            "noSintering:R",
            fi.time().timeName(),
            fi.mesh()
        ),
        fi.mesh(),
        dimensionedScalar(inv(dimTime*dimLength), 0)
    );

    return fvm::Su(R, fractal_.fld());
}


// ************************************************************************* //

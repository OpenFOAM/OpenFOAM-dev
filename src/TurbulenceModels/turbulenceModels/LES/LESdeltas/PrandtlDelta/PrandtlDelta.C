/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "PrandtlDelta.H"
#include "wallDist.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(PrandtlDelta, 0);
    addToRunTimeSelectionTable(LESdelta, PrandtlDelta, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::PrandtlDelta::calcDelta()
{
    delta_ = min
    (
        static_cast<const volScalarField&>(geometricDelta_()),
        (kappa_/Cdelta_)*wallDist(mesh_).y()
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PrandtlDelta::PrandtlDelta
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    LESdelta(name, mesh),
    geometricDelta_(LESdelta::New(name, mesh, dict.subDict(type() + "Coeffs"))),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    Cdelta_
    (
        dict.subDict(type() + "Coeffs").lookupOrDefault<scalar>("Cdelta", 0.158)
    )
{
    calcDelta();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::PrandtlDelta::read(const dictionary& dict)
{
    const dictionary& coeffDict(dict.subDict(type() + "Coeffs"));

    geometricDelta_().read(coeffDict);
    dict.readIfPresent<scalar>("kappa", kappa_);
    coeffDict.readIfPresent<scalar>("Cdelta", Cdelta_);
    calcDelta();
}


void Foam::PrandtlDelta::correct()
{
    geometricDelta_().correct();

    if (mesh_.changing())
    {
        calcDelta();
    }
}


// ************************************************************************* //

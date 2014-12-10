/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(PrandtlDelta, 0);
addToRunTimeSelectionTable(LESdelta, PrandtlDelta, dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void PrandtlDelta::calcDelta()
{
    delta_ = min
    (
        static_cast<const volScalarField&>(geometricDelta_()),
        (kappa_/Cdelta_)*wallDist(mesh_).y()
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

PrandtlDelta::PrandtlDelta
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dd
)
:
    LESdelta(name, mesh),
    geometricDelta_(LESdelta::New(name, mesh, dd.subDict(type() + "Coeffs"))),
    kappa_(dd.lookupOrDefault<scalar>("kappa", 0.41)),
    Cdelta_
    (
        dd.subDict(type() + "Coeffs").lookupOrDefault<scalar>("Cdelta", 0.158)
    )
{
    calcDelta();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void PrandtlDelta::read(const dictionary& d)
{
    const dictionary& dd(d.subDict(type() + "Coeffs"));

    geometricDelta_().read(dd);
    d.readIfPresent<scalar>("kappa", kappa_);
    dd.readIfPresent<scalar>("Cdelta", Cdelta_);
    calcDelta();
}


void PrandtlDelta::correct()
{
    geometricDelta_().correct();

    if (mesh_.changing())
    {
        calcDelta();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

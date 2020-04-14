/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "maxDeltaxyz.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{
    defineTypeNameAndDebug(maxDeltaxyz, 0);
    addToRunTimeSelectionTable(LESdelta, maxDeltaxyz, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::LESModels::maxDeltaxyz::calcDelta()
{
    const fvMesh& mesh = momentumTransportModel_.mesh();

    label nD = mesh.nGeometricD();

    const cellList& cells = mesh.cells();
    scalarField hmax(cells.size());

    forAll(cells,celli)
    {
        scalar deltaMaxTmp = 0.0;
        const labelList& cFaces = mesh.cells()[celli];
        const point& centrevector = mesh.cellCentres()[celli];

        forAll(cFaces, cFacei)
        {
            label facei = cFaces[cFacei];
            const point& facevector = mesh.faceCentres()[facei];
            scalar tmp = mag(facevector - centrevector);
            if (tmp > deltaMaxTmp)
            {
                deltaMaxTmp = tmp;
            }
        }

        hmax[celli] = deltaCoeff_*deltaMaxTmp;
    }

    if (nD == 3)
    {
        delta_.primitiveFieldRef() = hmax;
    }
    else if (nD == 2)
    {
        WarningInFunction
            << "Case is 2D, LES is not strictly applicable\n"
            << endl;

        delta_.primitiveFieldRef() = hmax;
    }
    else
    {
        FatalErrorInFunction
            << "Case is not 3D or 2D, LES is not applicable"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LESModels::maxDeltaxyz::maxDeltaxyz
(
    const word& name,
    const momentumTransportModel& turbulence,
    const dictionary& dict
)
:
    LESdelta(name, turbulence),
    deltaCoeff_
    (
        dict.optionalSubDict(type() + "Coeffs").lookupOrDefault<scalar>
        (
            "deltaCoeff",
            1
        )
    )
{
    calcDelta();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::LESModels::maxDeltaxyz::read(const dictionary& dict)
{
    const dictionary& coeffsDict(dict.optionalSubDict(type() + "Coeffs"));

    coeffsDict.readIfPresent<scalar>("deltaCoeff", deltaCoeff_);

    calcDelta();
}


void Foam::LESModels::maxDeltaxyz::correct()
{
    if (momentumTransportModel_.mesh().changing())
    {
        calcDelta();
    }
}


// ************************************************************************* //

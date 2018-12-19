/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "IDDESDelta.H"
#include "addToRunTimeSelectionTable.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{
    defineTypeNameAndDebug(IDDESDelta, 0);
    addToRunTimeSelectionTable(LESdelta, IDDESDelta, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::LESModels::IDDESDelta::calcDelta()
{
    const volScalarField& hmax = hmax_;
    const fvMesh& mesh = turbulenceModel_.mesh();

    // Wall-normal vectors
    const volVectorField& n = wallDist::New(mesh).n();

    tmp<volScalarField> tfaceToFacenMax
    (
        volScalarField::New
        (
            "faceToFaceMax",
            mesh,
            dimensionedScalar(dimLength, 0)
        )
    );

    scalarField& faceToFacenMax = tfaceToFacenMax.ref().primitiveFieldRef();

    const cellList& cells = mesh.cells();
    const vectorField& faceCentres = mesh.faceCentres();

    forAll(cells, celli)
    {
        scalar maxDelta = 0.0;
        const labelList& cFaces = cells[celli];
        const vector nci = n[celli];

        forAll(cFaces, cFacei)
        {
            label facei = cFaces[cFacei];
            const point& fci = faceCentres[facei];

            forAll(cFaces, cFacej)
            {
                label facej = cFaces[cFacej];
                const point& fcj = faceCentres[facej];
                scalar ndfc = nci & (fcj - fci);

                if (ndfc > maxDelta)
                {
                    maxDelta = ndfc;
                }
            }
        }

        faceToFacenMax[celli] = maxDelta;
    }


    label nD = mesh.nGeometricD();

    if (nD == 2)
    {
        WarningInFunction
            << "Case is 2D, LES is not strictly applicable" << nl
            << endl;
    }
    else if (nD != 3)
    {
        FatalErrorInFunction
            << "Case must be either 2D or 3D" << exit(FatalError);
    }

    delta_.primitiveFieldRef() =
        min
        (
            max
            (
                max
                (
                    Cw_*wallDist::New(mesh).y(),
                    Cw_*hmax
                ),
                tfaceToFacenMax
            ),
            hmax
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LESModels::IDDESDelta::IDDESDelta
(
    const word& name,
    const turbulenceModel& turbulence,
    const dictionary& dict
)
:
    LESdelta(name, turbulence),
    hmax_
    (
        IOobject::groupName("hmax", turbulence.U().group()),
        turbulence,
        dict
    ),
    Cw_
    (
        dict.optionalSubDict(type() + "Coeffs").lookupOrDefault<scalar>
        (
            "Cw",
            0.15
        )
    )
{
    calcDelta();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::LESModels::IDDESDelta::read(const dictionary& dict)
{
    const dictionary& coeffsDict(dict.optionalSubDict(type() + "Coeffs"));

    coeffsDict.readIfPresent<scalar>("Cw", Cw_);

    calcDelta();
}


void Foam::LESModels::IDDESDelta::correct()
{
    if (turbulenceModel_.mesh().changing())
    {
        hmax_.correct();
        calcDelta();
    }
}


// ************************************************************************* //

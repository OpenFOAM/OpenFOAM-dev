/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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

#include "interfaceHeight.H"
#include "fvMesh.H"
#include "interpolation.H"
#include "IOmanip.H"
#include "meshSearch.H"
#include "lineCellFace.H"
#include "Time.H"
#include "uniformDimensionedFields.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(interfaceHeight, 0);
    addToRunTimeSelectionTable(functionObject, interfaceHeight, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::interfaceHeight::writePositions()
{
    const uniformDimensionedVectorField& g =
        mesh_.lookupObject<uniformDimensionedVectorField>("g");
    const vector gHat = g.value()/mag(g.value());

    const volScalarField& alpha =
        mesh_.lookupObject<volScalarField>(alphaName_);

    autoPtr<interpolation<scalar>>
        interpolator
        (
            interpolation<scalar>::New(interpolationScheme_, alpha)
        );

    if (Pstream::master())
    {
        writeTime(file(HEIGHT_FILE));
        writeTime(file(POSITION_FILE));
    }

    forAll(locations_, li)
    {
        // Create a set along a ray projected in the direction of gravity
        const sampledSets::lineCellFace set
        (
            "",
            mesh_,
            meshSearch(mesh_),
            "xyz",
            locations_[li] + gHat*mesh_.bounds().mag(),
            locations_[li] - gHat*mesh_.bounds().mag()
        );

        // Find the height of the location above the boundary
        scalar hLB = set.size() ? - gHat & (locations_[li] - set[0]) : - vGreat;
        reduce(hLB, maxOp<scalar>());

        // Calculate the integrals of length and length*alpha along the sampling
        // line. The latter is equal to the equivalent length with alpha equal
        // to one.
        scalar sumLength = 0, sumLengthAlpha = 0;
        for(label si = 0; si < set.size() - 1; ++ si)
        {
            if (set.segments()[si] != set.segments()[si+1])
            {
                continue;
            }

            const vector& p0 = set[si], p1 = set[si+1];
            const label c0 = set.cells()[si], c1 = set.cells()[si+1];
            const label f0 = set.faces()[si], f1 = set.faces()[si+1];
            const scalar a0 = interpolator->interpolate(p0, c0, f0);
            const scalar a1 = interpolator->interpolate(p1, c1, f1);

            const scalar l = - gHat & (p1 - p0);
            sumLength += l;
            sumLengthAlpha += l*(a0 + a1)/2;
        }

        reduce(sumLength, sumOp<scalar>());
        reduce(sumLengthAlpha, sumOp<scalar>());

        // Write out
        if (Pstream::master())
        {
            // Interface heights above the boundary and location
            const scalar hIB =
                liquid_ ? sumLengthAlpha : sumLength - sumLengthAlpha;
            const scalar hIL = hIB - hLB;

            // Position of the interface
            const point p = locations_[li] - gHat*hIL;

            const Foam::Omanip<int> w = valueWidth(1);

            file(HEIGHT_FILE) << w << hIB << w << hIL;
            file(POSITION_FILE) << '(' << w << p.x() << w << p.y()
                << valueWidth() << p.z() << ") ";
        }
    }

    if (Pstream::master())
    {
        file(HEIGHT_FILE).endl();
        file(POSITION_FILE).endl();
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::interfaceHeight::writeFileHeader(const label i)
{
    forAll(locations_, li)
    {
        writeHeaderValue
        (
            file(i),
            "Location " + Foam::name(li),
            locations_[li]
        );
    }

    switch (fileID(i))
    {
        case HEIGHT_FILE:
            writeHeaderValue
            (
                file(i),
                "hB",
                "Interface height above the boundary"
            );
            writeHeaderValue
            (
                file(i),
                "hL",
                "Interface height above the location"
            );
            break;
        case POSITION_FILE:
            writeHeaderValue(file(i), "p", "Interface position");
            break;
    }

    const Foam::Omanip<int> w = valueWidth(1);

    writeCommented(file(i), "Location");
    forAll(locations_, li)
    {
        switch (fileID(i))
        {
            case HEIGHT_FILE:
                file(i) << w << li << w << ' ';
                break;
            case POSITION_FILE:
                file(i) << w << li << w << ' ' << w << ' ' << "  ";
                break;
        }
    }
    file(i).endl();

    writeCommented(file(i), "Time");
    forAll(locations_, li)
    {
        switch (fileID(i))
        {
            case HEIGHT_FILE:
                file(i) << w << "hB" << w << "hL";
                break;
            case POSITION_FILE:
                file(i) << w << "p" << w << ' ' << w << ' ' << "  ";
                break;
        }
    }
    file(i).endl();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::interfaceHeight::interfaceHeight
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(mesh_, name),
    alphaName_("alpha"),
    liquid_(true),
    locations_(),
    interpolationScheme_("cellPoint")
{
    read(dict);

    wordList names(2);
    names[HEIGHT_FILE] = "height";
    names[POSITION_FILE] = "position";

    resetNames(wordList(names));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::interfaceHeight::~interfaceHeight()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::interfaceHeight::read(const dictionary& dict)
{
    dict.readIfPresent("alpha", alphaName_);
    dict.readIfPresent("liquid", liquid_);
    dict.lookup("locations") >> locations_;
    dict.readIfPresent("interpolationScheme", interpolationScheme_);

    return true;
}


bool Foam::functionObjects::interfaceHeight::execute()
{
    return true;
}


bool Foam::functionObjects::interfaceHeight::end()
{
    return true;
}


bool Foam::functionObjects::interfaceHeight::write()
{
    logFiles::write();

    writePositions();

    return true;
}


// ************************************************************************* //

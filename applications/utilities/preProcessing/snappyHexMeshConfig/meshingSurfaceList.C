/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "meshingSurfaceList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::meshingSurfaceList::mergeBoundingBoxes
(
    boundBox& bb1,
    const boundBox& bb2
)
{
    if (bb1.volume() == 0)
    {
        bb1 = bb2;
        return;
    }

    point& min1 = bb1.min();
    point& max1 = bb1.max();
    const point& min2 = bb2.min();
    const point& max2 = bb2.max();

    forAll(min1, i)
    {
        min1[i] = Foam::min(min1[i], min2[i]);
        max1[i] = Foam::max(max1[i], max2[i]);
    }
}


void Foam::meshingSurfaceList::swapExternalIndexZero(const label index)
{
    if (index == 0)
    {
        return;
    }

    autoPtr<meshingSurface> s0Ptr(set(0, nullptr));
    autoPtr<meshingSurface> sIndexPtr(set(index, nullptr));

    set(index, s0Ptr.ptr());
    set(0, sIndexPtr.ptr());
}


bool Foam::meshingSurfaceList::regionsValid
(
    const wordList& specifiedRegions,
    const wordList& regions,
    const word& opt
)
{
    if (specifiedRegions.empty())
    {
        return false;
    }

    forAll(specifiedRegions, s)
    {
        bool match(false);

        forAll(regions, r)
        {
            match = (specifiedRegions[s] == regions[r]) || match;
            if (match)
            {
                break;
            }
        }

        if (!match)
        {
            FatalErrorInFunction
                << "Region '"<< specifiedRegions[s]
                << "' specified with the '" << opt
                << "' option" << nl
                << "does not match any regions in the external surface"
                << exit(FatalError);
        }
    }

    return true;
}


void Foam::meshingSurfaceList::setSurfaceTypes
(
    const List<word>& surfaces,
    const surfaceType& type
)
{
    forAll(surfaces, i)
    {
        const word surface = fileName(surfaces[i]).lessExt();
        const word surfaceType = meshingSurface::surfaceTypeNames[type];
        bool match = false;

        // Check if the surface name matches a surface name
        forAll(*this, surfi)
        {
            if (surface == operator[](surfi).name())
            {
                // For cellZone and rotatingZone, ensure surface is closed
                if
                (
                    !operator[](surfi).closed()
                 && (
                        type == surfaceType::rotatingZone
                     || type == surfaceType::cellZone
                    )
                )
                {
                    FatalErrorInFunction
                        << "Argument to '-" << surfaceType
                        << "' contains the " << surfaceType << " '"
                        << surface << "'" << nl
                        << "which is not a closed surface."
                        << exit(FatalError);
                }

                operator[](surfi).type() = type;
                match = true;
                break;
            }
        }

        if (!match)
        {
            FatalErrorInFunction
                << "Argument to '-" << surfaceType
                << "' contains the " << surfaceType << " '"
                << surface << "'" << nl
                << "which does not correspond to any surfaces."
                << exit(FatalError);
        }
    }
}


void Foam::meshingSurfaceList::setRotatingZoneBounds()
{
    forAll(*this, surfi)
    {
        if (operator[](surfi).type() == surfaceType::rotatingZone)
        {
            mergeBoundingBoxes(rzbb_, operator[](surfi).bb());
        }
    }
}


void Foam::meshingSurfaceList::identifyCellZones()
{
    forAll(*this, i)
    {
        if
        (
            operator[](i).type() == surfaceType::external
         || operator[](i).type() == surfaceType::rotatingZone
         || operator[](i).nParts() != 1
         || !operator[](i).closed()
        )
        {
            continue;
        }

        forAll(*this, j)
        {
            if (i == j)
            {
                continue;
            }

            if (operator[](i).bb().contains(operator[](j).bb()))
            {
                operator[](i).type() = surfaceType::cellZone;
                continue;
            }
        }
    }
}


void Foam::meshingSurfaceList::reportWordList(const wordList& wl)
{
    forAll(wl, i)
    {
        Info<< wl[i];
        if (i != wl.size() - 1)
        {
            Info<< ", ";
        }
    }

    Info<< endl;
}


void Foam::meshingSurfaceList::reportSurfaces()
{
    Info<< "Case contains the following surface geometry files:"
        << endl;

    forAll(*this, i)
    {
        Info<< "+ " << operator[](i).name() << nl
            << "  + File: " << operator[](i).file() << nl
            << "  + Bounding box: " << operator[](i).bb() << nl
            << "  + " << (operator[](i).closed() ? "Closed" : "Open")
            << " surface" << endl;

        switch (operator[](i).type())
        {
            case surfaceType::external:
            {
                Info << "  + External boundary surface" << endl;

                if (!operator[](i).inletRegions().empty())
                {
                    Info<< "  + Inlet regions: ";
                    reportWordList(operator[](i).inletRegions());
                }

                if (!operator[](i).outletRegions().empty())
                {
                    Info << "  + Outlet regions: ";
                    reportWordList(operator[](i).outletRegions());
                }

                break;
            }

            case surfaceType::wall:
            {
                Info << "  + Wall boundary surface" << endl;
                break;
            }

            case surfaceType::cellZone:
            {
                Info << "  + Cell zone surface" << endl;
                break;
            }

            case surfaceType::rotatingZone:
            {
                Info << "  + Rotating zone surface" << endl;
                break;
            }

            case surfaceType::baffle:
            {
                Info << "  + Baffle wall surface" << endl;
                break;
            }
        }
    }

    Info<< endl;
}


void Foam::meshingSurfaceList::setBounds(const boundBox& bb)
{
    if (bb.contains(bb_))
    {
        Info<< "Specifed bounding box contains the overall bounding box"
            << endl;
    }
    else
    {
        WarningInFunction
            << "Specified bounding box does not contain the overall "
            << "bounding box"
            << endl;
    }

    if (bb.mag() > 1.5*bb_.mag())
    {
        Info<< "Specified bounding box is > (1.5 * overall "
            << "bounding box)" << nl
            << "**Assuming this is an external flow**"
            << endl;

        operator[](0).type() = surfaceType::wall;
    }

    Info<< endl;

    bb_ = bb;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meshingSurfaceList::meshingSurfaceList
(
    const Time& time,
    const fileNameList& surfaces,
    const wordList& cellZones,
    const wordList& rotatingZones,
    const wordList& baffles,
    const boundBox& bb,
    const wordList& specifiedInletRegions,
    const wordList& specifiedOutletRegions
)
:
    PtrList<meshingSurface>(),
    bb_(),
    rzbb_()
{
    // Load all the surfaces and construct the bounding box
    forAll(surfaces, i)
    {
        append(new meshingSurface(surfaces[i], time));

        mergeBoundingBoxes(bb_, last().bb());
    }

    label externalID = 0;
    forAll(surfaces, i)
    {
        // Test for external surface inflates the bounding box by a small factor
        // to account for imprecise surface geometry files
        boundBox bbInflate = operator[](i).bb();
        bbInflate.inflate(1e-4);

        if
        (
            operator[](i).closed()
         && operator[](i).nParts() == 1
         && bbInflate.contains(bb_)
        )
        {
            externalID = i;
            operator[](i).type() = surfaceType::external;

            // Override any inlet and outlet regions on external boundary
            // specified by the '-inletRegions' and '-outletRegions' options
            const wordList& regions = operator[](i).regions();

            wordList& inletRegions = operator[](i).inletRegions();
            if (regionsValid(specifiedInletRegions, regions, "-inletRegions"))
            {
                inletRegions = specifiedInletRegions;
            }

            wordList& outletRegions = operator[](i).outletRegions();
            if (regionsValid(specifiedOutletRegions, regions, "-outletRegions"))
            {
                outletRegions = specifiedOutletRegions;
            }

            // If inletRegions and outletRegions are both empty, set "template"
            // names
            if (inletRegions.empty() && outletRegions.empty())
            {
                inletRegions.append("<inletRegion>");
                outletRegions.append("<outletRegion>");
            }
        }
    }

    swapExternalIndexZero(externalID);

    if (!rotatingZones.empty())
    {
        setSurfaceTypes(rotatingZones, surfaceType::rotatingZone);
    }

    if (!baffles.empty())
    {
        setSurfaceTypes(baffles, surfaceType::baffle);
    }

    setRotatingZoneBounds();

    if (cellZones.empty())
    {
        identifyCellZones();
    }
    else
    {
        setSurfaceTypes(cellZones, surfaceType::cellZone);
    }

    reportSurfaces();

    if (bb.volume() != 0)
    {
        setBounds(bb);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::meshingSurfaceList::~meshingSurfaceList()
{}


// ************************************************************************* //

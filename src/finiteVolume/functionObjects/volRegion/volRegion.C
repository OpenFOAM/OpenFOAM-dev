/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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

#include "volRegion.H"
#include "volMesh.H"
#include "globalMeshData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(volRegion, 0);
}
}

template<>
const char*
Foam::NamedEnum
<
    Foam::functionObjects::volRegion::regionTypes,
    2
>::names[] = {"cellZone", "all"};

const Foam::NamedEnum
<
    Foam::functionObjects::volRegion::regionTypes,
    2
> Foam::functionObjects::volRegion::regionTypeNames_;


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::volRegion::writeFileHeader
(
    const writeFile& wf,
    Ostream& file
)
{
    wf.writeCommented(file, "Region");
    file<< setw(1) << ':' << setw(1) << ' '
        << regionTypeNames_[regionType_] << " " << regionName_ << endl;
    wf.writeHeaderValue(file, "Cells", nCells_);
    wf.writeHeaderValue(file, "Volume", V_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::volRegion::volRegion
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    regionType_
    (
        dict.found("regionType")
      ? regionTypeNames_.read(dict.lookup("regionType"))
      : vrtAll
    ),
    regionID_(-1)
{
    read(dict);

    // Cache integral properties of the region for writeFileHeader
    nCells_ = nCells();
    V_ = V();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::volRegion::~volRegion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::volRegion::read
(
    const dictionary& dict
)
{
    switch (regionType_)
    {
        case vrtCellZone:
        {
            dict.lookup("name") >> regionName_;

            regionID_ = mesh_.cellZones().findZoneID(regionName_);

            if (regionID_ < 0)
            {
                FatalIOErrorInFunction(dict)
                    << "Unknown cell zone name: " << regionName_
                    << ". Valid cell zones are: " << mesh_.cellZones().names()
                    << exit(FatalIOError);
            }

            if (nCells() == 0)
            {
                FatalIOErrorInFunction(dict)
                    << regionTypeNames_[regionType_]
                    << "(" << regionName_ << "):" << nl
                    << "    Region has no cells"
                    << exit(FatalIOError);
            }

            break;
        }

        case vrtAll:
        {
            break;
        }

        default:
        {
            FatalIOErrorInFunction(dict)
                << "Unknown region type. Valid region types are:"
                << regionTypeNames_
                << exit(FatalIOError);
        }
    }

    return true;
}


const Foam::labelList& Foam::functionObjects::volRegion::cellIDs() const
{
    if (regionType_ == vrtAll)
    {
        return labelList::null();
    }
    else
    {
        return mesh_.cellZones()[regionID_];
    }
}


Foam::label Foam::functionObjects::volRegion::nCells() const
{
    if (regionType_ == vrtAll)
    {
        return mesh_.globalData().nTotalCells();
    }
    else
    {
        return returnReduce(cellIDs().size(), sumOp<label>());
    }
}


Foam::scalar Foam::functionObjects::volRegion::V() const
{
    if (regionType_ == vrtAll)
    {
        return gSum(mesh_.V());
    }
    else
    {
        return gSum(scalarField(mesh_.V(), cellIDs()));
    }
}


// ************************************************************************* //

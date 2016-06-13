/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
#include "fvMesh.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace fieldValues
{
    defineTypeNameAndDebug(volRegion, 0);
    addToRunTimeSelectionTable(fieldValue, volRegion, dictionary);
    addToRunTimeSelectionTable(functionObject, volRegion, dictionary);
}
}
}

template<>
const char*
Foam::NamedEnum
<
    Foam::functionObjects::fieldValues::volRegion::regionTypes,
    2
>::names[] = {"cellZone", "all"};

template<>
const char*
Foam::NamedEnum
<
    Foam::functionObjects::fieldValues::volRegion::operationType,
    11
>::names[] =
{
    "none",
    "sum",
    "sumMag",
    "average",
    "weightedAverage",
    "volAverage",
    "weightedVolAverage",
    "volIntegrate",
    "min",
    "max",
    "CoV"
};

const Foam::NamedEnum
<
    Foam::functionObjects::fieldValues::volRegion::regionTypes,
    2
> Foam::functionObjects::fieldValues::volRegion::regionTypeNames_;

const Foam::NamedEnum
<
    Foam::functionObjects::fieldValues::volRegion::operationType,
    11
> Foam::functionObjects::fieldValues::volRegion::operationTypeNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::fieldValues::volRegion::setCellZoneCells()
{
    switch (regionType_)
    {
        case stCellZone:
        {
            dict().lookup("name") >> regionName_;

            label zoneId = mesh().cellZones().findZoneID(regionName_);

            if (zoneId < 0)
            {
                FatalErrorInFunction
                    << "Unknown cell zone name: " << regionName_
                    << ". Valid cell zones are: " << mesh().cellZones().names()
                    << nl << exit(FatalError);
            }

            cellId_ = mesh().cellZones()[zoneId];
            nCells_ = returnReduce(cellId_.size(), sumOp<label>());
            break;
        }

        case stAll:
        {
            cellId_ = identity(mesh().nCells());
            nCells_ = returnReduce(cellId_.size(), sumOp<label>());
            break;
        }

        default:
        {
            FatalErrorInFunction
               << "Unknown region type. Valid region types are:"
                << regionTypeNames_ << nl << exit(FatalError);
        }
    }

    if (debug)
    {
        Pout<< "Selected region size = " << cellId_.size() << endl;
    }
}


Foam::scalar Foam::functionObjects::fieldValues::volRegion::volume() const
{
    return gSum(filterField(mesh().V()));
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::fieldValues::volRegion::initialise
(
    const dictionary& dict
)
{
    setCellZoneCells();

    if (nCells_ == 0)
    {
        FatalErrorInFunction
            << type() << " " << name() << ": "
            << regionTypeNames_[regionType_] << "(" << regionName_ << "):" << nl
            << "    Region has no cells" << exit(FatalError);
    }

    volume_ = volume();

    Info<< type() << " " << name() << ":"
        << regionTypeNames_[regionType_] << "(" << regionName_ << "):" << nl
        << "    total cells  = " << nCells_ << nl
        << "    total volume = " << volume_
        << nl << endl;

    if (dict.readIfPresent("weightField", weightFieldName_))
    {
        Info<< "    weight field = " << weightFieldName_;
    }

    Info<< nl << endl;
}


void Foam::functionObjects::fieldValues::volRegion::writeFileHeader
(
    const label i
)
{
    writeCommented(file(), "Region type : ");
    file() << regionTypeNames_[regionType_] << " " << regionName_ << endl;
    writeCommented(file(), "Cells  : ");
    file() << nCells_ << endl;
    writeCommented(file(), "Volume : ");
    file() << volume_ << endl;

    writeCommented(file(), "Time");
    if (writeVolume_)
    {
        file() << tab << "Volume";
    }

    forAll(fields_, i)
    {
        file()
            << tab << operationTypeNames_[operation_]
            << "(" << fields_[i] << ")";
    }

    file() << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldValues::volRegion::volRegion
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldValue(name, runTime, dict, typeName),
    regionType_(regionTypeNames_.read(dict.lookup("regionType"))),
    operation_(operationTypeNames_.read(dict.lookup("operation"))),
    nCells_(0),
    cellId_(),
    weightFieldName_("none"),
    writeVolume_(dict.lookupOrDefault("writeVolume", false))
{
    if (!isA<fvMesh>(obr_))
    {
        FatalErrorInFunction
            << "objectRegistry is not an fvMesh" << exit(FatalError);
    }

    read(dict);
}


Foam::functionObjects::fieldValues::volRegion::volRegion
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    fieldValue(name, obr, dict, typeName),
    regionType_(regionTypeNames_.read(dict.lookup("regionType"))),
    operation_(operationTypeNames_.read(dict.lookup("operation"))),
    nCells_(0),
    cellId_(),
    weightFieldName_("none"),
    writeVolume_(dict.lookupOrDefault("writeVolume", false))
{
    if (!isA<fvMesh>(obr_))
    {
        FatalErrorInFunction
            << "objectRegistry is not an fvMesh" << exit(FatalError);
    }

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fieldValues::volRegion::~volRegion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fieldValues::volRegion::read
(
    const dictionary& dict
)
{
    fieldValue::read(dict);

    // No additional info to read
    initialise(dict);

    return true;
}


bool Foam::functionObjects::fieldValues::volRegion::write()
{
    fieldValue::write();

    if (Pstream::master())
    {
        writeTime(file());
    }

    if (writeVolume_)
    {
        volume_ = volume();
        if (Pstream::master())
        {
            file() << tab << volume_;
        }
        Log << "    total volume = " << volume_ << endl;
    }

    forAll(fields_, i)
    {
        const word& fieldName = fields_[i];
        bool processed = false;

        processed = processed || writeValues<scalar>(fieldName);
        processed = processed || writeValues<vector>(fieldName);
        processed = processed || writeValues<sphericalTensor>(fieldName);
        processed = processed || writeValues<symmTensor>(fieldName);
        processed = processed || writeValues<tensor>(fieldName);

        if (!processed)
        {
            WarningInFunction
                << "Requested field " << fieldName
                << " not found in database and not processed"
                << endl;
        }
    }

    if (Pstream::master())
    {
        file()<< endl;
    }

    Log << endl;

    return true;
}


// ************************************************************************* //

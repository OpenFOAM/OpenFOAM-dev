/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "cellSource.H"
#include "fvMesh.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* NamedEnum<fieldValues::cellSource::sourceType, 2>::names[] =
    {
        "cellZone",
        "all"
    };


    template<>
    const char* NamedEnum<fieldValues::cellSource::operationType, 11>::names[] =
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

    namespace fieldValues
    {
        defineTypeNameAndDebug(cellSource, 0);
        addToRunTimeSelectionTable(fieldValue, cellSource, dictionary);
    }
}


const Foam::NamedEnum<Foam::fieldValues::cellSource::sourceType, 2>
    Foam::fieldValues::cellSource::sourceTypeNames_;

const Foam::NamedEnum<Foam::fieldValues::cellSource::operationType, 11>
    Foam::fieldValues::cellSource::operationTypeNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fieldValues::cellSource::setCellZoneCells()
{
    switch (source_)
    {
        case stCellZone:
        {
            dict().lookup("sourceName") >> sourceName_;

            label zoneId = mesh().cellZones().findZoneID(sourceName_);

            if (zoneId < 0)
            {
                FatalErrorInFunction
                    << "Unknown cell zone name: " << sourceName_
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
               << "Unknown source type. Valid source types are:"
                << sourceTypeNames_ << nl << exit(FatalError);
        }
    }

    if (debug)
    {
        Pout<< "Selected source size = " << cellId_.size() << endl;
    }
}


Foam::scalar Foam::fieldValues::cellSource::volume() const
{
    return gSum(filterField(mesh().V()));
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fieldValues::cellSource::initialise(const dictionary& dict)
{
    setCellZoneCells();

    if (nCells_ == 0)
    {
        WarningInFunction
            << type() << " " << name_ << ": "
            << sourceTypeNames_[source_] << "(" << sourceName_ << "):" << nl
            << "    Source has no cells - deactivating" << endl;

        active_ = false;
        return;
    }

    volume_ = volume();

    Info<< type() << " " << name_ << ":"
        << sourceTypeNames_[source_] << "(" << sourceName_ << "):" << nl
        << "    total cells  = " << nCells_ << nl
        << "    total volume = " << volume_
        << nl << endl;

    if (dict.readIfPresent("weightField", weightFieldName_))
    {
        Info<< "    weight field = " << weightFieldName_;
    }

    Info<< nl << endl;
}


void Foam::fieldValues::cellSource::writeFileHeader(const label i)
{
    writeCommented(file(), "Source : ");
    file() << sourceTypeNames_[source_] << " " << sourceName_ << endl;
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

Foam::fieldValues::cellSource::cellSource
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    fieldValue(name, obr, dict, typeName, loadFromFiles),
    source_(sourceTypeNames_.read(dict.lookup("source"))),
    operation_(operationTypeNames_.read(dict.lookup("operation"))),
    nCells_(0),
    cellId_(),
    weightFieldName_("none"),
    writeVolume_(dict.lookupOrDefault("writeVolume", false))
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fieldValues::cellSource::~cellSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fieldValues::cellSource::read(const dictionary& dict)
{
    fieldValue::read(dict);

    if (active_)
    {
        // no additional info to read
        initialise(dict);
    }
}


void Foam::fieldValues::cellSource::write()
{
    fieldValue::write();

    if (active_)
    {
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
            if (log_) Info<< "    total volume = " << volume_ << endl;
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

        if (log_) Info<< endl;
    }
}


// ************************************************************************* //

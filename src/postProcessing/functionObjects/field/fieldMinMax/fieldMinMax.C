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

#include "fieldMinMax.H"
#include "fieldTypes.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fieldMinMax, 0);

    template<>
    const char* NamedEnum
    <
        fieldMinMax::modeType,
        2
    >::names[] =
    {
        "magnitude",
        "component"
    };
}


const Foam::NamedEnum<Foam::fieldMinMax::modeType, 2>
Foam::fieldMinMax::modeTypeNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fieldMinMax::fieldMinMax
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    functionObjectFile(obr, name, typeName),
    name_(name),
    obr_(obr),
    active_(true),
    log_(true),
    location_(true),
    mode_(mdMag),
    fieldSet_()
{
    // Check if the available mesh is an fvMesh otherise deactivate
    if (!isA<fvMesh>(obr_))
    {
        active_ = false;
        WarningIn
        (
            "fieldMinMax::fieldMinMax"
            "("
                "const word&, "
                "const objectRegistry&, "
                "const dictionary&, "
                "const bool"
            ")"
        )   << "No fvMesh available, deactivating " << name_
            << endl;
    }

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fieldMinMax::~fieldMinMax()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fieldMinMax::read(const dictionary& dict)
{
    if (active_)
    {
        log_ = dict.lookupOrDefault<Switch>("log", true);
        location_ = dict.lookupOrDefault<Switch>("location", true);

        mode_ = modeTypeNames_[dict.lookupOrDefault<word>("mode", "magnitude")];
        dict.lookup("fields") >> fieldSet_;
    }
}


void Foam::fieldMinMax::writeFileHeader(const label i)
{
    OFstream& file = this->file();

    writeHeader(file, "Field minima and maxima");
    writeCommented(file, "Time");

    if (location_)
    {
        writeTabbed(file, "field");

        writeTabbed(file, "min");
        writeTabbed(file, "location(min)");

        if (Pstream::parRun())
        {
            writeTabbed(file, "processor");
        }

        writeTabbed(file, "max");
        writeTabbed(file, "location(max)");

        if (Pstream::parRun())
        {
            writeTabbed(file, "processor");
        }
    }
    else
    {
        forAll(fieldSet_, fieldI)
        {
            writeTabbed(file, "min(" + fieldSet_[fieldI] + ')');
            writeTabbed(file, "max(" + fieldSet_[fieldI] + ')');
        }
    }

    file<< endl;
}


void Foam::fieldMinMax::execute()
{
    // Do nothing - only valid on write
}


void Foam::fieldMinMax::end()
{
    // Do nothing - only valid on write
}


void Foam::fieldMinMax::timeSet()
{
    // Do nothing - only valid on write
}


void Foam::fieldMinMax::write()
{
    if (active_)
    {
        functionObjectFile::write();

        if (!location_)
        {
            file()<< obr_.time().value();
        }

        Info(log_)<< type() << " " << name_ <<  " output:" << nl;

        forAll(fieldSet_, fieldI)
        {
            calcMinMaxFields<scalar>(fieldSet_[fieldI], mdCmpt);
            calcMinMaxFields<vector>(fieldSet_[fieldI], mode_);
            calcMinMaxFields<sphericalTensor>(fieldSet_[fieldI], mode_);
            calcMinMaxFields<symmTensor>(fieldSet_[fieldI], mode_);
            calcMinMaxFields<tensor>(fieldSet_[fieldI], mode_);
        }

        if (!location_)
        {
            file()<< endl;
        }

        Info(log_)<< endl;
    }
}



// ************************************************************************* //

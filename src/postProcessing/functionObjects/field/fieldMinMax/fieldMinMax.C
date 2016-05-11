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

#include "fieldMinMax.H"
#include "fieldTypes.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fieldMinMax, 0);
}
}

template<>
const char* Foam::NamedEnum
<
    Foam::functionObjects::fieldMinMax::modeType,
    2
>::names[] = {"magnitude", "component"};

const Foam::NamedEnum
<
    Foam::functionObjects::fieldMinMax::modeType,
    2
> Foam::functionObjects::fieldMinMax::modeTypeNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldMinMax::fieldMinMax
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    functionObjectFiles(obr, name, typeName),
    name_(name),
    obr_(obr),
    log_(true),
    location_(true),
    mode_(mdMag),
    fieldSet_()
{
    if (!isA<fvMesh>(obr))
    {
        FatalErrorInFunction
            << "objectRegistry is not an fvMesh" << exit(FatalError);
    }

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fieldMinMax::~fieldMinMax()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::functionObjects::fieldMinMax::read(const dictionary& dict)
{
    log_ = dict.lookupOrDefault<Switch>("log", true);
    location_ = dict.lookupOrDefault<Switch>("location", true);

    mode_ = modeTypeNames_[dict.lookupOrDefault<word>("mode", "magnitude")];
    dict.lookup("fields") >> fieldSet_;
}


void Foam::functionObjects::fieldMinMax::writeFileHeader(const label i)
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
        forAll(fieldSet_, fieldi)
        {
            writeTabbed(file, "min(" + fieldSet_[fieldi] + ')');
            writeTabbed(file, "max(" + fieldSet_[fieldi] + ')');
        }
    }

    file<< endl;
}


void Foam::functionObjects::fieldMinMax::execute()
{}


void Foam::functionObjects::fieldMinMax::end()
{}


void Foam::functionObjects::fieldMinMax::timeSet()
{}


void Foam::functionObjects::fieldMinMax::write()
{
    functionObjectFiles::write();

    if (!location_) writeTime(file());
    if (log_) Info<< type() << " " << name_ <<  " output:" << nl;

    forAll(fieldSet_, fieldi)
    {
        calcMinMaxFields<scalar>(fieldSet_[fieldi], mdCmpt);
        calcMinMaxFields<vector>(fieldSet_[fieldi], mode_);
        calcMinMaxFields<sphericalTensor>(fieldSet_[fieldi], mode_);
        calcMinMaxFields<symmTensor>(fieldSet_[fieldi], mode_);
        calcMinMaxFields<tensor>(fieldSet_[fieldi], mode_);
    }

    if (!location_) file()<< endl;
    if (log_) Info<< endl;
}



// ************************************************************************* //

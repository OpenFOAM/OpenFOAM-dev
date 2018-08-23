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

#include "fieldMinMax.H"
#include "fieldTypes.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fieldMinMax, 0);
    addToRunTimeSelectionTable(functionObject, fieldMinMax, dictionary);
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


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldMinMax::fieldMinMax
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    location_(true),
    mode_(modeType::mag),
    fieldSet_()
{
    read(dict);
    resetName(typeName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fieldMinMax::~fieldMinMax()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fieldMinMax::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    location_ = dict.lookupOrDefault<Switch>("location", true);

    mode_ = modeTypeNames_[dict.lookupOrDefault<word>("mode", "magnitude")];
    dict.lookup("fields") >> fieldSet_;

    return true;
}


bool Foam::functionObjects::fieldMinMax::execute()
{
    return true;
}


bool Foam::functionObjects::fieldMinMax::write()
{
    logFiles::write();

    if (Pstream::master() && !location_) writeTime(file());
    Log << type() << " " << name() <<  " write:" << nl;

    forAll(fieldSet_, fieldi)
    {
        calcMinMaxFields<scalar>(fieldSet_[fieldi], modeType::cmpt);
        calcMinMaxFields<vector>(fieldSet_[fieldi], mode_);
        calcMinMaxFields<sphericalTensor>(fieldSet_[fieldi], mode_);
        calcMinMaxFields<symmTensor>(fieldSet_[fieldi], mode_);
        calcMinMaxFields<tensor>(fieldSet_[fieldi], mode_);
    }

    if (Pstream::master() && !location_) file() << endl;
    Log << endl;

    return true;
}



// ************************************************************************* //

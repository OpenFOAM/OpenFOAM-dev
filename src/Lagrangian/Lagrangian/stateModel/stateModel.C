/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "stateModel.H"
#include "Time.H"
#include "timeIOdictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(stateModel::stateModelWriter, 0);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::IOobject Foam::stateModel::io
(
    const word& name,
    const objectRegistry& db,
    const bool processor,
    const bool write
)
{
    return
        IOobject
        (
            name + "State",
            db.time().name(),
            processor ? "processorUniform" : "uniform",
            db,
            write ? IOobject::NO_READ : IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::stateModel::stateModelWriter::stateModelWriter
(
    const stateModel& model,
    const bool processor
)
:
    regIOobject(io(model.name(), model.db(), processor, true)),
    model_(model),
    processor_(processor)
{}


Foam::stateModel::stateModel()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::stateModel::~stateModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::stateModel::stateModelWriter::writeData(Ostream& os) const
{
    if (!processor_)
    {
        model_.writeState(os);
    }

    if (processor_ == Pstream::parRun())
    {
        model_.writeProcessorState(os);
    }

    return os.good();
}


bool Foam::stateModel::writeState(const bool write) const
{
    if (write)
    {
        bool allOk = true;

        const bool ok = stateModel::stateModelWriter(*this, false).write();
        allOk = allOk && ok;

        if (Pstream::parRun())
        {
            const bool ok = stateModel::stateModelWriter(*this, true).write();
            allOk = allOk && ok;
        }

        return allOk;
    }
    else
    {
        return false;
    }
}


void Foam::stateModel::writeState(Ostream& os) const
{}


void Foam::stateModel::writeProcessorState(Ostream& os) const
{}


Foam::dictionary Foam::stateModel::stateDict
(
    const word& name,
    const objectRegistry& db
)
{
    timeIOdictionary stateDict(io(name, db, false, false));

    if (Pstream::parRun())
    {
        stateDict.merge(timeIOdictionary(io(name, db, true, false)));
    }

    return stateDict;
}


// ************************************************************************* //

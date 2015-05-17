/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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

#include "residuals.H"
#include "volFields.H"
#include "dictionary.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(residuals, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::residuals::residuals
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
    fieldSet_()
{
    // Check if the available mesh is an fvMesh otherwise deactivate
    if (!isA<fvMesh>(obr_))
    {
        active_ = false;
        WarningIn
        (
            "residuals::residuals"
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

Foam::residuals::~residuals()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::residuals::read(const dictionary& dict)
{
    if (active_)
    {
        dict.lookup("fields") >> fieldSet_;
    }
}


void Foam::residuals::writeFileHeader(const label i)
{
    if (Pstream::master())
    {
        writeHeader(file(), "Residuals");
        writeCommented(file(), "Time");

        forAll(fieldSet_, fieldI)
        {
            writeTabbed(file(), fieldSet_[fieldI]);
        }

        file() << endl;
    }
}


void Foam::residuals::execute()
{
    // Do nothing - only valid on write
}


void Foam::residuals::end()
{
    // Do nothing - only valid on write
}


void Foam::residuals::timeSet()
{
    // Do nothing - only valid on write
}


void Foam::residuals::write()
{
    if (active_)
    {
        functionObjectFile::write();

        if (Pstream::master())
        {
            file()<< obr_.time().value();

            forAll(fieldSet_, fieldI)
            {
                const word& fieldName = fieldSet_[fieldI];

                writeResidual<scalar>(fieldName);
                writeResidual<vector>(fieldName);
                writeResidual<sphericalTensor>(fieldName);
                writeResidual<symmTensor>(fieldName);
                writeResidual<tensor>(fieldName);
            }

            file() << endl;
        }
    }
}

// ************************************************************************* //

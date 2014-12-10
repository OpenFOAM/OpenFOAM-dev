/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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

#include "calcMag.H"
#include "volFields.H"
#include "dictionary.H"
#include "calcMag.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(calcMag, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::calcMag::calcMag
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    active_(true),
    fieldName_("undefined-fieldName"),
    resultName_("undefined-resultName")
{
    // Check if the available mesh is an fvMesh, otherwise deactivate
    if (!isA<fvMesh>(obr_))
    {
        active_ = false;
        WarningIn
        (
            "calcMag::calcMag"
            "("
                "const word&, "
                "const objectRegistry&, "
                "const dictionary&, "
                "const bool"
            ")"
        )   << "No fvMesh available, deactivating." << nl
            << endl;
    }

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::calcMag::~calcMag()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::calcMag::read(const dictionary& dict)
{
    if (active_)
    {
        dict.lookup("fieldName") >> fieldName_;
        dict.lookup("resultName") >> resultName_;

        if (resultName_ == "none")
        {
            resultName_ = "mag(" + fieldName_ + ")";
        }
    }
}


void Foam::calcMag::execute()
{
    if (active_)
    {
        bool processed = false;

        calc<scalar>(fieldName_, resultName_, processed);
        calc<vector>(fieldName_, resultName_, processed);
        calc<sphericalTensor>(fieldName_, resultName_, processed);
        calc<symmTensor>(fieldName_, resultName_, processed);
        calc<tensor>(fieldName_, resultName_, processed);

        if (!processed)
        {
            WarningIn("void Foam::calcMag::write()")
                << "Unprocessed field " << fieldName_ << endl;
        }
    }
}


void Foam::calcMag::end()
{
    if (active_)
    {
        execute();
    }
}


void Foam::calcMag::timeSet()
{
    // Do nothing
}


void Foam::calcMag::write()
{
    if (active_)
    {
        if (obr_.foundObject<regIOobject>(resultName_))
        {
            const regIOobject& field =
                obr_.lookupObject<regIOobject>(resultName_);

            Info<< type() << " " << name_ << " output:" << nl
                << "    writing field " << field.name() << nl << endl;

            field.write();
        }
    }
}


// ************************************************************************* //

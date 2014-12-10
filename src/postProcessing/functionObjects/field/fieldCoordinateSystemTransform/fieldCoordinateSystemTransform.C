/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "fieldCoordinateSystemTransform.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(fieldCoordinateSystemTransform, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fieldCoordinateSystemTransform::fieldCoordinateSystemTransform
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
    fieldSet_(),
    coordSys_(obr, dict)
{
    // Check if the available mesh is an fvMesh otherise deactivate
    if (isA<fvMesh>(obr_))
    {
        read(dict);

        Info<< type() << " " << name_ << ":" << nl
            << "   Applying transformation from global Cartesian to local "
            << coordSys_ << nl << endl;
    }
    else
    {
        active_ = false;
        WarningIn
        (
            "fieldCoordinateSystemTransform::fieldCoordinateSystemTransform"
            "("
                "const word&, "
                "const objectRegistry&, "
                "const dictionary&, "
                "const bool"
            ")"
        )   << "No fvMesh available, deactivating " << name_
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fieldCoordinateSystemTransform::~fieldCoordinateSystemTransform()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fieldCoordinateSystemTransform::read(const dictionary& dict)
{
    if (active_)
    {
        dict.lookup("fields") >> fieldSet_;
    }
}


void Foam::fieldCoordinateSystemTransform::execute()
{
    if (active_)
    {
        Info<< type() << " " << name_ << " output:" << nl;

        forAll(fieldSet_, fieldI)
        {
            // If necessary load field
            transform<scalar>(fieldSet_[fieldI]);
            transform<vector>(fieldSet_[fieldI]);
            transform<sphericalTensor>(fieldSet_[fieldI]);
            transform<symmTensor>(fieldSet_[fieldI]);
            transform<tensor>(fieldSet_[fieldI]);
        }
    }
}


void Foam::fieldCoordinateSystemTransform::end()
{
    if (active_)
    {
        execute();
    }
}


void Foam::fieldCoordinateSystemTransform::timeSet()
{
    // Do nothing
}


void Foam::fieldCoordinateSystemTransform::write()
{
    if (active_)
    {
        Info<< type() << " " << name_ << " output:" << nl;

        forAll(fieldSet_, fieldI)
        {
            const word fieldName = fieldSet_[fieldI] + ":Transformed";

            const regIOobject& field =
                obr_.lookupObject<regIOobject>(fieldName);

            Info<< "    writing field " << field.name() << nl;

            field.write();
        }

        Info<< endl;
    }
}


// ************************************************************************* //

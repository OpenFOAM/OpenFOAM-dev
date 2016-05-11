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

#include "fieldCoordinateSystemTransform.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fieldCoordinateSystemTransform, 0);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldCoordinateSystemTransform::
fieldCoordinateSystemTransform
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    fieldSet_(),
    coordSys_(obr, dict)
{
    if (!isA<fvMesh>(obr))
    {
        FatalErrorInFunction
            << "objectRegistry is not an fvMesh" << exit(FatalError);
    }

    read(dict);

    Info<< type() << " " << name_ << ":" << nl
        << "   Applying transformation from global Cartesian to local "
        << coordSys_ << nl << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fieldCoordinateSystemTransform::
~fieldCoordinateSystemTransform()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::functionObjects::fieldCoordinateSystemTransform::read
(
    const dictionary& dict
)
{
    dict.lookup("fields") >> fieldSet_;
}


void Foam::functionObjects::fieldCoordinateSystemTransform::execute()
{
    Info<< type() << " " << name_ << " output:" << nl;

    forAll(fieldSet_, fieldi)
    {
        // If necessary load field
        transform<scalar>(fieldSet_[fieldi]);
        transform<vector>(fieldSet_[fieldi]);
        transform<sphericalTensor>(fieldSet_[fieldi]);
        transform<symmTensor>(fieldSet_[fieldi]);
        transform<tensor>(fieldSet_[fieldi]);
    }
}


void Foam::functionObjects::fieldCoordinateSystemTransform::end()
{
    execute();
}


void Foam::functionObjects::fieldCoordinateSystemTransform::timeSet()
{}


void Foam::functionObjects::fieldCoordinateSystemTransform::write()
{
    Info<< type() << " " << name_ << " output:" << nl;

    forAll(fieldSet_, fieldi)
    {
        const word fieldName = fieldSet_[fieldi] + ":Transformed";

        const regIOobject& field =
            obr_.lookupObject<regIOobject>(fieldName);

        Info<< "    writing field " << field.name() << nl;

        field.write();
    }

    Info<< endl;
}


// ************************************************************************* //

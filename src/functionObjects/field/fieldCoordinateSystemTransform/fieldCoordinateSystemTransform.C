/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fieldCoordinateSystemTransform, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        fieldCoordinateSystemTransform,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldCoordinateSystemTransform::
fieldCoordinateSystemTransform
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fields_(),
    coordSys_(coordinateSystem::New(mesh_, dict)())
{
    read(dict);

    Log << type() << " " << name << ":" << nl
        << "   Applying transformation from global Cartesian to local "
        << coordSys_ << nl << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fieldCoordinateSystemTransform::
~fieldCoordinateSystemTransform()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::word
Foam::functionObjects::fieldCoordinateSystemTransform::transformFieldName
(
    const word& fieldName
) const
{
    return fieldName + ":Transformed";
}


bool Foam::functionObjects::fieldCoordinateSystemTransform::read
(
    const dictionary& dict
)
{
    fvMeshFunctionObject::read(dict);

    dict.lookup("fields") >> fields_;

    return true;
}


Foam::wordList
Foam::functionObjects::fieldCoordinateSystemTransform::fields() const
{
    return fields_;
}


bool Foam::functionObjects::fieldCoordinateSystemTransform::execute()
{
    forAll(fields_, fieldi)
    {
        transform<scalar>(fields_[fieldi]);
        transform<vector>(fields_[fieldi]);
        transform<sphericalTensor>(fields_[fieldi]);
        transform<symmTensor>(fields_[fieldi]);
        transform<tensor>(fields_[fieldi]);
    }

    return true;
}


bool Foam::functionObjects::fieldCoordinateSystemTransform::write()
{
    forAll(fields_, fieldi)
    {
        writeObject(transformFieldName(fields_[fieldi]));
    }

    return true;
}


// ************************************************************************* //

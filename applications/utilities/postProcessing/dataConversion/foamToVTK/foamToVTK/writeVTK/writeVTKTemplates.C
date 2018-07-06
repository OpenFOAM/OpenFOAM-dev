/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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

#include "writeVTK.H"
#include "objectRegistry.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class GeoField>
Foam::UPtrList<const GeoField>
Foam::functionObjects::writeVTK::lookupFields() const
{
    DynamicList<word> allNames(obr_.toc().size());
    forAll(objectNames_, i)
    {
        wordList names(obr_.names<GeoField>(objectNames_[i]));

        if (names.size())
        {
            allNames.append(names);
        }
    }

    UPtrList<const GeoField> fields(allNames.size());

    forAll(allNames, i)
    {
        const GeoField& field = obr_.lookupObject<GeoField>(allNames[i]);
        Info<< "    Writing " << GeoField::typeName
            << " field " << field.name() << endl;
        fields.set(i, &field);
    }

    return fields;
}


// ************************************************************************* //

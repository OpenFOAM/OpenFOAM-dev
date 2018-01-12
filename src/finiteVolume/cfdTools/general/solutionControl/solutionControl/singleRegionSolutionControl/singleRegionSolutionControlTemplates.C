/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::singleRegionSolutionControl::storePrevIterTypeFields() const
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    HashTable<fieldType*>
        flds(mesh_.objectRegistry::lookupClass<fieldType>());

    forAllIter(typename HashTable<fieldType*>, flds, iter)
    {
        fieldType& fld = *iter();

        const word& fName = fld.name();

        size_t prevIterField = fName.find("PrevIter");

        if (prevIterField == word::npos && mesh_.relaxField(fName))
        {
            fld.storePrevIter();
        }
    }
}


// ************************************************************************* //

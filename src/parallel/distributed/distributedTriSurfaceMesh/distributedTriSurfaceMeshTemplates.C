/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "distributedTriSurfaceMesh.H"
#include "triSurfaceFields.H"
#include "distributionMap.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//template<class Type>
//void Foam::distributedTriSurfaceMesh::getField
//(
//    const word& fieldName,
//    const List<pointIndexHit>& info,
//    List<Type>& values
//) const
//{
//    typedef DimensionedField<Type, triSurfaceGeoMesh> DimensionedSurfField;
//
//
//    // Get query data (= local index of triangle)
//    // ~~~~~~~~~~~~~~
//
//    labelList triangleIndex(info.size());
//    autoPtr<distributionMap> mapPtr
//    (
//        calcLocalQueries
//        (
//            info,
//            triangleIndex
//        )
//    );
//    const distributionMap& map = mapPtr();
//
//
//    // Do my tests
//    // ~~~~~~~~~~~
//
//    const DimensionedSurfField& fld = lookupObject<DimensionedSurfField>
//    (
//        fieldName
//    );
//    const triSurface& s = static_cast<const triSurface&>(*this);
//
//    values.setSize(triangleIndex.size());
//
//    forAll(triangleIndex, i)
//    {
//        label triI = triangleIndex[i];
//        values[i] = fld[triI];
//    }
//
//
//    // Send back results
//    // ~~~~~~~~~~~~~~~~~
//
//    map.reverseDistribute(info.size(), values);
//}


template<class Type>
void Foam::distributedTriSurfaceMesh::distributeFields
(
    const distributionMap& map
)
{
    typedef DimensionedField<Type, triSurfaceGeoMesh> DimensionedSurfField;

    HashTable<DimensionedSurfField*> fields
    (
        objectRegistry::lookupClass<DimensionedSurfField>()
    );

    for
    (
        typename HashTable<DimensionedSurfField*>::iterator fieldIter =
            fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        DimensionedSurfField& field = *fieldIter();

        label oldSize = field.size();

        map.distribute(field);

        if (debug)
        {
            Info<< "Mapped " << field.typeName << ' ' << field.name()
                << " from size " << oldSize << " to size " << field.size()
                << endl;
        }
    }
}


// ************************************************************************* //

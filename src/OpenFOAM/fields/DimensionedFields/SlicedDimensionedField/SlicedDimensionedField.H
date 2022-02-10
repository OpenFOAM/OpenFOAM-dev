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

Class
    Foam::SlicedDimensionedField

Description
    Specialisation of DimensionedField which holds a slice of a given complete
    field in such a form that it acts as a DimensionedField.

    The destructor is wrapped to avoid deallocation of the storage of the
    complete fields when this is destroyed.

SourceFiles
    SlicedDimensionedField.C

\*---------------------------------------------------------------------------*/

#ifndef SlicedDimensionedField_H
#define SlicedDimensionedField_H

#include "DimensionedField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class SlicedDimensionedField Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class GeoMesh>
class SlicedDimensionedField
:
    public DimensionedField<Type, GeoMesh>
{
public:

    // Public Typedefs

        //- Type of mesh on which this SlicedDimensionedField is instantiated
        typedef typename GeoMesh::Mesh Mesh;


    // Constructors

        //- Construct from components and field to slice
        SlicedDimensionedField
        (
            const IOobject&,
            const Mesh&,
            const dimensionSet&,
            const Field<Type>& iField
        );


    //- Destructor
    ~SlicedDimensionedField();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SlicedDimensionedField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

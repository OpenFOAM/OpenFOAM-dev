/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2026 OpenFOAM Foundation
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
    Foam::GeometricFieldListSlicer

Description
    Class to provide list slices to different parts of a geometric field

See also
    Foam::FieldListSlice

\*---------------------------------------------------------------------------*/

#ifndef GeometricFieldListSlicer_H
#define GeometricFieldListSlicer_H

#include "DimensionedFieldListSlicer.H"
#include "GeometricField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                  Class GeometricFieldListSlicer Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class GeometricFieldListSlicer
:
    public DimensionedFieldListSlicer<Type>
{
    // Private Data

        //- Slices into the patch field data. Indexed by patch, then by field.
        List<List<UList<Type>>> patchFields_;


public:

    // Constructors

        //- Construct null. Set later.
        inline GeometricFieldListSlicer()
        {}

        //- Construct from a mesh and a list of fields
        template<class GeoMesh, class GeoField>
        inline GeometricFieldListSlicer
        (
            const GeoMesh& mesh,
            const UPtrList<GeoField>& geoFields
        )
        {
            set(mesh, geoFields);
        }


    // Member Functions

        //- Set the field pointers
        template<class GeoMesh, class GeoField>
        inline void set
        (
            const GeoMesh& mesh,
            const UPtrList<GeoField>& geoFields
        )
        {
            DimensionedFieldListSlicer<Type>::set(geoFields);

            List<List<UList<Type>>> tpatchFields(mesh.boundary().size());
            patchFields_.transfer(tpatchFields);

            forAll(mesh.boundary(), patchi)
            {
                FieldListSlicer<Type>::set
                (
                    [&](const label geoFieldi) -> decltype(auto)
                    {
                        return geoFields[geoFieldi].boundaryField()[patchi];
                    },
                    geoFields.size(),
                    patchFields_[patchi]
                );
            }
        }

        //- Return a slice for a patch element
        inline FieldListSlice<Type> patchSlice
        (
            const label patchi,
            const label patchElementi
        ) const
        {
            return FieldListSlice<Type>(patchFields_[patchi], patchElementi);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

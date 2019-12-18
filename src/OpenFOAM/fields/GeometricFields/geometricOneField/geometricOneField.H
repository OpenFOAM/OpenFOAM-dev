/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::geometricOneField

Description
    A class representing the concept of a GeometricField of 1 used to avoid
    unnecessary manipulations for objects which are known to be one at
    compile-time.

    Used for example as the density argument to a function written for
    compressible to be used for incompressible flow.

\*---------------------------------------------------------------------------*/

#ifndef geometricOneField_H
#define geometricOneField_H

#include "oneFieldField.H"
#include "dimensionSet.H"
#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class geometricOneField Declaration
\*---------------------------------------------------------------------------*/

class geometricOneField
:
    public one
{

public:

    // Public Typedefs

        typedef oneField Internal;
        typedef oneField Patch;
        typedef oneFieldField Boundary;
        typedef one cmptType;

    // Constructors

        //- Construct null
        geometricOneField()
        {}


    // Member Operators

        inline const dimensionSet& dimensions() const;

        inline one operator[](const label) const;

        inline oneField field() const;

        inline oneField oldTime() const;

        inline Internal operator()() const;

        inline Internal v() const;

        inline typename Internal::FieldType primitiveField() const;

        inline Boundary boundaryField() const;
};


inline const geometricOneField& operator*
(
    const geometricOneField&,
    const geometricOneField&
);

inline const geometricOneField::Internal& operator*
(
    const geometricOneField::Internal&,
    const geometricOneField&
);

inline const geometricOneField::Internal& operator*
(
    const geometricOneField&,
    const geometricOneField::Internal&
);

inline const geometricOneField& operator/
(
    const geometricOneField&,
    const geometricOneField&
);

inline const geometricOneField::Internal& operator/
(
    const geometricOneField::Internal&,
    const geometricOneField&
);

inline const geometricOneField::Internal& operator/
(
    const geometricOneField&,
    const geometricOneField::Internal&
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "geometricOneFieldI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

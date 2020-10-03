/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
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
    Foam::cyclicACMIPointPatchField

Description
    Cyclic ACMI front and back plane patch field

SourceFiles
    cyclicACMIPointPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicACMIPointPatchField_H
#define cyclicACMIPointPatchField_H

#include "cyclicAMIPointPatchField.H"
#include "cyclicACMIPointPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class cyclicACMIPointPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class cyclicACMIPointPatchField
:
    public cyclicAMIPointPatchField<Type>
{
public:

    //- Runtime type information
    TypeName(cyclicACMIPointPatch::typeName_());


    // Constructors

        //- Inherit parent constructors
        using cyclicAMIPointPatchField<Type>::cyclicAMIPointPatchField;

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<Type>> clone
        (
            const DimensionedField<Type, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchField<Type>>
            (
                new cyclicACMIPointPatchField<Type>
                (
                    *this, iF
                )
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2019 OpenFOAM Foundation
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
    Foam::cyclicRepeatAMIFvPatchField

Description
    This boundary condition enforces an repeating condition between a pair of
    boundaries, whereby communication between the patches is performed using
    an arbitrary mesh interface (AMI) interpolation.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            cyclicRepeatAMI;
    }
    \endverbatim

Note
    The outer boundary of the patch pairs must be similar, i.e. if the owner
    patch is transformed to the neighbour patch, the outer perimeter of each
    patch should be identical (or very similar).

See also
    Foam::AMIInterpolation

SourceFiles
    cyclicRepeatAMIFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicRepeatAMIFvPatchField_H
#define cyclicRepeatAMIFvPatchField_H

#include "cyclicAMIFvPatchField.H"
#include "cyclicRepeatAMIFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class cyclicRepeatAMIFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class cyclicRepeatAMIFvPatchField
:
    public cyclicAMIFvPatchField<Type>
{
public:

    //- Runtime type information
    TypeName(cyclicRepeatAMIFvPatch::typeName_());


    // Constructors

        //- Inherit parent constructors
        using cyclicAMIFvPatchField<Type>::cyclicAMIFvPatchField;

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new cyclicRepeatAMIFvPatchField<Type>(*this)
            );
        }

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new cyclicRepeatAMIFvPatchField<Type>(*this, iF)
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

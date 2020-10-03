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
    Foam::cyclicACMIFvPatchField

Description
    This boundary condition enforces a cyclic condition between a pair of
    boundaries, whereby communication between the patches is performed using
    an arbitrarily coupled mesh interface (ACMI) interpolation.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            cyclicACMI;
    }
    \endverbatim

See also
    Foam::AMIInterpolation

SourceFiles
    cyclicACMIFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicACMIFvPatchField_H
#define cyclicACMIFvPatchField_H

#include "cyclicAMIFvPatchField.H"
#include "cyclicACMIFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class cyclicACMIFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class cyclicACMIFvPatchField
:
    public cyclicAMIFvPatchField<Type>
{
public:

    //- Runtime type information
    TypeName(cyclicACMIFvPatch::typeName_());


    // Constructors

        //- Inherit parent constructors
        using cyclicAMIFvPatchField<Type>::cyclicAMIFvPatchField;

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new cyclicACMIFvPatchField<Type>(*this, iF)
            );
        }


    // Member Functions

        // Access

            //- Return the fv patch
            const cyclicACMIFvPatch& cyclicACMIPatch() const
            {
                return refCast<const cyclicACMIFvPatch>(this->patch());
            }


        // Evaluation functions

            //- Return true if this patch field fixes a value
            //  Needed to check if a level has to be specified while solving
            //  Poissons equations
            virtual bool fixesValue() const
            {
                const scalarField& mask =
                    cyclicACMIPatch().cyclicACMIPatch().mask();

                if (gMax(mask) > 1e-5)
                {
                    // regions connected
                    return false;
                }
                else
                {
                    // fully separated
                    return nonOverlapPatchField().fixesValue();
                }
            }

            //- Return reference to non-overlapping patchField
            const fvPatchField<Type>& nonOverlapPatchField() const;

            //- Manipulate matrix
            virtual void manipulateMatrix(fvMatrix<Type>& matrix);

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "cyclicACMIFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

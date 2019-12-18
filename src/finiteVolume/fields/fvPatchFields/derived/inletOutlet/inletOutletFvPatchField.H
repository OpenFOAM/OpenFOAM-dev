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
    Foam::inletOutletFvPatchField

Description
    This boundary condition provides a generic outflow condition, with
    specified inflow for the case of return flow.

Usage
    \table
        Property     | Description             | Required    | Default value
        phi          | Flux field name         | no          | phi
        inletValue   | Inlet value for reverse flow | yes    |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            inletOutlet;
        phi             phi;
        inletValue      uniform 0;
        value           uniform 0;
    }
    \endverbatim

    The mode of operation is determined by the sign of the flux across the
    patch faces.

Note
    Sign conventions:
    - Positive flux (out of domain): apply zero-gradient condition
    - Negative flux (into of domain): apply the "inletValue" fixed-value

See also
    Foam::mixedFvPatchField
    Foam::zeroGradientFvPatchField
    Foam::outletInletFvPatchField

SourceFiles
    inletOutletFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef inletOutletFvPatchField_H
#define inletOutletFvPatchField_H

#include "mixedFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class inletOutletFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class inletOutletFvPatchField
:
    public mixedFvPatchField<Type>
{

protected:

    // Protected data

        //- Name of flux field
        word phiName_;


public:

    //- Runtime type information
    TypeName("inletOutlet");


    // Constructors

        //- Construct from patch and internal field
        inletOutletFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        inletOutletFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given inletOutletFvPatchField onto a new patch
        inletOutletFvPatchField
        (
            const inletOutletFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        inletOutletFvPatchField
        (
            const inletOutletFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new inletOutletFvPatchField<Type>(*this)
            );
        }

        //- Copy constructor setting internal field reference
        inletOutletFvPatchField
        (
            const inletOutletFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new inletOutletFvPatchField<Type>(*this, iF)
            );
        }


    // Member Functions

        // Attributes

            //- Return true: this patch field is altered by assignment
            virtual bool assignable() const
            {
                return true;
            }


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member Operators

        virtual void operator=(const fvPatchField<Type>& pvf);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "inletOutletFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::zeroInletOutletFvPatchField

Description
    This boundary condition provides a generic outflow condition

    Zero gradient is applied for outflow and a zero value applied for reverse
    flow.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            zeroInletOutlet;
    }
    \endverbatim

See also
    Foam::inletOutletFvPatchField
    Foam::zeroGradientFvPatchField

SourceFiles
    zeroInletOutletFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef zeroInletOutletFvPatchField_H
#define zeroInletOutletFvPatchField_H

#include "inletOutletFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class zeroInletOutletFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class zeroInletOutletFvPatchField
:
    public inletOutletFvPatchField<Type>
{

public:

    //- Runtime type information
    TypeName("zeroInletOutlet");


    // Constructors

        using inletOutletFvPatchField<Type>::inletOutletFvPatchField;

        //- Construct from patch, internal field and dictionary
        zeroInletOutletFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new zeroInletOutletFvPatchField<Type>(*this, iF)
            );
        }


    // Member Functions

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "zeroInletOutletFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

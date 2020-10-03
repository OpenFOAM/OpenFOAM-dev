/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2020 OpenFOAM Foundation
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
    Foam::waveInletOutletFvPatchField

Description
    This boundary condition provides an inlet-outlet condition with differing
    inlet values on either side of a wave interface. All the wave modelling
    parameters are obtained from a centrally registered waveSuperposition
    class.

Usage
    \table
        Property        | Description                | Req'd? | Default
        phi             | Name of the flux field     | no     | phi
        inletValueAbove | inlet value above the wave | no     | None
        inletValueBelow | inlet value below the wave | no     | None
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type        waveInletOutlet;
        phi         phi;
        inletValueAbove 0.01;
        inletValueBelow table ((0 0.01) (10 0.1));
    }
    \endverbatim

See also
    Foam::waveSuperposition

SourceFiles
    waveInletOutletFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef waveInletOutletFvPatchField_H
#define waveInletOutletFvPatchField_H

#include "mixedFvPatchField.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class waveInletOutletFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class waveInletOutletFvPatchField
:
    public mixedFvPatchField<Type>
{
    // Private Data

        //- Inlet value above the wave
        const autoPtr<Function1<Type>> inletValueAbove_;

        //- Inlet value below the wave
        const autoPtr<Function1<Type>> inletValueBelow_;

        //- Name of flux field
        const word phiName_;


public:

    //- Runtime type information
    TypeName("waveInletOutlet");


    // Constructors

        //- Construct from patch and internal field
        waveInletOutletFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        waveInletOutletFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given field onto a new patch
        waveInletOutletFvPatchField
        (
            const waveInletOutletFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        waveInletOutletFvPatchField
        (
            const waveInletOutletFvPatchField<Type>&
        ) = delete;

        //- Copy constructor setting internal field reference
        waveInletOutletFvPatchField
        (
            const waveInletOutletFvPatchField<Type>&,
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
                new waveInletOutletFvPatchField<Type>
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "waveInletOutletFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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
    Foam::advectiveFvPatchField

Description
    This boundary condition provides an advective outflow condition, based on
    solving DDt(W, field) = 0 at the boundary where \c W is the wave velocity
    and \c field is the field to which this boundary condition is applied.

    The standard (Euler, backward, CrankNicolson, localEuler) time schemes are
    supported.  Additionally an optional mechanism to relax the value at
    the boundary to a specified far-field value is provided which is
    switched on by specifying the relaxation length-scale \c lInf and the
    optionally time-varying far-field value \c fieldInf.

    The flow/wave speed \c (w) at the outlet is provided by the virtual function
    advectionSpeed() the default implementation of which requires the name of
    the flux field \c (phi) and optionally the density \c (rho) if the
    mass-flux rather than the volumetric-flux is given.

    The flow/wave speed at the outlet can be changed by deriving a specialised
    BC from this class and over-riding advectionSpeed()  e.g. in
    waveTransmissiveFvPatchField the advectionSpeed() calculates and returns
    the flow-speed plus the acoustic wave speed creating an acoustic wave
    transmissive boundary condition.

Usage
    \table
        Property     | Description              | Required    | Default value
        phi          | Flux field name          | no          | phi
        rho          | Density field name       | no          | rho
        fieldInf     | Far field value Function1           | no |
        lInf         | Far field value relaxation distance | no |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            advective;
        phi             phi;
    }
    \endverbatim

    Note:
        If \c lInf is specified, \c fieldInf is required.
        \c rho is only required in the case of a mass-based flux.

SourceFiles
    advectiveFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef advectiveFvPatchField_H
#define advectiveFvPatchField_H

#include "mixedFvPatchField.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class advectiveFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class advectiveFvPatchField
:
    public mixedFvPatchField<Type>
{
protected:

    // Private Data

        //- Name of the flux transporting the field
        word phiName_;

        //- Name of the density field used to normalise the mass flux
        //  if necessary
        word rhoName_;

        //- Far-field value
        autoPtr<Function1<Type>> fieldInf_;

        //- Relaxation length-scale
        scalar lInf_;


public:

    //- Runtime type information
    TypeName("advective");


    // Constructors

        //- Construct from patch, internal field and dictionary
        advectiveFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, fvMesh>&,
            const dictionary&
        );

        //- Construct by mapping given advectiveFvPatchField
        //  onto a new patch
        advectiveFvPatchField
        (
            const advectiveFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, fvMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        advectiveFvPatchField(const advectiveFvPatchField&) = delete;

        //- Copy constructor setting internal field reference
        advectiveFvPatchField
        (
            const advectiveFvPatchField&,
            const DimensionedField<Type, fvMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, fvMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new advectiveFvPatchField<Type>(*this, iF)
            );
        }


    // Member Functions

        using fvPatchField<Type>::db;
        using fvPatchField<Type>::time;
        using fvPatchField<Type>::patch;
        using fvPatchField<Type>::patchInternalField;
        using fvPatchField<Type>::internalField;

        // Evaluation functions

            //- Calculate and return the advection speed at the boundary
            virtual tmp<scalarField> advectionSpeed() const;

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

        //- Write data not value
        void writeData(Ostream&) const;

        //- Write
        virtual void write(Ostream&) const;


    // Member Operators

        //- Inherit assignment
        using mixedFvPatchField<Type>::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "advectiveFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

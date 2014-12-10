/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::porousBafflePressureFvPatchField

Group
    grpCoupledBoundaryConditions

Description
    This boundary condition provides a jump condition, using the \c cyclic
    condition as a base.

    The porous baffle introduces a pressure jump defined by:

        \f[
            \Delta p = -(D \mu U + 0.5 I \rho |U|^2 )L
        \f]

    where

    \vartable
        p      | pressure [Pa]
        \rho   | density [kg/m3]
        \mu    | laminar viscosity [Pa s]
        I      | inertial coefficient
        D      | Darcy coefficient
        L      | porous media length in the flow direction
    \endvartable


    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        patchType    | underlying patch type should be \c cyclic| yes |
        D            | Darcy coefficient       | yes         |
        I            | inertial coefficient    | yes         |
        L            | porous media length in the flow direction | yes |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            porousBafflePressure;
        patchType       cyclic;
        jump            uniform 0;
        I               1000000;
        D               0.001;
        L               0.1;
        value           uniform 0;
    }
    \endverbatim

Note
     The underlying \c patchType should be set to \c cyclic

SourceFiles
    porousBafflePressureFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef porousBafflePressureFvPatchField_H
#define porousBafflePressureFvPatchField_H

#include "fixedJumpFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class porousBafflePressureFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class porousBafflePressureFvPatchField
:
    public fixedJumpFvPatchField<Type>
{
    // Private data

        //- Darcy pressure loss coefficient
        scalar D_;

        //- Inertia pressure lost coefficient
        scalar I_;

        //- Porous media length
        scalar length_;


public:

    //- Runtime type information
    TypeName("porousBafflePressure");


    // Constructors

        //- Construct from patch and internal field
        porousBafflePressureFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        porousBafflePressureFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given porousBafflePressureFvPatchField
        // onto a new patch
        porousBafflePressureFvPatchField
        (
            const porousBafflePressureFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        porousBafflePressureFvPatchField
        (
            const porousBafflePressureFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new porousBafflePressureFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        porousBafflePressureFvPatchField
        (
            const porousBafflePressureFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new porousBafflePressureFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


//- Specialisation of the jump-condition for the pressure
template<>
void porousBafflePressureFvPatchField<scalar>::updateCoeffs();


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "porousBafflePressureFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

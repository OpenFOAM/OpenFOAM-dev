/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::atmBoundaryLayerInletEpsilonFvPatchScalarField

Group
    grpRASBoundaryConditions grpInletBoundaryConditions

Description
    This boundary condition specifies an inlet value for the turbulence
    dissipation, \f$\epsilon\f$ (\c epsilon), appropriate for atmospheric
    boundary layers (ABL), and designed to be used in conjunction with the
    \c ABLInletVelocity inlet velocity boundary condition.

        \f[
            \epsilon = \frac{(U^*)^3}{K(z - z_g + z_0)}
        \f]

    where
    \vartable
        U^*     | frictional velocity
        K       | Karman's constant
        z       | vertical co-ordinate [m]
        z_0     | surface roughness length [m]
        z_g     | minimum vlaue in z direction [m]
    \endvartable

    and:

        \f[
            U^* = K \frac{U_{ref}}{ln\left(\frac{Z_{ref} + z_0}{z_0}\right)}
        \f]

    where:
    \vartable
        U_{ref} | reference velocity at \f$Z_{ref}\f$ [m/s]
        Z_{ref} | reference height [m]
    \endvartable

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        z            | vertical co-ordinate [m] | yes        |
        kappa        | Karman's constanat      | no          | 0.41
        Uref         | reference velocity [m/s] | yes        |
        Href         | reference height [m]    | yes         |
        z0           | surface roughness length [m] | yes    |
        zGround      | minimum z co-ordinate [m] | yes       |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            atmBoundaryLayerInletEpsilon;
        z               1.0;
        kappa           0.41;
        Uref            1.0;
        Href            0.0;
        z0              uniform 0.0;
        zGround         uniform 0.0;
    }
    \endverbatim

    Reference:
        D.M. Hargreaves and N.G. Wright,  "On the use of the k-epsilon model
        in commercial CFD software to model the neutral atmospheric boundary
        layer", Journal of Wind Engineering and Industrial Aerodynamics
        95(2007), pp 355-369.

SourceFiles
    atmBoundaryLayerInletEpsilonFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef atmBoundaryLayerInletEpsilonFvPatchScalarField_H
#define atmBoundaryLayerInletEpsilonFvPatchScalarField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class atmBoundaryLayerInletEpsilonFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class atmBoundaryLayerInletEpsilonFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Direction of the z-coordinate
        vector z_;

        //- Von Karman constant
        const scalar kappa_;

        //- Reference velocity
        const scalar Uref_;

        //- Reference height
        const scalar Href_;

        //- Surface roughness length
        scalarField z0_;

        //- Minimum co-ordinate value in z direction
        scalarField zGround_;

        //- Frictional velocity
        scalarField Ustar_;


public:

    //- Runtime type information
    TypeName("atmBoundaryLayerInletEpsilon");


    // Constructors

        //- Construct from patch and internal field
        atmBoundaryLayerInletEpsilonFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        atmBoundaryLayerInletEpsilonFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  atmBoundaryLayerInletEpsilonFvPatchScalarField onto a new patch
        atmBoundaryLayerInletEpsilonFvPatchScalarField
        (
            const atmBoundaryLayerInletEpsilonFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new atmBoundaryLayerInletEpsilonFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        atmBoundaryLayerInletEpsilonFvPatchScalarField
        (
            const atmBoundaryLayerInletEpsilonFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new atmBoundaryLayerInletEpsilonFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return max value
            const scalarField& Ustar() const
            {
                return Ustar_;
            }

            //- Return z direction
            const vector& z() const
            {
                return z_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

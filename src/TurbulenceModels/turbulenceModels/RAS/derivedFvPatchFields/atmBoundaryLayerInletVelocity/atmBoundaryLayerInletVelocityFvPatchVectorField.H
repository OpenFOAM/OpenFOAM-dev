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
    Foam::atmBoundaryLayerInletVelocityFvPatchVectorField

Group
    grpRASBoundaryConditions grpInletBoundaryConditions

Description
    This boundary condition specifies a velocity inlet profile appropriate
    for atmospheric boundary layers (ABL).  The profile is derived from the
    friction velocity, flow direction and the direction of the parabolic
    co-ordinate \c z.

        \f[
            U = \frac{U^*}{K} ln\left(\frac{z - z_g + z_0}{z_0}\right)
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

    Reference:
        D.M. Hargreaves and N.G. Wright,  "On the use of the k-epsilon model
        in commercial CFD software to model the neutral atmospheric boundary
        layer", Journal of Wind Engineering and Industrial Aerodynamics
        95(2007), pp 355-369.

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        n            | flow direction          | yes         |
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
        type            atmBoundaryLayerInletVelocity;
        n               (0 1 0);
        z               1.0;
        kappa           0.41;
        Uref            1.0;
        Href            0.0;
        z0              uniform 0.0;
        zGround         uniform 0.0;
    }
    \endverbatim

Note
    D.M. Hargreaves and N.G. Wright recommend Gamma epsilon in the
    k-epsilon model should be changed from 1.3 to 1.11 for consistency.
    The roughness height (Er) is given by Er = 20 z0 following the same
    reference.

SourceFiles
    atmBoundaryLayerInletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef atmBoundaryLayerInletVelocityFvPatchVectorField_H
#define atmBoundaryLayerInletVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
       Class atmBoundaryLayerInletVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class atmBoundaryLayerInletVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Frictional velocity
        scalarField Ustar_;

        //- Flow direction
        vector n_;

        //- Direction of the z-coordinate
        vector z_;

        //- Surface roughness lenght
        scalarField z0_;

        //- Von Karman constant
        const scalar kappa_;

        //- Reference velocity
        const scalar Uref_;

        //- Reference hight
        const scalar Href_;

        //- Minimum corrdinate value in z direction
        scalarField zGround_;


public:

    //- Runtime type information
    TypeName("atmBoundaryLayerInletVelocity");


    // Constructors

        //- Construct from patch and internal field
        atmBoundaryLayerInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        atmBoundaryLayerInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        // atmBoundaryLayerInletVelocityFvPatchVectorField onto a new patch
        atmBoundaryLayerInletVelocityFvPatchVectorField
        (
            const atmBoundaryLayerInletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new atmBoundaryLayerInletVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        atmBoundaryLayerInletVelocityFvPatchVectorField
        (
            const atmBoundaryLayerInletVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new atmBoundaryLayerInletVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return Ustar
            const scalarField& Ustar() const
            {
                return Ustar_;
            }

            //- Return flow direction
            const vector& n() const
            {
                return n_;
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
                const fvPatchVectorField&,
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

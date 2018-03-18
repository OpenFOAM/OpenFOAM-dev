/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2016 OpenFOAM Foundation
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
    Foam::atmBoundaryLayer

Group
    grpRASBoundaryConditions grpInletBoundaryConditions

Description
    This class provides functions to evaluate the velocity and turbulence
    distributions appropriate for atmospheric boundary layers (ABL).

    The profile is derived from the friction velocity, flow direction and
    "vertical" direction:

        \f[
            U = \frac{U^*}{\kappa} ln\left(\frac{z - z_g + z_0}{z_0}\right)
        \f]

        \f[
            k = \frac{(U^*)^2}{\sqrt{C_mu}}
        \f]

        \f[
            \epsilon = \frac{(U^*)^3}{\kappa(z - z_g + z_0)}
        \f]

    where
    \vartable
        U^*     | Friction velocity
        \kappa  | von Karman's constant
        C_mu    | Turbulence viscosity coefficient
        z       | Vertical coordinate
        z_0     | Surface roughness height [m]
        z_g     | Minimum z-coordinate [m]
    \endvartable
    and
        \f[
            U^* = \kappa\frac{U_{ref}}{ln\left(\frac{Z_{ref} + z_0}{z_0}\right)}
        \f]
    where
    \vartable
        U_{ref} | Reference velocity at \f$Z_{ref}\f$ [m/s]
        Z_{ref} | Reference height [m]
    \endvartable

    Use in the atmBoundaryLayerInletVelocity, atmBoundaryLayerInletK and
    atmBoundaryLayerInletEpsilon boundary conditions.

    Reference:
        D.M. Hargreaves and N.G. Wright,  "On the use of the k-epsilon model
        in commercial CFD software to model the neutral atmospheric boundary
        layer", Journal of Wind Engineering and Industrial Aerodynamics
        95(2007), pp 355-369.

Usage
    \table
        Property     | Description                      | Required  | Default
        flowDir      | Flow direction                   | yes       |
        zDir         | Vertical direction               | yes       |
        kappa        | von Karman's constant            | no        | 0.41
        Cmu          | Turbulence viscosity coefficient | no        | 0.09
        Uref         | Reference velocity [m/s]         | yes       |
        Zref         | Reference height [m]             | yes       |
        z0           | Surface roughness height [m]     | yes       |
        zGround      | Minimum z-coordinate [m]         | yes       |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    ground
    {
        type            atmBoundaryLayerInletVelocity;
        flowDir         (1 0 0);
        zDir            (0 0 1);
        Uref            10.0;
        Zref            20.0;
        z0              uniform 0.1;
        zGround         uniform 0.0;
    }
    \endverbatim

Note
    D.M. Hargreaves and N.G. Wright recommend Gamma epsilon in the
    k-epsilon model should be changed from 1.3 to 1.11 for consistency.
    The roughness height (Er) is given by Er = 20 z0 following the same
    reference.

SourceFiles
    atmBoundaryLayer.C

\*---------------------------------------------------------------------------*/

#ifndef atmBoundaryLayer_H
#define atmBoundaryLayer_H

#include "fvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
       Class atmBoundaryLayer Declaration
\*---------------------------------------------------------------------------*/

class atmBoundaryLayer
{
    // Private data

        //- Flow direction
        vector flowDir_;

        //- Direction of the z-coordinate
        vector zDir_;

        //- Von Karman constant
        const scalar kappa_;

        //- Turbulent viscosity coefficient
        const scalar Cmu_;

        //- Reference velocity
        const scalar Uref_;

        //- Reference height
        const scalar Zref_;

        //- Surface roughness height
        scalarField z0_;

        //- Minimum coordinate value in z direction
        scalarField zGround_;

        //- Friction velocity
        scalarField Ustar_;


public:

    // Constructors

        //- Construct null
        atmBoundaryLayer();

        //- Construct from the coordinates field and dictionary
        atmBoundaryLayer(const vectorField& p, const dictionary&);

        //- Construct by mapping given
        // atmBoundaryLayer onto a new patch
        atmBoundaryLayer
        (
            const atmBoundaryLayer&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        atmBoundaryLayer(const atmBoundaryLayer&);


    // Member functions

        // Access

            //- Return flow direction
            const vector& flowDir() const
            {
                return flowDir_;
            }

            //- Return z-direction
            const vector& zDir() const
            {
                return zDir_;
            }

            //- Return friction velocity
            const scalarField& Ustar() const
            {
                return Ustar_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            void rmap(const atmBoundaryLayer&, const labelList&);


        // Evaluate functions

            //- Return the velocity distribution for the ATM
            tmp<vectorField> U(const vectorField& p) const;

            //- Return the turbulent kinetic energy distribution for the ATM
            tmp<scalarField> k(const vectorField& p) const;

            //- Return the turbulent dissipation rate distribution for the ATM
            tmp<scalarField> epsilon(const vectorField& p) const;


        //- Write
        void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

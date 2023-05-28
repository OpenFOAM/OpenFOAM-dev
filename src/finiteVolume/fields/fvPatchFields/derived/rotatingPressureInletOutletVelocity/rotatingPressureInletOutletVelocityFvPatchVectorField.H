/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::rotatingPressureInletOutletVelocityFvPatchVectorField

Description
    This velocity inlet/outlet boundary condition is applied to patches in a
    rotating frame where the pressure is specified.  A zero-gradient is applied
    for outflow (as defined by the flux); for inflow, the velocity is obtained
    from the flux with a direction normal to the patch faces.

Usage
    \table
        Property     | Description             | Required    | Default value
        phi          | flux field name         | no          | phi
        tangentialVelocity | tangential velocity field | no  |
        origin       | origin of rotation in Cartesian co-ordinates | yes|
        axis         | axis of rotation        | yes         |
        omega        | angular velocity of the frame [rad/s] | no |
        rpm          | angular velocity of the frame [rpm]   | no |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            rotatingPressureInletOutletVelocity;
        phi             phi;
        tangentialVelocity uniform (0 0 0);
        origin          (0 0 0);
        axis            (0 0 1);
        rpm             100;
    }
    \endverbatim

    The \c omega or \c rpm entries are Function1 of time, see Foam::Function1s.

    Note:
        Sign conventions:
      - positive flux (out of domain): apply zero-gradient condition
      - negative flux (into of domain): derive from the flux in the patch-normal
        direction

See also
    Foam::pressureInletOutletVelocityFvPatchVectorField
    Foam::Function1s

SourceFiles
    rotatingPressureInletOutletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef rotatingPressureInletOutletVelocityFvPatchVectorField_H
#define rotatingPressureInletOutletVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "pressureInletOutletVelocityFvPatchVectorField.H"
#include "omega.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
    Class rotatingPressureInletOutletVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class rotatingPressureInletOutletVelocityFvPatchVectorField
:
    public pressureInletOutletVelocityFvPatchVectorField
{
    // Private Data

        //- Origin of the rotation
        vector origin_;

        //- Axis of the rotation
        vector axis_;

        //- Angular velocity of the frame (rad/sec)
        Function1s::omega omega_;


    // Private Member Functions

        //- Calculate the tangentialVelocity from omega
        void calcTangentialVelocity();


public:

    //- Runtime type information
    TypeName("rotatingPressureInletOutletVelocity");


    // Constructors

        //- Construct from patch, internal field and dictionary
        rotatingPressureInletOutletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  rotatingPressureInletOutletVelocityFvPatchVectorField
        //  onto a new patch
        rotatingPressureInletOutletVelocityFvPatchVectorField
        (
            const rotatingPressureInletOutletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        rotatingPressureInletOutletVelocityFvPatchVectorField
        (
            const rotatingPressureInletOutletVelocityFvPatchVectorField&
        ) = delete;

        //- Copy constructor setting internal field reference
        rotatingPressureInletOutletVelocityFvPatchVectorField
        (
            const rotatingPressureInletOutletVelocityFvPatchVectorField&,
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
                new rotatingPressureInletOutletVelocityFvPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

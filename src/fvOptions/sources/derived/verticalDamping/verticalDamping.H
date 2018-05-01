/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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
    Foam::fv::verticalDamping

Description
    This fvOption applies an explicit damping force to components of the vector
    field in the direction of gravity. Its intended purpose is to damp the
    vertical motions of an interface in the region approaching an outlet so that
    no reflections are generated.

    Damping is achieved by applying a force to the momentum equation
    proportional to the momentum of the flow in the direction of gravity. The
    constant of proportionality is given by a coefficient \f$\lambda\f$ which
    has units of inverse-time. In the absence of any other forces this would
    generate an exponential decay of the vertical velocity.

        \f[
            \frac{d (m u_z)}{d t} = - \lambda m u_z
        \f]
        \f[
            u_z = u_{z0} e^{- \lambda t}
        \f]

    The coefficient \f$\lambda\f$ should be set based on the desired level of
    damping and the residence time of a perturbation through the damping zone.
    For example, if waves moving at 2 [m/s] are travelling through a damping
    zone 8 [m] in length, then the residence time is 4 [s]. If it is deemed
    necessary to damp for 5 time-scales, then \f$\lambda\f$ should be set to
    equal 5/(4 [s]) = 1.2 [1/s].

Usage
    Example usage:
    \verbatim
    verticalDamping1
    {
        type            verticalDamping;

        selectionMode   cellZone;
        cellZone        nearOutlet;

        lambda          [0 0 -1 0 0 0 0] 1; // Damping coefficient

        timeStart       0;
        duration        1e6;
    }
    \endverbatim

    Example usage with graduated onset:
    \verbatim
    verticalDamping1
    {
        type            verticalDamping;

        selectionMode   all;

        // Define the line along which to apply the graduation
        origin          (1200 0 0);
        direction       (1 0 0);

        // Or, define multiple lines
        // origins         ((1200 0 0) (1200 -300 0) (1200 300 0));
        // directions      ((1 0 0) (0 -1 0) (0 1 0));

        ramp
        {
            type        halfCosineRamp;
            start       0;
            duration    600;
        }

        lambda          [0 0 -1 0 0 0 0] 1; // Damping coefficient

        timeStart       0;
        duration        1e6;
    }
    \endverbatim

SourceFiles
    verticalDamping.C

\*---------------------------------------------------------------------------*/

#ifndef verticalDamping_H
#define verticalDamping_H

#include "cellSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type> class Function1;

namespace fv
{

/*---------------------------------------------------------------------------*\
                   Class verticalDamping Declaration
\*---------------------------------------------------------------------------*/

class verticalDamping
:
    public cellSetOption
{
private:

    // Private data

        //- Damping coefficient [1/s]
        dimensionedScalar lambda_;

        //- The ramping function
        autoPtr<Function1<scalar>> ramp_;

        //- Origins of the ramping coordinate
        vectorField origins_;

        //- Directions of increasing ramping coordinate
        vectorField directions_;


    // Private Member Functions

        //- Source term to momentum equation
        void add
        (
            const volVectorField& alphaRhoU,
            fvMatrix<vector>& eqn,
            const label fieldi
        );


public:

    //- Runtime type information
    TypeName("verticalDamping");


    // Constructors

        //- Construct from components
        verticalDamping
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~verticalDamping()
    {}


    // Member Functions

         // Add explicit and implicit contributions

            //- Source term to momentum equation
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const label fieldi
            );

            //- Source term to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldi
            );

            //- Source term to phase momentum equation
            virtual void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldi
            );


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

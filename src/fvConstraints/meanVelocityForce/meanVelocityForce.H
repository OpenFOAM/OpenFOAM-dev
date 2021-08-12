/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    Foam::fv::meanVelocityForce

Description
    Calculates and applies the force necessary to maintain the specified mean
    velocity.

Usage
    Example usage:
    \verbatim
    meanVelocityConstraint
    {
        type            meanVelocityForce;

        selectionMode   all;

        U               U;          // Name of velocity field

        Ubar            (10.0 0 0); // Desired mean velocity
        relaxation      0.2;        // Optional relaxation factor
    }
    \endverbatim

Note
    Currently only handles kinematic pressure (incompressible solvers).

SourceFiles
    meanVelocityForce.C

\*---------------------------------------------------------------------------*/

#ifndef meanVelocityForce_H
#define meanVelocityForce_H

#include "fvConstraint.H"
#include "fvCellSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class meanVelocityForce Declaration
\*---------------------------------------------------------------------------*/

class meanVelocityForce
:
    public fvConstraint
{
    // Private data

        //- The set of cells the fvConstraint applies to
        fvCellSet set_;

        //- Name of the velocity field
        word UName_;

        //- Average velocity
        vector Ubar_;

        //- Relaxation factor
        scalar relaxation_;

        //- Pressure gradient before correction
        mutable scalar gradP0_;

        //- Change in pressure gradient
        mutable scalar dGradP_;

        //- Matrix 1/A coefficients field pointer
        mutable autoPtr<volScalarField> rAPtr_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs();

        //- Calculate and return the magnitude of the mean velocity
        //  averaged over the selected cellSet
        virtual scalar magUbarAve(const volVectorField& U) const;

        //- Write the pressure gradient to file (for restarts etc)
        void writeProps(const scalar gradP) const;

        //- Correct driving force for a constant mass flow rate
        void update(fvMatrix<vector>& eqn) const;


public:

    //- Runtime type information
    TypeName("meanVelocityForce");


    // Constructors

        //- Construct from explicit source name and mesh
        meanVelocityForce
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construction
        meanVelocityForce(const meanVelocityForce&) = delete;


    // Member Functions

        //- Return the average velocity
        const vector& Ubar() const
        {
            return Ubar_;
        }

        //- Return the list of fields constrained by the fvConstraint
        virtual wordList constrainedFields() const;

        //- Add the momentum source and set the 1/A coefficient
        virtual bool constrain
        (
            fvMatrix<vector>& eqn,
            const word& fieldName
        ) const;

        //- Correct the pressure gradient
        virtual bool constrain(volVectorField& U) const;

        //- Update for mesh changes
        virtual void updateMesh(const mapPolyMesh&);

        //- Update for mesh motion
        virtual bool movePoints();

        //- Read source dictionary
        virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const meanVelocityForce&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

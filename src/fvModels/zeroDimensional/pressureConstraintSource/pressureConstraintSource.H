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
    Foam::fv::zeroDimensional::pressureConstraintSource

Description
    Zero-dimensional pressure constraint source.

    This model facilitates specification of a constant or time-varying
    pressure. It uses the fluid’s compressibility to determine the amount of
    mass necessary to add or remove to achieve the desired pressure. Iteration
    may be necessary to converge the constraint in the case of non-linear
    equations of state.

    Properties are added or removed with their current value. The model
    therefore represents a uniform expansion or contraction in infinite space.

    This only works for a compressible fluid. In the case of an incompressible
    fluid, a standard fixed value fvConstraint should be used instead to bypass
    solution of the pressure equation alltogether.

Usage
    Example usage:
    \verbatim
    {
        type            pressureConstraintSource;
        p               1e5;
    }
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef pressureConstraintSource_H
#define pressureConstraintSource_H

#include "constraintSource.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
namespace zeroDimensional
{

/*---------------------------------------------------------------------------*\
                Class pressureConstraintSource Declaration
\*---------------------------------------------------------------------------*/

class pressureConstraintSource
:
    public constraintSource
{
    // Private data

        //- The pressure value
        autoPtr<Function1<scalar>> p_;


    // Private member functions

        //- Non-virtual read
        void readCoeffs();

        //- Return the mass source
        tmp<volScalarField::Internal> dmdt() const;


public:

    //- Runtime type information
    TypeName("pressureConstraintSource");


    // Constructors

        //- Construct from dictionary
        pressureConstraintSource
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~pressureConstraintSource();


    // Member Functions

        //- Read dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace zeroDimensional
} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

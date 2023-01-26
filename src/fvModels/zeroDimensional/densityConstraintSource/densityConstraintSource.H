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
    Foam::fv::zeroDimensional::densityConstraintSource

Description
    Zero-dimensional density constraint source.

    This model facilitates specification of a time-varying density (or,
    equivalently, specific volume) be using a Function1. It adds or removes
    mass appropriately to achieve the desired density.

    Fluid properties are added or removed with their current value. This model
    therefore represents a uniform expansion or contraction in an infinite
    space.

    Note that in the case of a constant density this model does not do
    anything. A zero-dimension case naturally has both constant mass and
    volume.

    This only works for a compressible fluid.

Usage
    Example usage:
    \verbatim
    {
        type            densityConstraintSource;
        rho
        {
            type            scale;
            values
            (
                (0 1.16)
                (1 1.16)
                (1.1 2.02)
                (10 2.02)
            );
        }
    }
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef densityConstraintSource_H
#define densityConstraintSource_H

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
                Class densityConstraintSource Declaration
\*---------------------------------------------------------------------------*/

class densityConstraintSource
:
    public constraintSource
{
    // Private data

        //- The density value
        autoPtr<Function1<scalar>> rho_;


    // Private member functions

        //- Non-virtual read
        void readCoeffs();

        //- Return the mass source
        virtual tmp<volScalarField::Internal> dmdt() const;


public:

    //- Runtime type information
    TypeName("densityConstraintSource");


    // Constructors

        //- Construct from dictionary
        densityConstraintSource
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~densityConstraintSource();


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

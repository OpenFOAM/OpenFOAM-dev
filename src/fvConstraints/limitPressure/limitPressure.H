/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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
    Foam::fv::limitPressure

Description
    Limits the specified pressure field to be between specified minimum and
    maximum limits.

Usage
    Example usage:
    \verbatim
    limitp
    {
        type            limitPressure;

        // p               p_rgh; // Optional entry to specify the pressure

        min             200;
        max             500;
    }
    \endverbatim

SourceFiles
    limitPressure.C

\*---------------------------------------------------------------------------*/

#ifndef limitPressure_H
#define limitPressure_H

#include "fvConstraint.H"
#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class limitPressure Declaration
\*---------------------------------------------------------------------------*/

class limitPressure
:
    public fvConstraint
{
    // Private data

        //- Pressure field name, default = p
        word pName_;

        //- Minimum pressure limit
        dimensionedScalar pMin_;

        //- Maximum pressure limit
        dimensionedScalar pMax_;

        //- Pressure min-limited
        bool limitMinP_;

        //- Pressure max-limited
        bool limitMaxP_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs();


public:

    //- Runtime type information
    TypeName("limitPressure");


    // Constructors

        //- Construct from components
        limitPressure
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construction
        limitPressure(const limitPressure&) = delete;


    //- Destructor
    virtual ~limitPressure()
    {}


    // Member Functions

        //- Return the list of fields constrained by the fvConstraint
        virtual wordList constrainedFields() const;

        //- Constrain the energy field
        virtual bool constrain(volScalarField& he) const;

        //- Read dictionary
        virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const limitPressure&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

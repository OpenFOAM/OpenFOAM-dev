/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2021 OpenFOAM Foundation
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
    Foam::fv::limitTemperature

Description
    Correction for temperature to apply limits between minimum and maximum
    values.

Usage
    Example usage:
    \verbatim
    limitT
    {
        type            limitTemperature;

        selectionMode   all;

        phase           gas; // Optional phase name

        min             200;
        max             500;
    }
    \endverbatim

SourceFiles
    limitTemperature.C

\*---------------------------------------------------------------------------*/

#ifndef limitTemperature_H
#define limitTemperature_H

#include "cellSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class limitTemperature Declaration
\*---------------------------------------------------------------------------*/

class limitTemperature
:
    public cellSetOption
{
    // Private data

        //- Minimum temperature limit [K]
        scalar Tmin_;

        //- Maximum temperature limit [K]
        scalar Tmax_;

        //- Optional phase name
        word phaseName_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs();


public:

    //- Runtime type information
    TypeName("limitTemperature");


    // Constructors

        //- Construct from components
        limitTemperature
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construction
        limitTemperature(const limitTemperature&) = delete;


    //- Destructor
    virtual ~limitTemperature()
    {}


    // Member Functions

        //- Return the list of fields corrected by the option
        virtual wordList correctedFields() const;

        //- Correct the energy field
        virtual void correct(volScalarField& he) const;

        //- Read dictionary
        virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const limitTemperature&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

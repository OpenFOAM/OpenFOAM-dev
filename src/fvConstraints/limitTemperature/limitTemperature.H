/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2024 OpenFOAM Foundation
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
    Limits the temperature to be between minimum and maximum values.

Usage
    Example usage:
    \verbatim
    limitT
    {
        type            limitTemperature;

        select          all;

        // field           T;   // Optional energy/temperature field name
                                // Set to T for compressibleVoF

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

#include "fvConstraint.H"
#include "fvCellSet.H"

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
    public fvConstraint
{
    // Private data

        //- The set of cells the fvConstraint applies to
        fvCellSet set_;

        //- Minimum temperature limit [K]
        scalar Tmin_;

        //- Maximum temperature limit [K]
        scalar Tmax_;

        //- Optional energy/temperature field name
        word fieldName_;

        //- Optional phase name
        word phaseName_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs(const dictionary& dict);


public:

    //- Runtime type information
    TypeName("limitTemperature");


    // Constructors

        //- Construct from components
        limitTemperature
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        limitTemperature(const limitTemperature&) = delete;


    //- Destructor
    virtual ~limitTemperature()
    {}


    // Member Functions

        //- Return the list of fields constrained by the fvConstraint
        virtual wordList constrainedFields() const;

        //- Constrain the energy field
        virtual bool constrain(volScalarField& he) const;

        //- Update for mesh motion
        virtual bool movePoints();

        //- Update topology using the given map
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);

        //- Redistribute or update using the given distribution map
        virtual void distribute(const polyDistributionMap&);

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

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::fv::bXiTimedIgnition

Description
    Abstract base-class for timed ignition models
    for the Weller b-Xi combustion models

SourceFiles
    bXiTimedIgnition.C

\*---------------------------------------------------------------------------*/

#ifndef bXiTimedIgnition_H
#define bXiTimedIgnition_H

#include "bXiIgnition.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class bXiTimedIgnition Declaration
\*---------------------------------------------------------------------------*/

class bXiTimedIgnition
:
    public bXiIgnition
{
protected:

    // Protected Data

        //- Ignition start time
        dimensionedScalar start_;

        //- Ignition duration
        dimensionedScalar duration_;

        //- Ignition repeat period
        dimensionedScalar period_;

        //- Combustion duration
        //  after which ignited is set false and the composition reset
        dimensionedScalar combustionDuration_;

        //- Set true when the composition is reset
        //  at the end of the combustion duration
        //  and reset false at the beginning of ignition
        mutable bool reset_;


private:

    // Private Member Functions

        //- Non-virtual read
        void readCoeffs(const dictionary& dict);

        //- Return the time relative to the ignition time
        //  for the current period
        scalar ignRelTime(const scalar t) const;


public:

    //- Runtime type information
    TypeName("bXiTimedIgnition");


    // Constructors

        //- Construct from explicit source name and mesh
        bXiTimedIgnition
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        bXiTimedIgnition
        (
            const bXiTimedIgnition&
        ) = delete;


    // Member Functions

        // Checks

            //- Return the list of fields for which the option adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;

            //- Return true during the specified ignition duration
            virtual bool igniting(const dimensionedScalar duration) const;

            //- Return true during the ignition duration
            virtual bool igniting() const;

            //- Return true during the combustion duration
            virtual bool ignited() const;


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const bXiTimedIgnition&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

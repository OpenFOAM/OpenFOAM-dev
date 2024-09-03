/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::fv::multiCycleConstantbXiIgnition

Description
    Simple constant rate ignition model for the Weller b-Xi combustion models

Usage
    Example usage:
    \verbatim
    multiCycleConstantbXiIgnition1
    {
        type            multiCycleConstantbXiIgnition;

        select          cellZone;
        cellZone        ignition;

        start           0;
        duration        0.003;
        strength        2;

        XiCorr
        {
            type            spherical;
            select          all;
        }
    }
    \endverbatim

    Where:
    \table
        Property     | Description                 | Required   | Default value
        select       | Ignition cell set selection | yes        |
        start        | Ignition start time         | yes        |
        duration     | Ignition duration           | yes        |
        strength     | Ignition strength [1/s]     | yes        |
        XiCorr       | Flame-wrinkling correction  | yes        |
    \endtable

See also
    Foam::XiCorrModel
    Foam::XiCorrModels::planar
    Foam::XiCorrModels::cylindrical
    Foam::XiCorrModels::spherical

SourceFiles
    multiCycleConstantbXiIgnition.C

\*---------------------------------------------------------------------------*/

#ifndef multiCycleConstantbXiIgnition_H
#define multiCycleConstantbXiIgnition_H

#include "constantbXiIgnition.H"
#include "XiCorrModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                     Class multiCycleConstantbXiIgnition Declaration
\*---------------------------------------------------------------------------*/

class multiCycleConstantbXiIgnition
:
    public constantbXiIgnition
{
    // Private Data

        //- Ignition repeat period
        dimensionedScalar period_;

        //- Combustion duration
        //  after which ignited is set false and the composition reset
        dimensionedScalar combustionDuration_;

        //- Set true when the composition is reset
        //  at the end of the combustion duration
        //  and reset false at the beginning of ignition
        mutable bool reset_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs(const dictionary& dict);

        //- Return the time relative to the ignition time
        //  for the current period
        scalar ignRelTime(const scalar t) const;


public:

    //- Runtime type information
    TypeName("multiCycleConstantbXiIgnition");


    // Constructors

        //- Construct from explicit source name and mesh
        multiCycleConstantbXiIgnition
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        multiCycleConstantbXiIgnition
        (
            const multiCycleConstantbXiIgnition&
        ) = delete;


    // Member Functions

        // Checks

            //- Return true during the ignition duration
            virtual bool igniting() const;

            //- Return true during the combustion duration
            virtual bool ignited() const;


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const multiCycleConstantbXiIgnition&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

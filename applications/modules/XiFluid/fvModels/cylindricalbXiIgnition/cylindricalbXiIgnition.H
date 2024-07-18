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
    Foam::fv::cylindricalbXiIgnition

Description
    Simple constant rate ignition model for the Weller b-Xi combustion models

Usage
    Example usage:
    \verbatim
    cylindricalbXiIgnition1
    {
        type            cylindricalbXiIgnition;

        select          cellZone;
        cellZone        ignition;

        start           0;
        duration        0.003;
        strength        2;

        XiCorr
        {
            select          all;
            thickness       0.001;
            cylinderFraction 0.25;
        }
    }
    \endverbatim

    Where:
    \table
        Property         | Description              | Required  | Default value
        start            | Ignition start time      | yes       |
        duration         | Ignition duration        | yes       |
        strength         | Ignition strength [1/s]  | yes       |
        thickness        | Thickness of the 2D mesh | yes       |
        cylinderFraction | Kernel cylinder fraction | no        | 1
    \endtable

SourceFiles
    cylindricalbXiIgnition.C

\*---------------------------------------------------------------------------*/

#ifndef cylindricalbXiIgnition_H
#define cylindricalbXiIgnition_H

#include "constantbXiIgnition.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                     Class cylindricalbXiIgnition Declaration
\*---------------------------------------------------------------------------*/

class cylindricalbXiIgnition
:
    public constantbXiIgnition
{
    // Private Data

        //- Thickness of the 2D mesh
        dimensionedScalar thickness_;

        //- Kernel cylinder fraction
        //  If the ignition is adjacent to n symmetry planes: 1/2^n
        dimensionedScalar cylinderFraction_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs();

        //- Return the area of the ignition kernel
        //  calculated from the given kernel volume
        virtual dimensionedScalar Ak(const dimensionedScalar& Vk) const;


public:

    //- Runtime type information
    TypeName("cylindricalbXiIgnition");


    // Constructors

        //- Construct from explicit source name and mesh
        cylindricalbXiIgnition
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        cylindricalbXiIgnition
        (
            const cylindricalbXiIgnition&
        ) = delete;


    // Member Functions

        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const cylindricalbXiIgnition&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

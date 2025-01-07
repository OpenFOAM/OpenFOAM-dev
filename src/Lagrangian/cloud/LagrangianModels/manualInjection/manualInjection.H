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
    Foam::manualInjection

Description
    Manual injection model. This injects particles instantaneously at a list of
    positions specified in a separate file.

Usage
    \table
        Property    | Description                     | Required? | Default
        file        | Name of the file containing the positions | yes |
        units       | Units of the positions in the file | no     | [m]
        time        | The time at which to inject     | no        | 0
    \endtable

    Example specification:
    \verbatim
    <LagrangianModelName>
    {
        type        manualInjection;
        file        "manualInjectionPositions";
        units       [mm];
    }
    \endverbatim

SourceFiles
    manualInjection.C

\*---------------------------------------------------------------------------*/

#ifndef manualInjection_H
#define manualInjection_H

#include "LagrangianInjection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Lagrangian
{

/*---------------------------------------------------------------------------*\
                          Class manualInjection Declaration
\*---------------------------------------------------------------------------*/

class manualInjection
:
    public LagrangianInjection
{
private:

    // Private Data

        //- File containing positions
        fileName fileName_;

        //- Optional unit conversion for the positions
        unitConversion units_;

        //- The time at which to inject
        scalar time_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs(const dictionary& modelDict);


public:

    //- Runtime type information
    TypeName("manualInjection");


    // Constructors

        //- Construct from components
        manualInjection
        (
            const word& name,
            const LagrangianMesh& mesh,
            const dictionary& modelDict,
            const dictionary& stateDict
        );


    // Member Functions

        //- Create new elements in the Lagrangian mesh
        LagrangianSubMesh modify
        (
            LagrangianMesh& mesh,
            const LagrangianSubMesh&
        ) const;


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& modelDict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Lagrangian
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

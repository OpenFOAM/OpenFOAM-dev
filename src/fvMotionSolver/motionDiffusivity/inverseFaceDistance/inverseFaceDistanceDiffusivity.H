/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::inverseFaceDistanceDiffusivity

Description
    Inverse distance to the given patches motion diffusivity.

SourceFiles
    inverseFaceDistanceDiffusivity.C

\*---------------------------------------------------------------------------*/

#ifndef inverseFaceDistanceDiffusivity_H
#define inverseFaceDistanceDiffusivity_H

#include "motionDiffusivity.H"
#include "wordReList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class inverseFaceDistanceDiffusivity Declaration
\*---------------------------------------------------------------------------*/

class inverseFaceDistanceDiffusivity
:
    public motionDiffusivity
{
    // Private Data

        //- Patches selected to base the distance on
        wordReList patchNames_;


public:

    //- Runtime type information
    TypeName("inverseFaceDistance");


    // Constructors

        //- Construct for the given fvMesh and data Istream
        inverseFaceDistanceDiffusivity(const fvMesh& mesh, Istream& mdData);

        //- Disallow default bitwise copy construction
        inverseFaceDistanceDiffusivity
        (
            const inverseFaceDistanceDiffusivity&
        ) = delete;


    //- Destructor
    virtual ~inverseFaceDistanceDiffusivity();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const inverseFaceDistanceDiffusivity&) = delete;

        //- Return diffusivity field
        virtual tmp<surfaceScalarField> operator()() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

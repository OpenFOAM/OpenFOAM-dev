/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::wallDist

Description
    Specialisation of patchDist for wall distance calculation

SourceFiles
    wallDist.C

\*---------------------------------------------------------------------------*/

#ifndef wallDist_H
#define wallDist_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class wallDist Declaration
\*---------------------------------------------------------------------------*/

class wallDist
:
    public volScalarField
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        wallDist(const wallDist&);

        //- Disallow default bitwise assignment
        void operator=(const wallDist&);


public:

    // Constructors

        //- Construct from mesh
        wallDist(const fvMesh& mesh);


    //- Destructor
    virtual ~wallDist();


    // Member Functions

        const volScalarField& y() const
        {
            return *this;
        }

        //- Correct for mesh changes
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

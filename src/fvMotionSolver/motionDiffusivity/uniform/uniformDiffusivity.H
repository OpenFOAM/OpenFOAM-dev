/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::uniformDiffusivity

Description
    Uniform uniform finite volume mesh motion diffusivity.

SourceFiles
    uniformDiffusivity.C

\*---------------------------------------------------------------------------*/

#ifndef uniformDiffusivity_H
#define uniformDiffusivity_H

#include "motionDiffusivity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class uniformDiffusivity Declaration
\*---------------------------------------------------------------------------*/

class uniformDiffusivity
:
    public motionDiffusivity
{
public:

    //- Runtime type information
    TypeName("uniform");


    // Constructors

        //- Construct for the given fvMesh and data Istream
        uniformDiffusivity(const fvMesh& mesh, Istream& mdData);

        //- Disallow default bitwise copy construction
        uniformDiffusivity(const uniformDiffusivity&) = delete;


    //- Destructor
    virtual ~uniformDiffusivity();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const uniformDiffusivity&) = delete;

        //- Return diffusivity field
        virtual tmp<surfaceScalarField> operator()() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

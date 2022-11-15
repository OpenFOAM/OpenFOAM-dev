/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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
    Foam::twoPhases

Description
    Interface to two phases

SourceFiles
    twoPhases.C

\*---------------------------------------------------------------------------*/

#ifndef twoPhases_H
#define twoPhases_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class twoPhases Declaration
\*---------------------------------------------------------------------------*/

class twoPhases
{
public:

    //- Runtime type information
    TypeName("twoPhases");


    // Constructors

        //- Default construct
        twoPhases()
        {}


    //- Destructor
    virtual ~twoPhases()
    {}


    // Member Functions

        //- Access the mesh
        const fvMesh& mesh() const
        {
            return alpha1().mesh();
        }

        //- Return the index associated with a given phase name
        inline bool index(const word& name) const
        {
            if (alpha1().group() == name)
            {
                return false;
            }

            if (alpha2().group() == name)
            {
                return true;
            }

            FatalErrorInFunction
                << "Index of phase " << name << " requested from phases"
                << alpha1().group() << " and " << alpha2().group()
                << exit(FatalError);
            return false;
        }

        //- Return the volume-fraction of phase 1
        virtual const volScalarField& alpha1() const = 0;

        //- Return the volume-fraction of phase 2
        virtual const volScalarField& alpha2() const = 0;

        //- Return the volume-fraction of a given phase
        const volScalarField& alpha(const bool index) const
        {
            return index ? alpha2() : alpha1();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

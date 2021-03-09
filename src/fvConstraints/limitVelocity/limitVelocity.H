/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2021 OpenFOAM Foundation
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
    Foam::fv::limitVelocity

Description
    Limits the velocity magnitude to the specified \c max value.

Usage
    Example usage:
    \verbatim
    limitU
    {
        type            limitVelocity;

        selectionMode   all;

        max             100;
    }
    \endverbatim

SourceFiles
    limitVelocity.C

\*---------------------------------------------------------------------------*/

#ifndef limitVelocity_H
#define limitVelocity_H

#include "fvConstraint.H"
#include "fvCellSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class limitVelocity Declaration
\*---------------------------------------------------------------------------*/

class limitVelocity
:
    public fvConstraint
{
    // Private data

        //- The set of cells the fvConstraint applies to
        fvCellSet set_;

        //- Velocity field name, default = U
        word UName_;

        //- Maximum velocity magnitude
        scalar max_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs();


public:

    //- Runtime type information
    TypeName("limitVelocity");


    // Constructors

        //- Construct from components
        limitVelocity
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construction
        limitVelocity(const limitVelocity&) = delete;


    //- Destructor
    virtual ~limitVelocity()
    {}


    // Member Functions

        //- Return the list of fields constrained by the fvConstraint
        virtual wordList constrainedFields() const;

        //- Constrain the velocity field
        virtual void constrain(volVectorField& U) const;

        //- Update for mesh changes
        virtual void updateMesh(const mapPolyMesh&);

        //- Read dictionary
        virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const limitVelocity&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

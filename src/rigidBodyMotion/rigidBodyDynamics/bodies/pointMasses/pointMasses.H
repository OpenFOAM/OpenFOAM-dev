/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2026 OpenFOAM Foundation
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
    Foam::RBD::pointMasses

Description
    Specialisation of rigidBody to construct a pointMasses given the mass and
    lengths of the sides.

SourceFiles
    pointMassesI.H
    pointMasses.C

\*---------------------------------------------------------------------------*/

#ifndef RBD_pointMasses_H
#define RBD_pointMasses_H

#include "rigidBody.H"
#include "TableReader.H"
#include "pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{

/*---------------------------------------------------------------------------*\
                           Class pointMasses Declaration
\*---------------------------------------------------------------------------*/

class pointMasses
:
    public rigidBody
{
    // Private member data

        //- Table reader
        const autoPtr<TableReader<point, scalar>> reader_;

        //- Points
        pointField points_;

        //- Masses
        scalarField masses_;


public:

    //- Runtime type information
    TypeName("pointMasses");


    // Constructors

        //- Construct from dictionary
        inline pointMasses
        (
            const word& name,
            const dictionary& dict
        );

        //- Return clone of this pointMasses
        virtual autoPtr<rigidBody> clone() const;


    //- Destructor
    virtual ~pointMasses();


    // Member Functions

        //- Return the lengths of the sides of the pointMasses
        inline const vector& L() const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RBD
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

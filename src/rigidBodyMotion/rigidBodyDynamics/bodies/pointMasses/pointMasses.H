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

        //- Indices of the points when sorted in the three axes
        mutable PtrList<labelList> sortedPointis_;


    // Private member functions

        //- Return the n-th moments of the sections of the body between
        //  planes defined by distances along the given axis
        template<class Type, class Op>
        tmp<Field<Type>> sectionMuNs
        (
            const direction axis,
            const scalarField& distances,
            const Op& op
        ) const;


public:

    //- Runtime type information
    TypeName("pointMasses");


    // Constructors

        //- Construct from dictionary
        pointMasses
        (
            const word& name,
            const dictionary& dict
        );

        //- Copy construct
        pointMasses(const pointMasses&);

        //- Return clone of this pointMasses
        virtual autoPtr<rigidBody> clone() const;


    //- Destructor
    virtual ~pointMasses();


    // Member Functions

        //- Return the zeroth moments of the sections of the body between
        //  planes defined by distances along the given axis
        virtual tmp<scalarField> sectionMu0s
        (
            const direction axis,
            const scalarField& distances
        ) const;

        //- Return the first moments of the sections of the body between
        //  planes defined by distances along the given axis
        virtual tmp<vectorField> sectionMu1s
        (
            const direction axis,
            const scalarField& distances
        ) const;

        //- Return the second moments of the sections of the body between
        //  planes defined by distances along the given axis
        virtual tmp<symmTensorField> sectionMu2s
        (
            const direction axis,
            const scalarField& distances
        ) const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RBD
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

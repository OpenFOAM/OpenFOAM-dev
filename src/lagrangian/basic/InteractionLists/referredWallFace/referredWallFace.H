/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::referredWallFace

Description
    Storage for referred wall faces.  Stores patch index, face and
    associated points

SourceFiles
    referredWallFaceI.H
    referredWallFace.C

\*---------------------------------------------------------------------------*/

#ifndef referredWallFace_H
#define referredWallFace_H

#include "face.H"
#include "pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class referredWallFace;
Istream& operator>>(Istream&, referredWallFace&);
Ostream& operator<<(Ostream&, const referredWallFace&);


/*---------------------------------------------------------------------------*\
                        Class referredWallFace Declaration
\*---------------------------------------------------------------------------*/

class referredWallFace
:
    public face
{
    // Private data

        //- Points of face
        pointField pts_;

        //- Index of originating patch
        label patchI_;


public:

    // Constructors

        //- Construct null
        referredWallFace();

        //- Construct from components
        referredWallFace
        (
            const face& f,
            const pointField& pts,
            label patchI
        );

        //- Construct as copy
        referredWallFace(const referredWallFace&);


    //- Destructor
    ~referredWallFace();


    // Member Functions

        // Access

            //- Return access to the stored points
            inline const pointField& points() const;

            //- Return non-const access to the stored points
            inline pointField& points();

            //- Return access to the patch index
            inline label patchIndex() const;

            //- Return non-const access to the patch index
            inline label& patchIndex();


    // Member Operators

        bool operator==(const referredWallFace&) const;
        bool operator!=(const referredWallFace&) const;

    // IOstream Operators

        friend Istream& operator>>(Istream&, referredWallFace&);
        friend Ostream& operator<<(Ostream&, const referredWallFace&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "referredWallFaceI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

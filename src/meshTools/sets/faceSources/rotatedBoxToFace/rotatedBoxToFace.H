/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2024 OpenFOAM Foundation
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
    Foam::rotatedBoxToFace

Description
    A topoSetSource to select faces based on cell centres inside a rotated
    and/or skewed box.

    The box can be defined with an origin and three vectors; i, j, and k. The
    origin is one corner of the box, and the vectors are the edges connected to
    that corner.

    For example, the following defines a box rotated 45 degrees around the
    z-axis, with width and depth of 0.2, height of 200, and with a bottom left
    corner at (0.4 0.4 -100):

    \verbatim
        origin  (0.4 0.4 -100);
        i       (0.141421 0.141421 0);
        j       (-0.141421 0.141421 0);
        k       (0 0 200);
    \endverbatim

    Alternatively, the box can be defined using a non-rotated box and details
    of how it should be rotated. This syntax is triggered by the presence of
    the keyword "box". A standard bounding box is supplied, along with a centre
    of rotation and two vectors, n1 and n2. The rotation is taken to be that
    which transforms n1 onto n2.

    The above example can be equivalently specified in this alternative form as
    follows:

    \verbatim
        box     (0.4 0.4 -100) (0.6 0.6 100);
        centre  (0.4 0.4 0);
        n1      (1 0 0);
        n2      (1 1 0);
    \endverbatim

SourceFiles
    rotatedBoxToFace.C

\*---------------------------------------------------------------------------*/

#ifndef rotatedBoxToFace_H
#define rotatedBoxToFace_H

#include "topoSetSource.H"
#include "treeBoundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class rotatedBoxToFace Declaration
\*---------------------------------------------------------------------------*/

class rotatedBoxToFace
:
    public topoSetSource
{
    // Private Data

        //- Skewed box origin
        vector origin_;

        //- Skewed box edge vector
        vector i_;

        //- Skewed box edge vector
        vector j_;

        //- Skewed box edge vector
        vector k_;


    // Private Member Functions

        void combine(topoSet& set, const bool add) const;


public:

    //- Runtime type information
    TypeName("rotatedBoxToFace");


    // Constructors

        //- Construct from components
        rotatedBoxToFace
        (
            const polyMesh& mesh,
            const vector& origin,
            const vector& i,
            const vector& j,
            const vector& k
        );

        //- Construct from dictionary
        rotatedBoxToFace(const polyMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~rotatedBoxToFace();


    // Member Functions

        virtual sourceType setType() const
        {
            return FACESETSOURCE;
        }

        virtual void applyToSet
        (
            const topoSetSource::setAction action,
            topoSet&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

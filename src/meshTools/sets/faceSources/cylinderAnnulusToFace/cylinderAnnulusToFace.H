/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2024 OpenFOAM Foundation
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
    Foam::cylinderAnnulusToFace

Description
    A topoSetSource to select faces based on face centres inside a
    cylinder annulus.

SourceFiles
    cylinderAnnulusToFace.C

\*---------------------------------------------------------------------------*/

#ifndef cylinderAnnulusToFace_H
#define cylinderAnnulusToFace_H

#include "topoSetSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class cylinderAnnulusToFace Declaration
\*---------------------------------------------------------------------------*/

class cylinderAnnulusToFace
:
    public topoSetSource
{
    // Private Data

        //- First point on cylinder axis
        vector point1_;

        //- Second point on cylinder axis
        vector point2_;

        //- Outer Radius
        scalar outerRadius_;

        //- Inner Radius
        scalar innerRadius_;


    // Private Member Functions

        void combine(topoSet& set, const bool add) const;


public:

    //- Runtime type information
    TypeName("cylinderAnnulusToFace");


    // Constructors

        //- Construct from components
        cylinderAnnulusToFace
        (
            const polyMesh& mesh,
            const vector& point1,
            const vector& point2,
            const scalar outerRadius,
            const scalar innerRadius
        );

        //- Construct from dictionary
        cylinderAnnulusToFace
        (
            const polyMesh& mesh,
            const dictionary& dict
        );


    // Destructor
    virtual ~cylinderAnnulusToFace();


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

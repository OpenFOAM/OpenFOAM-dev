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
    Foam::setToPointZone

Description
    A topoSetSource to select points based on usage in a pointSet.

SourceFiles
    setToPointZone.C

\*---------------------------------------------------------------------------*/

#ifndef setToPointZone_H
#define setToPointZone_H

#include "topoSetSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class setToPointZone Declaration
\*---------------------------------------------------------------------------*/

class setToPointZone
:
    public topoSetSource
{
    // Private Data

        //- Add usage string
        static addToUsageTable usage_;

        //- Name of set to use
        word setName_;

public:

    //- Runtime type information
    TypeName("setToPointZone");

    // Constructors

        //- Construct from components
        setToPointZone
        (
            const polyMesh& mesh,
            const word& setName
        );

        //- Construct from dictionary
        setToPointZone
        (
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Construct from Istream
        setToPointZone
        (
            const polyMesh& mesh,
            Istream&
        );


    //- Destructor
    virtual ~setToPointZone();


    // Member Functions

        virtual sourceType setType() const
        {
            return POINTZONESOURCE;
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

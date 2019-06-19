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
    Foam::setUpdater

Description
    Keeps cell/face/vertex sets uptodate. Both the ones loaded and the ones
    on disk.

SourceFiles
    setUpdater.C
    setUpdaterTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef setUpdater_H
#define setUpdater_H

#include "polyMeshModifier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class setUpdater Declaration
\*---------------------------------------------------------------------------*/

class setUpdater
:
    public polyMeshModifier
{
    // Private Member Functions

        //- Updates all sets
        template<class Type>
        void updateSets(const mapPolyMesh& morphMap) const;


public:

    //- Runtime type information
    TypeName("setUpdater");


    // Constructors

        //- Construct from dictionary
        setUpdater
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyTopoChanger& mme
        );

        //- Disallow default bitwise copy construction
        setUpdater(const setUpdater&) = delete;


    //- Destructor
    virtual ~setUpdater();


    // Member Functions

        //- Check for topology change
        virtual bool changeTopology() const;

        //- Insert the layer addition/removal instructions
        //  into the topological change
        virtual void setRefinement(polyTopoChange&) const;

        //- Modify motion points to comply with the topological change
        virtual void modifyMotionPoints(pointField& motionPoints) const;

        //- Force recalculation of locally stored data on topological change
        virtual void updateMesh(const mapPolyMesh&);

        //- Write
        virtual void write(Ostream&) const;

        //- Write dictionary
        virtual void writeDict(Ostream&) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const setUpdater&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
    #include "setUpdaterTemplates.C"
#endif


#endif

// ************************************************************************* //

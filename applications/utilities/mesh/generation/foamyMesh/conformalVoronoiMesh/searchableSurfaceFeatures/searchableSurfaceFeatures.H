/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2019 OpenFOAM Foundation
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
    Foam::searchableSurfaceFeatures

Description
    Decorator that returns the features of a searchable surface.

SourceFiles
    searchableSurfaceFeatures.C

\*---------------------------------------------------------------------------*/

#ifndef searchableSurfaceFeatures_H
#define searchableSurfaceFeatures_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "searchableSurface.H"
#include "extendedFeatureEdgeMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class searchableSurfaceFeatures Declaration
\*---------------------------------------------------------------------------*/

class searchableSurfaceFeatures
{
    // Private Data

        const searchableSurface& surface_;

        const dictionary& dict_;


public:

    //- Runtime type information
    TypeName("searchableSurfaceFeatures");

    // Declare run-time constructor selection table

        // For the dictionary constructor
        declareRunTimeSelectionTable
        (
            autoPtr,
            searchableSurfaceFeatures,
            dict,
            (
                const searchableSurface& surface,
                const dictionary& dict
            ),
            (surface, dict)
        );


    // Constructors

        //- Construct from components
        searchableSurfaceFeatures
        (
            const searchableSurface& surface,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        searchableSurfaceFeatures(const searchableSurfaceFeatures&) = delete;

        //- Clone
        virtual autoPtr<searchableSurfaceFeatures> clone() const
        {
            NotImplemented;
            return autoPtr<searchableSurfaceFeatures>(nullptr);
        }


    // Selectors

        //- Return a reference to the selected searchableSurfaceFeatures
        static autoPtr<searchableSurfaceFeatures> New
        (
            const searchableSurface& surface,
            const dictionary& dict
        );


    //- Destructor
    virtual ~searchableSurfaceFeatures();


    // Member Functions

        //- Return a reference to the searchable surface
        const searchableSurface& surface() const
        {
            return surface_;
        }

        //- Return a reference to the dictionary
        const dictionary& dict() const
        {
            return dict_;
        }

        //- Return whether this searchable surface has features
        virtual bool hasFeatures() const = 0;

        //- Return an extendedFeatureEdgeMesh containing the features
        virtual autoPtr<extendedFeatureEdgeMesh> features() const
        {
            return autoPtr<extendedFeatureEdgeMesh>();
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const searchableSurfaceFeatures&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
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
    Foam::searchableBoxFeatures

Description

SourceFiles
    searchableBoxFeatures.C

\*---------------------------------------------------------------------------*/

#ifndef searchableBoxFeatures_H
#define searchableBoxFeatures_H

#include "searchableSurfaceFeatures.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class searchableBoxFeatures Declaration
\*---------------------------------------------------------------------------*/

class searchableBoxFeatures
:
    public searchableSurfaceFeatures
{
private:

    // Private Member Data

        //- Which side of the box to mesh
        const extendedFeatureEdgeMesh::sideVolumeType mode_;


public:

    //- Runtime type information
    TypeName("searchableBoxFeatures");


    // Constructors

        //- Construct from searchable surface and dictionary
        searchableBoxFeatures
        (
            const searchableSurface& surface,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        searchableBoxFeatures(const searchableBoxFeatures&) = delete;


    //- Destructor
    virtual ~searchableBoxFeatures();


    // Member Functions

        //- Return true for a searchable box having features
        virtual bool hasFeatures() const
        {
            return true;
        }

        //- Return an extendedFeatureEdgeMesh containing the features
        virtual autoPtr<extendedFeatureEdgeMesh> features() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const searchableBoxFeatures&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

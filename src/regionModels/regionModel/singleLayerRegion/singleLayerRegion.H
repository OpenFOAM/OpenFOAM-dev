/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::regionModels::singleLayerRegion

Description
    Base class for single layer region models

SourceFiles
    singleLayerRegion.C

\*---------------------------------------------------------------------------*/

#ifndef singleLayerRegion_H
#define singleLayerRegion_H

#include "regionModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{

/*---------------------------------------------------------------------------*\
                      Class singleLayerRegion Declaration
\*---------------------------------------------------------------------------*/

class singleLayerRegion
:
    public regionModel
{
    // Private member data

        // Region addressing

            //- Patch normal vectors
            volVectorField nHat_;

            //- Cell cross-sectional area magnitude [m^2]
            volScalarField::Internal magSf_;

            //- Cell layer volume/area [m]
            volScalarField VbyA_;


        // Addressing

            //- List of patch IDs opposite to internally coupled patches
            labelList passivePatchIDs_;


public:

    //- Runtime type information
    TypeName("regionModel");


    // Constructors

        //- Construct from mesh, region type and name
        singleLayerRegion
        (
            const fvMesh& mesh,
            const word& regionType,
            const word& modelName,
            bool readFields = true
        );

        //- Disallow default bitwise copy construction
        singleLayerRegion(const singleLayerRegion&) = delete;


    //- Destructor
    virtual ~singleLayerRegion();


    // Member Functions

        // Access

            // Region geometry

                //- Return the patch normal vectors
                const volVectorField& nHat() const;

                //- Return the face area magnitudes [m^2]
                const volScalarField::Internal& magSf() const;

                //- Return the cell layer volume/area [m]
                const volScalarField& VbyA() const;


            // Addressing

                //- Return the list of patch IDs opposite to internally
                //  coupled patches
                const labelList& passivePatchIDs() const;


        // Patch type information

            //- Return boundary types for mapped field patches
            //  Also maps internal field value
            //  Mapping region prescribed by underlying mapped poly patch
            template<class Type>
            wordList mappedFieldAndInternalPatchTypes() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const singleLayerRegion&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "singleLayerRegionTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

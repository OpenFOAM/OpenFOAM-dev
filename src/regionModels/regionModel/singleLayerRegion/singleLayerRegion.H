/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::singleLayerRegion

Description
    Base class for single layer region models

SourceFiles
    singleLayerRegion.C

\*---------------------------------------------------------------------------*/

#ifndef singleLayerRegion_H
#define singleLayerRegion_H

#include "regionModel.H"

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

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        singleLayerRegion(const singleLayerRegion&);

        //- Disallow default bitwise assignment
        void operator=(const singleLayerRegion&);

        //- Construct region mesh and fields
        void constructMeshObjects();

        //- Initialise the region
        void initialise();


protected:

    // Protected data

        // Region addressing

            //- Patch normal vectors
            autoPtr<volVectorField> nHatPtr_;

            //- Face area magnitudes / [m2]
            autoPtr<volScalarField> magSfPtr_;


        // Addressing

            //- List of patch IDs opposite to internally coupled patches
            labelList passivePatchIDs_;


    // Protected member functions

        //- Read control parameters from dictionary
        virtual bool read();


public:

    //- Runtime type information
    TypeName("regionModel");


    // Constructors

        //- Construct null
        singleLayerRegion(const fvMesh& mesh, const word& regionType);

        //- Construct from mesh, region type and name
        singleLayerRegion
        (
            const fvMesh& mesh,
            const word& regionType,
            const word& modelName,
            bool readFields = true
        );


    //- Destructor
    virtual ~singleLayerRegion();


    // Member Functions

        // Access

            // Region geometry

                //- Return the patch normal vectors
                virtual const volVectorField& nHat() const;

                //- Return the face area magnitudes / [m2]
                virtual const volScalarField& magSf() const;


            // Addressing

                //- Return the list of patch IDs opposite to internally
                //  coupled patches
                virtual const labelList& passivePatchIDs() const;


        // Patch type information

            //- Return boundary types for mapped field patches
            //  Also maps internal field value
            //  Mapping region prescribed by underlying mapped poly patch
            template<class Type>
            wordList mappedFieldAndInternalPatchTypes() const;

            //- Return boundary types for pushed mapped field patches
            //  Mapping region prescribed by underlying mapped poly patch
            template<class Type>
            wordList mappedPushedFieldPatchTypes() const;
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

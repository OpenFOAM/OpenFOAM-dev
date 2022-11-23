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
    Foam::singleLayerRegionModel

Description
    Base class for single layer region models

SourceFiles
    singleLayerRegionModel.C

\*---------------------------------------------------------------------------*/

#ifndef singleLayerRegionModel_H
#define singleLayerRegionModel_H

#include "fvMesh.H"
#include "timeIOdictionary.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class singleLayerRegionModel Declaration
\*---------------------------------------------------------------------------*/

class singleLayerRegionModel
:
    public IOdictionary
{

protected:

    // Protected data

        //- Reference to the primary mesh database
        const fvMesh& primaryMesh_;

        //- Reference to the time database
        const Time& time_;

        //- Information output
        Switch infoOutput_;

        //- Model name
        const word modelName_;

        //- Region name
        word regionName_;

        //- Region mesh
        fvMesh mesh_;

        //- Model coefficients dictionary
        dictionary coeffs_;

        //- Dictionary of output properties
        autoPtr<IOdictionary> outputPropertiesPtr_;


        // Addressing

            //- List of patch IDs on the primary region coupled to this region
            labelList primaryPatchIDs_;

            //- List of patch IDs internally coupled with the primary region
            labelList intCoupledPatchIDs_;

            //- List of patch IDs opposite to internally coupled patches
            labelList passivePatchIDs_;


        //- Patch normal vectors
        volVectorField nHat_;

        //- Cell cross-sectional area magnitude [m^2]
        volScalarField::Internal magSf_;

        //- Cell layer volume/area [m]
        volScalarField VbyA_;


    // Protected member functions

        //- Read control parameters from dictionary
        virtual bool read();


public:

    //- Runtime type information
    TypeName("singleLayerRegionModel");


    // Constructors

        //- Construct from mesh, region type and name
        singleLayerRegionModel
        (
            const fvMesh& primaryMesh,
            const word& regionType,
            const word& modelName,
            bool readFields = true
        );

        //- Disallow default bitwise copy construction
        singleLayerRegionModel(const singleLayerRegionModel&) = delete;


    //- Destructor
    virtual ~singleLayerRegionModel();


    // Member Functions

        // Access

            //- Return the reference to the primary mesh database
            inline const fvMesh& primaryMesh() const;

            //- Return the reference to the time database
            inline const Time& time() const;

            //- Return the information flag
            inline const Switch& infoOutput() const;

            //- Return the model name
            inline const word& modelName() const;

            //- Return the region mesh
            inline const fvMesh& mesh() const;

            //- Return the region mesh database for manipulation
            inline fvMesh& mesh();

            //- Return the model coefficients dictionary
            inline const dictionary& coeffs() const;

            //- Return the solution dictionary
            inline const dictionary& solution() const;

            //- Return const access to the output properties dictionary
            inline const IOdictionary& outputProperties() const;

            //- Return output properties dictionary
            inline IOdictionary& outputProperties();


            // Region geometry

                //- Return the patch normal vectors
                const volVectorField& nHat() const;

                //- Return the face area magnitudes [m^2]
                const volScalarField::Internal& magSf() const;

                //- Return the cell layer volume/area [m]
                const volScalarField& VbyA() const;


            // Addressing

                //- Return true if patchi on the local region is a coupled
                //  patch to the primary region
                inline bool isCoupledPatch(const label regionPatchi) const;

                //- Return true if patchi on the primary region is a coupled
                //  patch to the local region
                inline bool isRegionPatch(const label primaryPatchi) const;

                //- Return the list of patch IDs on the primary region coupled
                //  to this region
                inline const labelList& primaryPatchIDs() const;

                //- Return the list of patch IDs internally coupled with the
                //  primary region
                inline const labelList& intCoupledPatchIDs() const;

                //- Return region ID corresponding to primaryPatchID
                inline label regionPatchID(const label primaryPatchID) const;

                //- Return the list of patch IDs opposite to internally
                //  coupled patches
                const labelList& passivePatchIDs() const;


        // Patch type information

            //- Return boundary types for mapped field patches
            //  Also maps internal field value
            //  Mapping region prescribed by underlying mapped poly patch
            template<class Type>
            wordList mappedFieldAndInternalPatchTypes() const;


        // Helper

            //- Return the coupled patch ID paired with coupled patch
            //  regionPatchi
            label nbrCoupledPatchID
            (
                const singleLayerRegionModel& nbrRegion,
                const label regionPatchi
            ) const;

            //- Convert a local region field to the primary region
            template<class Type>
            void toPrimary
            (
                const label regionPatchi,
                Field<Type>& regionField
            ) const;

            //- Convert a primary region field to the local region
            template<class Type>
            void toRegion
            (
                const label regionPatchi,
                Field<Type>& primaryFieldField
            ) const;

            //- Return a primary patch field mapped the local region
            template<class Type>
            void toRegion
            (
                Field<Type>& rf,
                const typename GeometricField
                <
                    Type,
                    fvPatchField,
                    volMesh
                >::Boundary& pBf
            ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const singleLayerRegionModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "singleLayerRegionModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "singleLayerRegionModelTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

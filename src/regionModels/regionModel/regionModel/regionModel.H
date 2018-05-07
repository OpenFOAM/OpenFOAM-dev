/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::regionModels::regionModel

Description
    Base class for region models

SourceFiles
    regionModelI.H
    regionModel.C

\*---------------------------------------------------------------------------*/

#ifndef regionModel_H
#define regionModel_H

#include "IOdictionary.H"
#include "Switch.H"
#include "labelList.H"
#include "volFields.H"
#include "mappedPatchBase.H"
#include "regionModelFunctionObjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace regionModels
{

/*---------------------------------------------------------------------------*\
                         Class regionModel Declaration
\*---------------------------------------------------------------------------*/

class regionModel
:
    public IOdictionary
{

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        regionModel(const regionModel&);

        //- Disallow default bitwise assignment
        void operator=(const regionModel&);

        //- Construct region mesh and fields
        void constructMeshObjects();

        //- Initialise the region
        void initialise();


protected:

    // Protected data

        //- Reference to the primary mesh database
        const fvMesh& primaryMesh_;

        //- Reference to the time database
        const Time& time_;

        //- Active flag
        Switch active_;

        //- Active information output
        Switch infoOutput_;

        //- Model name
        const word modelName_;

        //- Pointer to the region mesh database
        autoPtr<fvMesh> regionMeshPtr_;

        //- Model coefficients dictionary
        dictionary coeffs_;

        //- Dictionary of output properties
        autoPtr<IOdictionary> outputPropertiesPtr_;


        // Addressing

            //- List of patch IDs on the primary region coupled to this region
            labelList primaryPatchIDs_;

            //- List of patch IDs internally coupled with the primary region
            labelList intCoupledPatchIDs_;


        //- Region name
        word regionName_;

        //- Region model function objects
        regionModelFunctionObjectList functions_;


        // Inter-region AMI interpolation caching

            //- List of region names this region is coupled to
            mutable wordList interRegionAMINames_;

            //- List of AMI objects per coupled region
            mutable PtrList<PtrList<AMIPatchToPatchInterpolation>>
                interRegionAMI_;


    // Protected member functions

        //- Read control parameters from dictionary
        virtual bool read();

        //- Read control parameters from dictionary
        virtual bool read(const dictionary& dict);

        //- Create or return a new inter-region AMI object
        virtual const AMIPatchToPatchInterpolation& interRegionAMI
        (
            const regionModel& nbrRegion,
            const label regionPatchi,
            const label nbrPatchi,
            const bool flip
        ) const;


public:

    //- Runtime type information
    TypeName("regionModel");


    // Constructors

        //- Construct null
        regionModel(const fvMesh& mesh, const word& regionType);

        //- Construct from mesh, region type and name
        regionModel
        (
            const fvMesh& mesh,
            const word& regionType,
            const word& modelName,
            bool readFields = true
        );

        //- Construct from mesh and name and dict
        regionModel
        (
            const fvMesh& mesh,
            const word& regionType,
            const word& modelName,
            const dictionary& dict,
            bool readFields = true
        );


    //- Destructor
    virtual ~regionModel();


    // Member Functions

        // Access

            //- Return the reference to the primary mesh database
            inline const fvMesh& primaryMesh() const;

            //- Return the reference to the time database
            inline const Time& time() const;

            //- Return the active flag
            inline const Switch& active() const;

            //- Return the information flag
            inline const Switch& infoOutput() const;

            //- Return the model name
            inline const word& modelName() const;

            //- Return the region mesh database
            inline const fvMesh& regionMesh() const;

            //- Return the region mesh database for manipulation
            inline fvMesh& regionMesh();

            //- Return the model coefficients dictionary
            inline const dictionary& coeffs() const;

            //- Return the solution dictionary
            inline const dictionary& solution() const;

            //- Return const access to the output properties dictionary
            inline const IOdictionary& outputProperties() const;

            //- Return output properties dictionary
            inline IOdictionary& outputProperties();


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


        // Helper

            //- Return the coupled patch ID paired with coupled patch
            //  regionPatchi
            label nbrCoupledPatchID
            (
                const regionModel& nbrRegion,
                const label regionPatchi
            ) const;

            //- Map patch field from another region model to local patch
            template<class Type>
            tmp<Foam::Field<Type>> mapRegionPatchField
            (
                const regionModel& nbrRegion,
                const label regionPatchi,
                const label nbrPatchi,
                const Field<Type>& nbrField,
                const bool flip = false
            ) const;

            //- Map patch field from another region model to local patch
            template<class Type>
            tmp<Field<Type>> mapRegionPatchField
            (
                const regionModel& nbrRegion,
                const word& fieldName,
                const label regionPatchi,
                const bool flip = false
            ) const;

            //- Map patch internal field from another region model to local
            //  patch
            template<class Type>
            tmp<Field<Type>> mapRegionPatchInternalField
            (
                const regionModel& nbrRegion,
                const word& fieldName,
                const label regionPatchi,
                const bool flip = false
            ) const;

            //- Convert a local region field to the primary region
            template<class Type>
            void toPrimary
            (
                const label regionPatchi,
                List<Type>& regionField
            ) const;

            //- Convert a primary region field to the local region
            template<class Type>
            void toRegion
            (
                const label regionPatchi,
                List<Type>& primaryFieldField
            ) const;

            //- Convert a local region field to the primary region with op
            template<class Type, class CombineOp>
            void toPrimary
            (
                const label regionPatchi,
                List<Type>& regionField,
                const CombineOp& cop
            ) const;

            //- Convert a primary region field to the local region with op
            template<class Type, class CombineOp>
            void toRegion
            (
                const label regionPatchi,
                List<Type>& primaryFieldField,
                const CombineOp& cop
            ) const;


        // Evolution

            //- Main driver routing to evolve the region - calls other evolves
            virtual void evolve();

            //- Pre-evolve region
            virtual void preEvolveRegion();

            //- Evolve the region
            virtual void evolveRegion();

            //- Post-evolve region
            virtual void postEvolveRegion();


        // I-O

            //- Provide some feedback
            virtual void info();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "regionModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "regionModelTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

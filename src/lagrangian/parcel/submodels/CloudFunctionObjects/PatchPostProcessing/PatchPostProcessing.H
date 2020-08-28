/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::PatchPostProcessing

Description
    Standard post-processing

SourceFiles
    PatchPostProcessing.C

\*---------------------------------------------------------------------------*/

#ifndef PatchPostProcessing_H
#define PatchPostProcessing_H

#include "CloudFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class PatchPostProcessing Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class PatchPostProcessing
:
    public CloudFunctionObject<CloudType>
{
    // Private Data

        typedef typename CloudType::particleType parcelType;

        //- Maximum number of parcels to store - set as a scalar for I/O
        scalar maxStoredParcels_;

        //- List of patch indices to post-process
        labelList patchIDs_;

        //- List of time for each data record
        List<DynamicList<scalar>> times_;

        //- List of output data per patch
        List<DynamicList<string>> patchData_;


    // Private Member Functions

        //- Returns local patchi if patch is in patchIds_ list
        label applyToPatch(const label globalPatchi) const;


protected:

    // Protected Member Functions

        //- Write post-processing info
        void write();


public:

    //- Runtime type information
    TypeName("patchPostProcessing");


    // Constructors

        //- Construct from dictionary
        PatchPostProcessing
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        PatchPostProcessing(const PatchPostProcessing<CloudType>& ppm);

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObject<CloudType>> clone() const
        {
            return autoPtr<CloudFunctionObject<CloudType>>
            (
                new PatchPostProcessing<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~PatchPostProcessing();


    // Member Functions

        // Access

            //- Return maximum number of parcels to store per patch
            inline label maxStoredParcels() const;

            //- Return const mapping from local to global patch ids
            inline const labelList& patchIDs() const;


        // Evaluation

            //- Post-patch hook
            virtual void postPatch
            (
                const parcelType& p,
                const polyPatch& pp,
                bool& keepParticle
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PatchPostProcessingI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PatchPostProcessing.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

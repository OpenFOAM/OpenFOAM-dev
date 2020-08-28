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
    Foam::FieldActivatedInjection

Description
    Injection at specified positions, with the conditions:

    For injection to be allowed
      \verbatim
      factor*referenceField[celli] >= thresholdField[celli]
      \endverbatim
      where:
        - \c referenceField is the field used to supply the look-up values
        - \c thresholdField supplies the values beyond which the injection is
            permitted.

    Limited to a user-supplied number of injections per injector location

SourceFiles
    FieldActivatedInjection.C

\*---------------------------------------------------------------------------*/

#ifndef FieldActivatedInjection_H
#define FieldActivatedInjection_H

#include "InjectionModel.H"
#include "distributionModel.H"
#include "volFieldsFwd.H"
#include "GlobalIOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class FieldActivatedInjection Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class FieldActivatedInjection
:
    public InjectionModel<CloudType>
{
    // Private Data

        // Model parameters

            //- Factor to apply to reference field
            const scalar factor_;

            //- Reference field
            const volScalarField& referenceField_;

            //- Threshold field
            const volScalarField& thresholdField_;


        // Injector properties

            //- Name of file containing positions data
            const word positionsFile_;

            //- Field of injector (x,y,z) positions
            GlobalIOField<vector> positions_;

            //- List of cell labels corresponding to injector positions
            labelList injectorCells_;

            //- List of tetFace labels corresponding to injector positions
            labelList injectorTetFaces_;

            //- List of tetPt labels corresponding to injector positions
            labelList injectorTetPts_;

            //- Number of parcels per injector
            const label nParcelsPerInjector_;

            //- List of number of parcels injected for each injector
            labelList nParcelsInjected_;


        // Parcel properties

            //- Initial parcel velocity
            const vector U0_;

            //- List of parcel diameters
            scalarList diameters_;

            //- Parcel size distribution model
            const autoPtr<distributionModel>
                sizeDistribution_;


public:

    //- Runtime type information
    TypeName("fieldActivatedInjection");


    // Constructors

        //- Construct from dictionary
        FieldActivatedInjection
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        FieldActivatedInjection(const FieldActivatedInjection<CloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<InjectionModel<CloudType>> clone() const
        {
            return autoPtr<InjectionModel<CloudType>>
            (
                new FieldActivatedInjection<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~FieldActivatedInjection();


    // Member Functions

        //- Set injector locations when mesh is updated
        virtual void updateMesh();

        //- Return the end-of-injection time
        scalar timeEnd() const;

        //- Number of parcels to introduce relative to SOI
        virtual label parcelsToInject(const scalar time0, const scalar time1);

        //- Volume of parcels to introduce relative to SOI
        virtual scalar volumeToInject(const scalar time0, const scalar time1);


        // Injection geometry

            //- Set the injection position and owner cell, tetFace and tetPt
            virtual void setPositionAndCell
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                vector& position,
                label& cellOwner,
                label& tetFacei,
                label& tetPti
            );

            //- Set the parcel properties
            virtual void setProperties
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                typename CloudType::parcelType& parcel
            );

            //- Flag to identify whether model fully describes the parcel
            virtual bool fullyDescribed() const;

            //- Return flag to identify whether or not injection of parcelI is
            //  permitted
            virtual bool validInjection(const label parcelI);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "FieldActivatedInjection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

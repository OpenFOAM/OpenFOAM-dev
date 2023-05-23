/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::CellZoneInjection

Description
    Injection positions specified by a particle number density within a cell
    set.

    User specifies:
      - Number density of particles in cell set (effective)
      - Initial parcel velocity

    Properties:
      - Parcel diameters obtained by PDF model
      - All parcels introduced at SOI

SourceFiles
    CellZoneInjection.C

\*---------------------------------------------------------------------------*/

#ifndef CellZoneInjection_H
#define CellZoneInjection_H

#include "InjectionModel.H"
#include "distribution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class CellZoneInjection Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class CellZoneInjection
:
    public InjectionModel<CloudType>
{
    // Private Data

        //- Name of cell zone
        const word cellZoneName_;

        //- Mass to inject
        const scalar massTotal_;

        //- Number density
        const scalar numberDensity_;

        //- Field of parcel coordinates
        List<barycentric> injectorCoordinates_;

        //- List of cell labels corresponding to injector coordinates
        labelList injectorCells_;

        //- List of tetFace labels corresponding to injector coordinates
        labelList injectorTetFaces_;

        //- List of tetPt labels corresponding to injector coordinates
        labelList injectorTetPts_;

        //- Field of parcel diameters
        scalarList diameters_;

        //- Initial parcel velocity
        const vector U0_;

        //- Parcel size distribution model
        const autoPtr<distribution> sizeDistribution_;


    // Private Member Functions

        //- Set the parcel injection positions
        void setPositions(const labelList& cellZoneCells);


public:

    //- Runtime type information
    TypeName("cellZoneInjection");


    // Constructors

        //- Construct from dictionary
        CellZoneInjection
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        CellZoneInjection(const CellZoneInjection<CloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<InjectionModel<CloudType>> clone() const
        {
            return autoPtr<InjectionModel<CloudType>>
            (
                new CellZoneInjection<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~CellZoneInjection();


    // Member Functions

        //- Set injector locations when mesh is updated
        virtual void topoChange();

        //- Return the end-of-injection time
        virtual scalar timeEnd() const;

        //- Number of parcels to introduce relative to SOI
        virtual label nParcelsToInject(const scalar time0, const scalar time1);

        //- Parcel mass to introduce relative to SOI
        virtual scalar massToInject(const scalar time0, const scalar time1);


        // Injection geometry

            //- Set the injection position and owner cell, tetFace and tetPt
            virtual void setPositionAndCell
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                barycentric& coordinates,
                label& celli,
                label& tetFacei,
                label& tetPti,
                label& facei
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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CellZoneInjection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

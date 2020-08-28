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
    Foam::HeatTransferModel

Description
    Templated heat transfer model class

SourceFiles
    HeatTransferModel.C
    HeatTransferModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef HeatTransferModel_H
#define HeatTransferModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "CloudSubModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class HeatTransferModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class HeatTransferModel
:
    public CloudSubModelBase<CloudType>
{
    // Private Data

        //- Apply Bird's correction to the htc
        const Switch BirdCorrection_;


public:

    //- Runtime type information
    TypeName("heatTransferModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        HeatTransferModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        HeatTransferModel(CloudType& owner);

        //- Construct from dictionary
        HeatTransferModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );

        //- Construct copy
        HeatTransferModel(const HeatTransferModel<CloudType>& htm);

        //- Construct and return a clone
        virtual autoPtr<HeatTransferModel<CloudType>> clone() const = 0;


    //- Destructor
    virtual ~HeatTransferModel();


    //- Selector
    static autoPtr<HeatTransferModel<CloudType>> New
    (
        const dictionary& dict,
        CloudType& owner
    );


    // Member Functions

        // Access

            //- Return the Bird htc correction flag
            const Switch& BirdCorrection() const;


        // Evaluation

            //- Nusselt number
            virtual scalar Nu
            (
                const scalar Re,
                const scalar Pr
            ) const = 0;

            //- Return heat transfer coefficient
            virtual scalar htc
            (
                const scalar dp,
                const scalar Re,
                const scalar Pr,
                const scalar kappa,
                const scalar NCpW
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeHeatTransferModel(CloudType)                                       \
                                                                               \
    typedef Foam::CloudType::thermoCloudType thermoCloudType;                  \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        Foam::HeatTransferModel<thermoCloudType>,                              \
        0                                                                      \
    );                                                                         \
    namespace Foam                                                             \
    {                                                                          \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            HeatTransferModel<thermoCloudType>,                                \
            dictionary                                                         \
        );                                                                     \
    }


#define makeHeatTransferModelType(SS, CloudType)                               \
                                                                               \
    typedef Foam::CloudType::thermoCloudType thermoCloudType;                  \
    defineNamedTemplateTypeNameAndDebug(Foam::SS<thermoCloudType>, 0);         \
                                                                               \
    Foam::HeatTransferModel<thermoCloudType>::                                 \
        adddictionaryConstructorToTable<Foam::SS<thermoCloudType>>             \
            add##SS##CloudType##thermoCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "HeatTransferModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

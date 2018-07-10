/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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
    Foam::IsotropyModel

Description
    Base class for collisional return-to-isotropy models.

SourceFiles
    IsotropyModel.C
    IsotropyModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef IsotropyModel_H
#define IsotropyModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "CloudSubModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

class TimeScaleModel;

/*---------------------------------------------------------------------------*\
                        Class IsotropyModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class IsotropyModel
:
    public CloudSubModelBase<CloudType>
{
protected:

        //- Time scale model
        autoPtr<TimeScaleModel> timeScaleModel_;


public:

    //- Runtime type information
    TypeName("isotropyModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        IsotropyModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        IsotropyModel(CloudType& owner);

        //- Construct from components
        IsotropyModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );

        //- Construct a copy
        IsotropyModel(const IsotropyModel<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<IsotropyModel<CloudType>> clone() const = 0;


    //- Destructor
    virtual ~IsotropyModel();


    //- Selector
    static autoPtr<IsotropyModel<CloudType>> New
    (
        const dictionary& dict,
        CloudType& owner
    );


    //- Member Functions

        //- Calculate velocities
        virtual void calculate() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeIsotropyModel(CloudType)                                           \
                                                                               \
    typedef Foam::CloudType::MPPICCloudType MPPICCloudType;                    \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        Foam::IsotropyModel<MPPICCloudType>,                                   \
        0                                                                      \
    );                                                                         \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            IsotropyModel<MPPICCloudType>,                                     \
            dictionary                                                         \
        );                                                                     \
    }


#define makeIsotropyModelType(SS, CloudType)                                   \
                                                                               \
    typedef Foam::CloudType::MPPICCloudType MPPICCloudType;                    \
    defineNamedTemplateTypeNameAndDebug                                        \
        (Foam::IsotropyModels::SS<MPPICCloudType>, 0);                         \
                                                                               \
    Foam::IsotropyModel<MPPICCloudType>::                                      \
        adddictionaryConstructorToTable                                        \
        <Foam::IsotropyModels::SS<MPPICCloudType>>                             \
            add##SS##CloudType##MPPICCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "IsotropyModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

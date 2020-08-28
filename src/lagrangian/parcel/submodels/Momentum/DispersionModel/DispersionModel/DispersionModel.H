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
    Foam::DispersionModel

Description

\*---------------------------------------------------------------------------*/

#ifndef DispersionModel_H
#define DispersionModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "CloudSubModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class DispersionModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class DispersionModel
:
    public CloudSubModelBase<CloudType>
{
public:

    //- Runtime type information
    TypeName("dispersionModel");


    // Declare runtime constructor selection table

    declareRunTimeSelectionTable
    (
        autoPtr,
        DispersionModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        DispersionModel(CloudType& owner);

        //- Construct from components
        DispersionModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );

        //- Construct copy
        DispersionModel(const DispersionModel<CloudType>& dm);

        //- Construct and return a clone
        virtual autoPtr<DispersionModel<CloudType>> clone() const = 0;


    //- Destructor
    virtual ~DispersionModel();


    //- Selector
    static autoPtr<DispersionModel<CloudType>> New
    (
        const dictionary& dict,
        CloudType& owner
    );


    // Member Functions

        //- Update (disperse particles)
        virtual vector update
        (
            const scalar dt,
            const label celli,
            const vector& U,
            const vector& Uc,
            vector& UTurb,
            scalar& tTurb
        ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeDispersionModel(CloudType)                                         \
                                                                               \
    typedef Foam::CloudType::momentumCloudType momentumCloudType;            \
    defineTemplateTypeNameAndDebug                                             \
        (Foam::DispersionModel<momentumCloudType>, 0);                        \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            DispersionModel<momentumCloudType>,                               \
            dictionary                                                         \
        );                                                                     \
    }

#define makeDispersionModelType(SS, CloudType)                                 \
                                                                               \
    typedef Foam::CloudType::momentumCloudType momentumCloudType;            \
    defineNamedTemplateTypeNameAndDebug(Foam::SS<momentumCloudType>, 0);      \
                                                                               \
    Foam::DispersionModel<momentumCloudType>::                                \
        adddictionaryConstructorToTable<Foam::SS<momentumCloudType>>          \
            add##SS##CloudType##momentumCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DispersionModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

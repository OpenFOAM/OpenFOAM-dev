/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    Foam::AtomisationModel

Description
    Templated atomisation model class

SourceFiles
    AtomisationModel.C
    AtomisationModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef AtomisationModel_H
#define AtomisationModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "CloudSubModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class AtomisationModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class AtomisationModel
:
    public CloudSubModelBase<CloudType>
{
public:

    //- Runtime type information
    TypeName("atomisationModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        AtomisationModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        AtomisationModel(CloudType& owner);

        //- Construct from dictionary
        AtomisationModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );

        //- Construct copy
        AtomisationModel(const AtomisationModel<CloudType>& am);

        //- Construct and return a clone
        virtual autoPtr<AtomisationModel<CloudType>> clone() const = 0;


    //- Destructor
    virtual ~AtomisationModel();


    //- Selector
    static autoPtr<AtomisationModel<CloudType>> New
    (
        const dictionary& dict,
        CloudType& owner
    );


    // Member Functions

        //- Average temperature calculation
        scalar Taverage(const scalar& Tliq, const scalar& Tc) const;

        //- Initial value of liquidCore
        virtual scalar initLiquidCore() const = 0;

        //- Flag to indicate if chi needs to be calculated
        virtual bool calcChi() const = 0;

        virtual void update
        (
            const scalar dt,
            scalar& d,
            scalar& liquidCore,
            scalar& tc,
            const scalar rho,
            const scalar mu,
            const scalar sigma,
            const scalar volFlowRate,
            const scalar rhoAv,
            const scalar Urel,
            const vector& pos,
            const vector& injectionPos,
            const scalar pAmbient,
            const scalar chi,
            Random& rndGen
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeAtomisationModel(CloudType)                                        \
                                                                               \
    typedef Foam::CloudType::sprayCloudType sprayCloudType;                    \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        Foam::AtomisationModel<sprayCloudType>,                                \
        0                                                                      \
    );                                                                         \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            AtomisationModel<sprayCloudType>,                                  \
            dictionary                                                         \
        );                                                                     \
    }


#define makeAtomisationModelType(SS, CloudType)                                \
                                                                               \
    typedef Foam::CloudType::sprayCloudType sprayCloudType;                    \
    defineNamedTemplateTypeNameAndDebug(Foam::SS<sprayCloudType>, 0);          \
                                                                               \
    Foam::AtomisationModel<sprayCloudType>::                                   \
        adddictionaryConstructorToTable<Foam::SS<sprayCloudType>>              \
            add##SS##CloudType##sprayCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "AtomisationModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

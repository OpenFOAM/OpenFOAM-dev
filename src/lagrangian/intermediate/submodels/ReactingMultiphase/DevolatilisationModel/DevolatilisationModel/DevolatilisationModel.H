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
    Foam::DevolatilisationModel

Description
    Templated devolatilisation model class

SourceFiles
    DevolatilisationModel.C
    DevolatilisationModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef DevolatilisationModel_H
#define DevolatilisationModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "CloudSubModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class DevolatilisationModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class DevolatilisationModel
:
    public CloudSubModelBase<CloudType>
{
protected:

    // Protected data

        //- Mass of lagrangian phase converted
        scalar dMass_;


public:

    //- Runtime type information
    TypeName("devolatilisationModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        DevolatilisationModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        DevolatilisationModel(CloudType& owner);

        //- Construct from dictionary
        DevolatilisationModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );

        //- Construct copy
        DevolatilisationModel(const DevolatilisationModel<CloudType>& dm);

        //- Construct and return a clone
        virtual autoPtr<DevolatilisationModel<CloudType> > clone() const
        {
            return autoPtr<DevolatilisationModel<CloudType> >
            (
                new DevolatilisationModel<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~DevolatilisationModel();


    //- Selector
    static autoPtr<DevolatilisationModel<CloudType> > New
    (
        const dictionary& dict,
        CloudType& owner
    );


    // Member Functions

        //- Update model
        virtual void calculate
        (
            const scalar dt,
            const scalar age,
            const scalar mass0,
            const scalar mass,
            const scalar T,
            const scalarField& YGasEff,
            const scalarField& YLiquidEff,
            const scalarField& YSolidEff,
            label& canCombust,
            scalarField& dMassDV
        ) const;

        //- Add to devolatilisation mass
        void addToDevolatilisationMass(const scalar dMass);


        // I-O

            //- Write injection info to stream
            virtual void info(Ostream& os);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeDevolatilisationModel(CloudType)                                  \
                                                                              \
    typedef CloudType::reactingMultiphaseCloudType                            \
        reactingMultiphaseCloudType;                                          \
    defineNamedTemplateTypeNameAndDebug                                       \
    (                                                                         \
        DevolatilisationModel<reactingMultiphaseCloudType>,                   \
        0                                                                     \
    );                                                                        \
    defineTemplateRunTimeSelectionTable                                       \
    (                                                                         \
        DevolatilisationModel<reactingMultiphaseCloudType>,                   \
        dictionary                                                            \
    );


#define makeDevolatilisationModelType(SS, CloudType)                          \
                                                                              \
    typedef CloudType::reactingMultiphaseCloudType                            \
        reactingMultiphaseCloudType;                                          \
    defineNamedTemplateTypeNameAndDebug(SS<reactingMultiphaseCloudType>, 0);  \
                                                                              \
    DevolatilisationModel<reactingMultiphaseCloudType>::                      \
        adddictionaryConstructorToTable<SS<reactingMultiphaseCloudType> >     \
        add##SS##CloudType##reactingMultiphaseCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "DevolatilisationModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

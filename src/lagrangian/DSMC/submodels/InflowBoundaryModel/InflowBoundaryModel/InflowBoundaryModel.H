/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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
    Foam::InflowBoundaryModel


Description
    Templated inflow boundary model class

SourceFiles
    InflowBoundaryModel.C
    InflowBoundaryModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef InflowBoundaryModel_H
#define InflowBoundaryModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class InflowBoundaryModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class InflowBoundaryModel
{
    // Private data

        //- The cloud dictionary
        const dictionary& dict_;

        // Reference to the owner cloud class
        CloudType& owner_;

        //- The coefficients dictionary
        const dictionary coeffDict_;


public:

    //- Runtime type information
    TypeName("InflowBoundaryModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        InflowBoundaryModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        InflowBoundaryModel(CloudType& owner);

        //- Construct from dictionary
        InflowBoundaryModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );


    //- Destructor
    virtual ~InflowBoundaryModel();


    //- Selector
    static autoPtr<InflowBoundaryModel<CloudType>> New
    (
        const dictionary& dict,
        CloudType& owner
    );


    // Access

        //- Return const access the owner cloud object
        inline const CloudType& owner() const;

        //- Return non-const access the owner cloud object for manipulation
        inline CloudType& owner();

        //- Return the owner cloud dictionary
        inline const dictionary& dict() const;

        //- Return the coefficients dictionary
        inline const dictionary& coeffDict() const;

    // Mapping

        //- Remap the particles to the correct cells following mesh change
        virtual void autoMap(const mapPolyMesh&)
        {}

    //- Introduce particles
    virtual void inflow() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeInflowBoundaryModel(CloudType)                                     \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(InflowBoundaryModel<CloudType>, 0);    \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        InflowBoundaryModel<CloudType>,                                        \
        dictionary                                                             \
    );


#define makeInflowBoundaryModelType(SS, CloudType)                             \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(SS<CloudType>, 0);                     \
                                                                               \
    InflowBoundaryModel<CloudType>::                                           \
        adddictionaryConstructorToTable<SS<CloudType>>                         \
            add##SS##CloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "InflowBoundaryModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::CompositionModel

Description
    Templated reacting parcel composition model class
    Consists of carrier species (via thermo package), and additional liquids
    and solids

SourceFiles
    CompositionModel.C
    CompositionModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef CompositionModel_H
#define CompositionModel_H

#include "CloudSubModelBase.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

#include "PtrList.H"
#include "SLGThermo.H"

#include "phasePropertiesList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class CompositionModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class CompositionModel
:
    public CloudSubModelBase<CloudType>
{
    // Private data

        //- Reference to the thermo database
        const SLGThermo& thermo_;

        //- List of phase properties
        phasePropertiesList phaseProps_;


public:

    //- Runtime type information
    TypeName("compositionModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        CompositionModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        CompositionModel(CloudType& owner);

        //- Construct from dictionary
        CompositionModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );

        //- Construct copy
        CompositionModel(const CompositionModel<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<CompositionModel<CloudType>> clone() const = 0;


    //- Destructor
    virtual ~CompositionModel();


    //- Selector
    static autoPtr<CompositionModel<CloudType>> New
    (
        const dictionary& dict,
        CloudType& owner
    );


    // Member Functions

        // Access

            //- Return the thermo database
            const SLGThermo& thermo() const;


            // Composition lists

                //- Return the carrier components (wrapper function)
                const basicSpecieMixture& carrier() const;

                //- Return the global (additional) liquids
                const liquidMixtureProperties& liquids() const;

                //- Return the global (additional) solids
                const solidMixtureProperties& solids() const;

                //- Return the list of phase properties
                const phasePropertiesList& phaseProps() const;

                //- Return the number of phases
                label nPhase() const;


            // Phase properties

                //- Return the list of phase type names
                //  If only 1 phase, return the component names of that phase
                const wordList& phaseTypes() const;

                //- Return the list of state labels (s), (l), (g) etc.
                const wordList& stateLabels() const;

                //- Return the list of component names for phaseI
                const wordList& componentNames(const label phaseI) const;

                //- Return global id of component cmptName in carrier thermo
                label carrierId
                (
                    const word& cmptName,
                    const bool allowNotFound = false
                ) const;

                //- Return local id of component cmptName in phase phaseI
                label localId
                (
                    const label phaseI,
                    const word& cmptName,
                    const bool allowNotFound = false
                ) const;

                //- Return carrier id of component given local id
                label localToCarrierId
                (
                    const label phaseI,
                    const label id,
                    const bool allowNotFound = false
                ) const;

                //- Return the list of phase phaseI mass fractions
                const scalarField& Y0(const label phaseI) const;

                //- Return the list of phase phaseI volume fractions fractions
                //  based on supplied mass fractions Y
                scalarField X(const label phaseI, const scalarField& Y) const;


            // Mixture properties

                //- Return the list of mixture mass fractions
                //  If only 1 phase, return component fractions of that phase
                virtual const scalarField& YMixture0() const = 0;

                // Indices of gas, liquid and solid phases in phase properties
                // list - returns -1 if not applicable

                    //- Gas id
                    virtual label idGas() const = 0;

                    //- Liquid id
                    virtual label idLiquid() const = 0;

                    //- Solid id
                    virtual label idSolid() const = 0;


        // Evaluation

            //- Return total enthalpy for the phase phaseI
            virtual scalar H
            (
                const label phaseI,
                const scalarField& Y,
                const scalar p,
                const scalar T
            ) const;

            //- Return sensible enthalpy for the phase phaseI
            virtual scalar Hs
            (
                const label phaseI,
                const scalarField& Y,
                const scalar p,
                const scalar T
            ) const;

            //- Return chemical enthalpy for the phase phaseI
            virtual scalar Hc
            (
                const label phaseI,
                const scalarField& Y,
                const scalar p,
                const scalar T
            ) const;

            //- Return specific heat caoacity for the phase phaseI
            virtual scalar Cp
            (
                const label phaseI,
                const scalarField& Y,
                const scalar p,
                const scalar T
            ) const;

            //- Return latent heat for the phase phaseI
            virtual scalar L
            (
                const label phaseI,
                const scalarField& Y,
                const scalar p,
                const scalar T
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeCompositionModel(CloudType)                                        \
                                                                               \
    typedef Foam::CloudType::reactingCloudType reactingCloudType;              \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        Foam::CompositionModel<reactingCloudType>,                             \
        0                                                                      \
    );                                                                         \
    namespace Foam                                                             \
    {                                                                          \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            CompositionModel<reactingCloudType>,                               \
            dictionary                                                         \
        );                                                                     \
    }


#define makeCompositionModelType(SS, CloudType)                                \
                                                                               \
    typedef Foam::CloudType::reactingCloudType reactingCloudType;              \
    defineNamedTemplateTypeNameAndDebug(Foam::SS<reactingCloudType>, 0);       \
                                                                               \
    Foam::CompositionModel<reactingCloudType>::                                \
        adddictionaryConstructorToTable<Foam::SS<reactingCloudType>>           \
            add##SS##CloudType##reactingCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CompositionModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

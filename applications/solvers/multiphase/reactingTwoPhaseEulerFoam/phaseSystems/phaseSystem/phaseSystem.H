/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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
    Foam::phaseSystem

Description
    Class to represent a system of phases and model interfacial transfers
    between them.

SourceFiles
    phaseSystem.C

\*---------------------------------------------------------------------------*/

#ifndef phaseSystem_H
#define phaseSystem_H

#include "IOdictionary.H"

#include "phaseModel.H"
#include "phasePair.H"
#include "orderedPhasePair.H"
#include "HashPtrTable.H"

#include "IOMRFZoneList.H"
#include "fvIOoptionList.H"

#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatricesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class blendingMethod;
template <class modelType> class BlendedInterfacialModel;
class surfaceTensionModel;
class aspectRatioModel;

/*---------------------------------------------------------------------------*\
                         Class phaseSystem Declaration
\*---------------------------------------------------------------------------*/

class phaseSystem
:
    public IOdictionary
{
public:

    // Public typedefs

        typedef
            HashPtrTable
            <
                volScalarField,
                phasePairKey,
                phasePairKey::hash
            >
            KdTable;

        typedef
            HashPtrTable
            <
                volScalarField,
                phasePairKey,
                phasePairKey::hash
            >
            VmTable;

        typedef
            HashPtrTable
            <
                fvVectorMatrix,
                word,
                string::hash
            >
            momentumTransferTable;

        typedef
            HashPtrTable
            <
                fvScalarMatrix,
                word,
                string::hash
            >
            heatTransferTable;

        typedef
            HashPtrTable
            <
                fvScalarMatrix,
                word,
                string::hash
            >
            massTransferTable;


protected:

    // Protected typedefs

        typedef
            HashTable<autoPtr<phaseModel>, word, word::hash>
            phaseModelTable;

        typedef
            HashTable<dictionary, phasePairKey, phasePairKey::hash>
            dictTable;

        typedef
            HashTable<autoPtr<phasePair>, phasePairKey, phasePairKey::hash>
            phasePairTable;

        typedef
            HashTable<autoPtr<blendingMethod>, word, word::hash>
            blendingMethodTable;

        typedef
            HashTable
            <
                autoPtr<surfaceTensionModel>,
                phasePairKey,
                phasePairKey::hash
            >
            surfaceTensionModelTable;

        typedef
            HashTable
            <
                autoPtr<aspectRatioModel>,
                phasePairKey,
                phasePairKey::hash
            >
            aspectRatioModelTable;


    // Protected data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Phase names
        wordList phaseNames_;

        //- Phase models
        phaseModelTable phaseModels_;

        //- Phase pairs
        phasePairTable phasePairs_;

        //- Total volumetric flux
        surfaceScalarField phi_;

        //- Rate of change of pressure
        volScalarField dpdt_;

        //- Optional MRF zones
        IOMRFZoneList MRF_;

        //- Optional FV-options
        mutable fv::IOoptionList fvOptions_;

        //- Blending methods
        blendingMethodTable blendingMethods_;


        // Sub Models

            //- Surface tension models
            surfaceTensionModelTable surfaceTensionModels_;

            //- Aspect ratio models
            aspectRatioModelTable aspectRatioModels_;


    // Protected member functions

        //- Generate the phases
        HashTable<autoPtr<phaseModel> > generatePhaseModels
        (
            const wordList& names
        ) const;

        //- Generate the mixture flux
        tmp<surfaceScalarField> generatePhi
        (
            const HashTable<autoPtr<phaseModel> >& phaseModels
        ) const;

        //- Generate pairs
        void generatePairs
        (
            const dictTable& modelDicts
        );

        //- Generate pairs and sub-model tables
        template<class modelType>
        void createSubModels
        (
            const dictTable& modelDicts,
            HashTable
            <
                autoPtr<modelType>,
                phasePairKey,
                phasePairKey::hash
            >& models
        );

        //- Generate pairs and sub-model tables
        template<class modelType>
        void generatePairsAndSubModels
        (
            const word& modelName,
            HashTable
            <
                autoPtr<modelType>,
                phasePairKey,
                phasePairKey::hash
            >& models
        );

        //- Generate pairs and blended sub-model tables
        template<class modelType>
        void generatePairsAndSubModels
        (
            const word& modelName,
            HashTable
            <
                autoPtr<BlendedInterfacialModel<modelType> >,
                phasePairKey,
                phasePairKey::hash
            >& models
        );

        //- Generate pairs and per-phase sub-model tables
        template<class modelType>
        void generatePairsAndSubModels
        (
            const word& modelName,
            HashTable
            <
                HashTable<autoPtr<modelType> >,
                phasePairKey,
                phasePairKey::hash
            >& models
        );


public:

    //- Runtime type information
    TypeName("phaseSystem");

    //- Default name of the phase properties dictionary
    static const word propertiesName;


    // Constructors

        //- Construct from fvMesh
        phaseSystem(const fvMesh& mesh);


    //- Destructor
    virtual ~phaseSystem();


    // Member Functions

        //- Return the mixture density
        tmp<volScalarField> rho() const;

        //- Return the mixture velocity
        tmp<volVectorField> U() const;

        //- Return the surface tension coefficient
        tmp<volScalarField> sigma(const phasePairKey& key) const;

        //- Return the drag coefficient
        virtual tmp<volScalarField> Kd(const phasePairKey& key) const = 0;

        //- Return the face drag coefficient
        virtual tmp<surfaceScalarField> Kdf(const phasePairKey& key) const = 0;

        //- Return the virtual mass coefficient
        virtual tmp<volScalarField> Vm(const phasePairKey& key) const = 0;

        //- Return the face virtual mass coefficient
        virtual tmp<surfaceScalarField> Vmf(const phasePairKey& key) const = 0;

        //- Return the combined force (lift + wall-lubrication)
        virtual tmp<volVectorField> F(const phasePairKey& key) const = 0;

        //- Return the combined face-force (lift + wall-lubrication)
        virtual tmp<surfaceScalarField> Ff(const phasePairKey& key) const = 0;

        //- Return the turbulent diffusivity
        //  Multiplies the phase-fraction gradient
        virtual tmp<volScalarField> D(const phasePairKey& key) const = 0;

        //- Return the interfacial mass flow rate
        virtual tmp<volScalarField> dmdt(const phasePairKey& key) const = 0;

        //- Return the momentum transfer matrices
        virtual autoPtr<momentumTransferTable> momentumTransfer() const = 0;

        //- Return the heat transfer matrices
        virtual autoPtr<heatTransferTable> heatTransfer() const = 0;

        //- Return the mass transfer matrices
        virtual autoPtr<massTransferTable> massTransfer() const = 0;

        //- Solve for the phase fractions
        virtual void solve();

        //- Correct the fluid properties other than the thermo and turbulence
        virtual void correct();

        //- Correct the kinematics
        virtual void correctKinematics();

        //- Correct the thermodynamics
        virtual void correctThermo();

        //- Correct the turbulence
        virtual void correctTurbulence();

        //- Read base phaseProperties dictionary
        virtual bool read();


        // Access

            //- Constant access the mesh
            inline const fvMesh& mesh() const;

            //- Constant access the mixture flux
            inline const surfaceScalarField& phi() const;

            //- Access the mixture flux
            inline surfaceScalarField& phi();

            //- Constant access the rate of change of the pressure
            inline const volScalarField& dpdt() const;

            //- Access the rate of change of the pressure
            inline volScalarField& dpdt();

            //- Return MRF zones
            inline const IOMRFZoneList& MRF() const;

            //- Optional FV-options
            inline fv::IOoptionList& fvOptions() const;

            //- Access a sub model between a phase pair
            template <class modelType>
            const modelType& lookupSubModel(const phasePair& key) const;

            //- Access a sub model between two phases
            template <class modelType>
            const modelType& lookupSubModel
            (
                const phaseModel& dispersed,
                const phaseModel& continuous
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "phaseSystemI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "phaseSystemTemplates.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

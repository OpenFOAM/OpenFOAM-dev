/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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
    Foam::radiation::mixtureFractionSoot

Description
    This soot model is purely an state model. The amount of soot produced is
    determined by a single step chemistry as :

        nuf Fuel + nuOx Ox = nuP P + nuSoot soot

    nuSoot is prescribed by the user.

    The single step chemistry used is read from the combustion.
    The soot is not considered into the thermodynamics of the system and it
    is not considered as an extra specie in the solver.

    The spacial distribution is given by the normalization of the first product
    on the rhs of the reaction by default or it can be added as input.

    The input dictionary reads like in the radiationProperties dictionary:

    sootModel mixtureFractionSoot<gasHThermoPhysics>;

    mixtureFractionSootCoeffs
    {
        nuSoot              0.015;
        Wsoot               12;
        mappingField        P;
    }

SourceFiles
    mixtureFractionSoot.C

\*---------------------------------------------------------------------------*/

#ifndef mixtureFractionSoot_H
#define mixtureFractionSoot_H

#include "interpolationLookUpTable.H"
#include "sootModel.H"
#include "HashTable.H"

#include "fluidThermo.H"
#include "singleStepReactingMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace radiation
{

/*---------------------------------------------------------------------------*\
                    Class mixtureFractionSoot Declaration
\*---------------------------------------------------------------------------*/
template<class ThermoType>
class mixtureFractionSoot
:
    public sootModel
{

    // Static functions

        //- Check mixture in thermo
        static const singleStepReactingMixture<ThermoType>& checkThermo
        (
            const fluidThermo&
        );


    // Private data

         //- Soot mass fraction
        volScalarField soot_;

        //- Soot model dictionary
        dictionary coeffsDict_;

        //- Soot yield
        scalar nuSoot_;

        //- Soot molecular weight
        scalar Wsoot_;

        //- Maximum soot mass concentration at stoichiometric
        scalar sootMax_;

        //- Name of the field mapping the soot
        word mappingFieldName_;

        //- Maximum value of the map field
        scalar mapFieldMax_;

        //- Thermo package
        const fluidThermo& thermo_;

        //- Auto Ptr to singleStepReactingMixture
        const singleStepReactingMixture<ThermoType>& mixture_;


public:

    //- Runtime type information
    TypeName("mixtureFractionSoot");


    // Constructors

        //- Construct from components
        mixtureFractionSoot
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const word& modelType
        );


    //- Destructor
    virtual ~mixtureFractionSoot();


    // Member Functions

        // Edit

            //- Main update/correction routine
            virtual void correct();


        // Access

            //- Return Ysoot
            const volScalarField& soot() const
            {
                return soot_;
            }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
} // End namespace radiation

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "mixtureFractionSoot.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

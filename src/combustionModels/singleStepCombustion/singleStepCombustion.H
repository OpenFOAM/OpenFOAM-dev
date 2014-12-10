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
    Foam::combustionModels::singleStepCombustion

Description
    Base class for combustion models using singleStepReactingMixture.

SourceFiles
    singleStepCombustion.C

\*---------------------------------------------------------------------------*/

#ifndef singleStepCombustion_H
#define singleStepCombustion_H

#include "singleStepReactingMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                    Class singleStepCombustion Declaration
\*---------------------------------------------------------------------------*/

template<class CombThermoType, class ThermoType>
class singleStepCombustion
:
    public CombThermoType
{
    // Private Member Functions

        //- Disallow copy construct
        singleStepCombustion(const singleStepCombustion&);

        //- Disallow default bitwise assignment
        void operator=(const singleStepCombustion&);


protected:

    // Protected data

        //- Pointer to singleStepReactingMixture mixture
        singleStepReactingMixture<ThermoType>* singleMixturePtr_;

        //- Fuel consumption rate
        volScalarField wFuel_;

        //- Semi-implicit (true) or explicit (false) treatment
        bool semiImplicit_;


public:

    // Constructors

        //- Construct from components
        singleStepCombustion(const word& modelType, const fvMesh& mesh);


    //- Destructor
    virtual ~singleStepCombustion();


    // Member Functions

        // Evolution

            //- Fuel consumption rate matrix
            virtual tmp<fvScalarMatrix> R(volScalarField& Y) const;

            //- Heat release rate calculated from fuel consumption rate matrix
            virtual tmp<volScalarField> dQ() const;

            //- Sensible enthalpy source term
            virtual tmp<volScalarField> Sh() const;


        // I-O

            //- Update properties from given dictionary
            virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "singleStepCombustion.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

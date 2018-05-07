/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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
    Foam::interfaceCompositionModels::Henry

Description
    Henry's law for gas solubiliy in liquid. The concentration of the dissolved
    species in the liquid is proportional to its partial pressure in the gas.
    The dimensionless constant of proportionality between concentrations on
    each side of the interface is \f$k\f$, and is given for each species.
    Mixing in the gas is assumed to be ideal.

SourceFiles
    Henry.C

\*---------------------------------------------------------------------------*/

#ifndef Henry_H
#define Henry_H

#include "InterfaceCompositionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace interfaceCompositionModels
{

/*---------------------------------------------------------------------------*\
                         Class Henry Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo, class OtherThermo>
class Henry
:
    public InterfaceCompositionModel<Thermo, OtherThermo>
{
    // Private data

        //- Solubility coefficients
        const scalarList k_;

        //- The remaining solvent species fraction
        volScalarField YSolvent_;


public:

    //- Runtime type information
    TypeName("Henry");

    // Constructors

        //- Construct from components
        Henry
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~Henry();


    // Member Functions

        //- Update the composition
        virtual void update(const volScalarField& Tf);

        //- The interface species fraction
        virtual tmp<volScalarField> Yf
        (
            const word& speciesName,
            const volScalarField& Tf
        ) const;

        //- The interface species fraction derivative w.r.t. temperature
        virtual tmp<volScalarField> YfPrime
        (
            const word& speciesName,
            const volScalarField& Tf
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace interfaceCompositionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Henry.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

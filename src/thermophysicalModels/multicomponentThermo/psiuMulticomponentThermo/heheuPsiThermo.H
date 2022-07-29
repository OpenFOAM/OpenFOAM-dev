/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::heheuMulticomponentThermo

Description
    Foam::heheuMulticomponentThermo

SourceFiles
    heheuMulticomponentThermo.C

\*---------------------------------------------------------------------------*/

#ifndef heheuPsiThermo_H
#define heheuPsiThermo_H

#include "heThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class heheuPsiThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BasicPsiThermo, class MixtureType>
class heheuPsiThermo
:
    public heThermo<BasicPsiThermo, MixtureType>
{
    // Private Data

        //- Unburnt temperature
        volScalarField Tu_;

        //- Unburnt energy
        volScalarField heu_;


    // Private Member Functions

        //- Calculate the thermo variables
        void calculate();


public:

    //- Runtime type information
    TypeName("heheuPsiThermo");


    // Constructors

        //- Construct from mesh and phase name
        heheuPsiThermo
        (
            const fvMesh&,
            const word& phaseName
        );

        //- Disallow default bitwise copy construction
        heheuPsiThermo
        (
            const heheuPsiThermo<BasicPsiThermo, MixtureType>&
        ) = delete;


    //- Destructor
    virtual ~heheuPsiThermo();


    // Member Functions

        //- Update properties
        virtual void correct();


        // Access to thermodynamic state variables.

            //- Unburnt gas enthalpy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& heu()
            {
                return heu_;
            }

            //- Unburnt gas enthalpy [J/kg]
            virtual const volScalarField& heu() const
            {
                return heu_;
            }

            //- Unburnt gas temperature [K]
            virtual const volScalarField& Tu() const
            {
                return Tu_;
            }


        // Fields derived from thermodynamic state variables

            //- Unburnt gas enthalpy for cell-set [J/kg]
            virtual tmp<scalarField> heu
            (
                const scalarField& T,
                const labelList& cells
            ) const;

            //- Unburnt gas enthalpy for patch [J/kg]
            virtual tmp<scalarField> heu
            (
                const scalarField& T,
                const label patchi
            ) const;


            //- Burnt gas temperature [K]
            virtual tmp<volScalarField> Tb() const;

            //- Unburnt gas compressibility [s^2/m^2]
            virtual tmp<volScalarField> psiu() const;

            //- Burnt gas compressibility [s^2/m^2]
            virtual tmp<volScalarField> psib() const;


        // Access to transport variables

            //- Dynamic viscosity of unburnt gas [kg/m/s]
            virtual tmp<volScalarField> muu() const;

            //- Dynamic viscosity of burnt gas [kg/m/s]
            virtual tmp<volScalarField> mub() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
    #include "heheuPsiThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

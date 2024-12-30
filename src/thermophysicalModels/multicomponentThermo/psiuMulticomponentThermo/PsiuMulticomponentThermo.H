/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::PsiuMulticomponentThermo

Description
    Thermo implementation based on compressibility with additional unburnt
    thermodynamic state

SourceFiles
    PsiuMulticomponentThermo.C

\*---------------------------------------------------------------------------*/

#ifndef PsiuMulticomponentThermo_H
#define PsiuMulticomponentThermo_H

#include "BasicThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class PsiuMulticomponentThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BaseThermo>
class PsiuMulticomponentThermo
:
    public BaseThermo
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

    // Constructors

        //- Construct from mesh and phase name
        PsiuMulticomponentThermo(const fvMesh&, const word& phaseName);

        //- Disallow default bitwise copy construction
        PsiuMulticomponentThermo
        (
            const PsiuMulticomponentThermo<BaseThermo>&
        ) = delete;


    //- Destructor
    virtual ~PsiuMulticomponentThermo();


    // Member Functions

        //- Update properties
        virtual void correct();


        // Mass fractions

            //- Return the residual fraction of fuel in the burnt mixture
            virtual tmp<volScalarField> fres() const;

            //- Reset the mixture to an unburnt state and update EGR
            virtual void reset();


        // Thermodynamic state

            //- Unburnt gas temperature [K]
            virtual const volScalarField& Tu() const
            {
                return Tu_;
            }

            //- Unburnt gas enthalpy [J/kg]
            virtual const volScalarField& heu() const
            {
                return heu_;
            }

            //- Unburnt gas enthalpy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& heu()
            {
                return heu_;
            }


        // Derived thermodynamic properties

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

            //- Standard enthalpy of reaction [J/kg]
            virtual tmp<volScalarField> hr() const;

            //- Unburnt gas compressibility [s^2/m^2]
            virtual tmp<volScalarField> psiu() const;

            //- Burnt gas compressibility [s^2/m^2]
            virtual tmp<volScalarField> psib() const;


        // Derived transport properties

            //- Dynamic viscosity of unburnt gas [kg/m/s]
            virtual tmp<volScalarField> muu() const;

            //- Dynamic viscosity of burnt gas [kg/m/s]
            virtual tmp<volScalarField> mub() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const PsiuMulticomponentThermo<BaseThermo>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
    #include "PsiuMulticomponentThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

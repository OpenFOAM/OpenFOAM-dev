/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::UBRhoFluidMulticomponentThermo

Description
    Thermo implementation based on compressibility with additional unburnt
    thermodynamic state

SourceFiles
    UBRhoFluidMulticomponentThermo.C

\*---------------------------------------------------------------------------*/

#ifndef UBRhoFluidMulticomponentThermo_H
#define UBRhoFluidMulticomponentThermo_H

#include "RhoFluidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class UBRhoFluidMulticomponentThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BaseThermo>
class UBRhoFluidMulticomponentThermo
:
    public RhoFluidThermo<BaseThermo>
{

public:

    // Constructors

        //- Construct from mesh and phase name
        UBRhoFluidMulticomponentThermo(const fvMesh&, const word& phaseName);

        //- Disallow default bitwise copy construction
        UBRhoFluidMulticomponentThermo
        (
            const UBRhoFluidMulticomponentThermo<BaseThermo>&
        ) = delete;


    //- Destructor
    virtual ~UBRhoFluidMulticomponentThermo();


    // Member Functions

        // Mass fractions

            //- Return the fuel-oxidant equivalence ratio
            virtual tmp<volScalarField> Phi() const;

            //- Reset the mixture to an unburnt state and update EGR
            virtual void reset
            (
                volScalarField& b,
                volScalarField& c,
                const PtrList<volScalarField>& Yb,
                const volScalarField& heb
            );

        //- Standard enthalpy of formation [J/kg]
        virtual tmp<volScalarField> hf() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=
        (
            const UBRhoFluidMulticomponentThermo<BaseThermo>&
        ) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
    #include "UBRhoFluidMulticomponentThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

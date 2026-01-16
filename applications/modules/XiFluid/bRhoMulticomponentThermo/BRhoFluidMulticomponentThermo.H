/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025-2026 OpenFOAM Foundation
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
    Foam::BRhoFluidMulticomponentThermo

Description
    Thermo implementation based on compressibility with additional unburnt
    thermodynamic state

SourceFiles
    BRhoFluidMulticomponentThermo.C

\*---------------------------------------------------------------------------*/

#ifndef BRhoFluidMulticomponentThermo_H
#define BRhoFluidMulticomponentThermo_H

#include "RhoFluidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class BRhoFluidMulticomponentThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BaseThermo>
class BRhoFluidMulticomponentThermo
:
    public RhoFluidThermo<BaseThermo>
{

public:

    // Constructors

        //- Construct from mesh and phase name
        BRhoFluidMulticomponentThermo(const fvMesh&, const word& phaseName);

        //- Disallow default bitwise copy construction
        BRhoFluidMulticomponentThermo
        (
            const BRhoFluidMulticomponentThermo<BaseThermo>&
        ) = delete;


    //- Destructor
    virtual ~BRhoFluidMulticomponentThermo();


    // Member Functions

        //- Standard enthalpy of formation [J/kg]
        virtual tmp<volScalarField> hf() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=
        (
            const BRhoFluidMulticomponentThermo<BaseThermo>&
        ) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
    #include "BRhoFluidMulticomponentThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

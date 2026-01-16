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
    Foam::URhoFluidMulticomponentThermo

Description
    Thermo implementation based on compressibility with additional unburnt
    thermodynamic state

SourceFiles
    URhoFluidMulticomponentThermo.C

\*---------------------------------------------------------------------------*/

#ifndef URhoFluidMulticomponentThermo_H
#define URhoFluidMulticomponentThermo_H

#include "RhoFluidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class URhoFluidMulticomponentThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BaseThermo>
class URhoFluidMulticomponentThermo
:
    public RhoFluidThermo<BaseThermo>
{

public:

    // Constructors

        //- Construct from mesh and phase name
        URhoFluidMulticomponentThermo(const fvMesh&, const word& phaseName);

        //- Disallow default bitwise copy construction
        URhoFluidMulticomponentThermo
        (
            const URhoFluidMulticomponentThermo<BaseThermo>&
        ) = delete;


    //- Destructor
    virtual ~URhoFluidMulticomponentThermo();


    // Member Functions

        //- Return the fuel-oxidant equivalence ratio
        virtual tmp<volScalarField> Phi() const;

        //- Return the prompt burnt gas composition
        virtual PtrList<volScalarField::Internal> prompt() const;

        //- Reset the mixture to an unburnt state and update EGR
        virtual void reset
        (
            volScalarField& b,
            volScalarField& c,
            const PtrList<volScalarField>& Yb,
            const volScalarField& heb
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=
        (
            const URhoFluidMulticomponentThermo<BaseThermo>&
        ) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
    #include "URhoFluidMulticomponentThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

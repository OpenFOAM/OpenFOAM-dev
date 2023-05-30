/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::heFluidMulticomponentThermo

Description
    Fluid multi-component thermo implementation

SourceFiles
    heFluidMulticomponentThermo.C

\*---------------------------------------------------------------------------*/

#ifndef heFluidMulticomponentThermo_H
#define heFluidMulticomponentThermo_H

#include "heThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class heFluidMulticomponentThermo Declaration
\*---------------------------------------------------------------------------*/

template<class HeThermo>
class heFluidMulticomponentThermo
:
    public HeThermo
{
public:

    // Constructors

        //- Construct from mesh and phase name
        heFluidMulticomponentThermo(const fvMesh&, const word& phaseName);

        //- Disallow default bitwise copy construction
        heFluidMulticomponentThermo
        (
            const heFluidMulticomponentThermo<HeThermo>&
        ) = delete;


    //- Destructor
    virtual ~heFluidMulticomponentThermo();


    // Member Functions

        // Specie transport properties

            //- Dynamic viscosity [kg/m/s]
            virtual scalar mui
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const;

            //- Dynamic viscosity [kg/m/s]
            virtual tmp<volScalarField> mui
            (
                const label speciei,
                const volScalarField& p,
                const volScalarField& T
            ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const heFluidMulticomponentThermo<HeThermo>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
    #include "heFluidMulticomponentThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

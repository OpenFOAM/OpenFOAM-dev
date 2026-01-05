/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::pureLagrangianThermo

Description
    Base-class for pure (i.e., single-component) Lagrangian thermodynamic
    models

SourceFiles
    pureLagrangianThermoI.H
    pureLagrangianThermo.C

\*---------------------------------------------------------------------------*/

#ifndef pureLagrangianThermo_H
#define pureLagrangianThermo_H

#include "basicLagrangianThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class pureLagrangianThermo Declaration
\*---------------------------------------------------------------------------*/

class pureLagrangianThermo
:
    virtual public basicLagrangianThermo
{
public:

    //- Destructor
    virtual ~pureLagrangianThermo();


    // Member Functions

        //- Get the slicer
        inline nil Yslicer() const;

        //- Get the composition of an element
        inline nil elementComposition(const nil, const label) const;

        //- Get the slicer for an injection
        inline nil Yslicer
        (
            const LagrangianInjection&,
            const LagrangianSubMesh&
        ) const;

        //- Get the composition of an element in an injection
        inline nil injectionElementComposition(const nil, const label) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pureLagrangianThermoI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

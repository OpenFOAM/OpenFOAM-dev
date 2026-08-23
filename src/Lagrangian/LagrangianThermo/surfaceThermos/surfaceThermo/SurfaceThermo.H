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
    Foam::SurfaceThermo

Description
    Lagrangian surface thermodynamic model implementation

SourceFiles
    SurfaceThermo.C

\*---------------------------------------------------------------------------*/

#ifndef SurfaceThermo_H
#define SurfaceThermo_H

#include "basicLagrangianSubThermo.H"
#include "dimensionSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class SurfaceThermoName Declaration
\*---------------------------------------------------------------------------*/

TemplateName(SurfaceThermo);


/*---------------------------------------------------------------------------*\
                        Class SurfaceThermo Declaration
\*---------------------------------------------------------------------------*/

template<class MixtureType, class BasicThermoType>
class SurfaceThermo
:
    public SurfaceThermoName,
    public MixtureType,
    public BasicThermoType
{
protected:

    // Protected Member Functions

        //- Return a LagrangianSubScalarField of the given property
        template<class Mixture, class Method, class ... Args>
        tmp<LagrangianSubScalarField> LagrangianSubScalarFieldProperty
        (
            const LagrangianSubMesh& subMesh,
            const word& psiName,
            const dimensionSet& psiDim,
            Mixture mixture,
            Method psiMethod,
            const Args& ... args
        ) const;


public:

    // Typedefs

        //- Mixture type
        typedef MixtureType mixtureType;

        //- Basic thermo
        typedef BasicThermoType basicThermoType;


    //- Disambiguate debug switch used by derivations
    using SurfaceThermoName::debug;


    // Constructors

        //- Construct from a reference to the far thermo
        SurfaceThermo(const basicLagrangianSubThermo& farThermo);

        //- Disallow default bitwise copy construction
        SurfaceThermo(const SurfaceThermo&) = delete;


    //- Destructor
    virtual ~SurfaceThermo();


    // Member Functions

        //- Return the name of the mixture
        virtual word mixtureName() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SurfaceThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

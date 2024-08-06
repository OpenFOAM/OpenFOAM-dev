/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2024 OpenFOAM Foundation
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
    Foam::solidThermophysicalTransportModels::isotropic

Description
    Solid thermophysical transport model for isotropic thermal conductivity

    This is the default transport model for solids and selected automatically if
    the thermophysicalTransport dictionary is not present in the constant or
    region directory.

SourceFiles
    isotropic.C

\*---------------------------------------------------------------------------*/

#ifndef isotropic_H
#define isotropic_H

#include "solidThermophysicalTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidThermophysicalTransportModels
{

/*---------------------------------------------------------------------------*\
                           Class isotropic Declaration
\*---------------------------------------------------------------------------*/

template<class SolidThermophysicalTransportModel>
class isotropic
:
    public SolidThermophysicalTransportModel
{

public:

    typedef typename SolidThermophysicalTransportModel::alphaField
        alphaField;

    //- Runtime type information
    TypeName("isotropic");


    // Constructors

        //- Construct from solid thermophysical properties
        isotropic
        (
            const alphaField& alpha,
            const solidThermo& thermo
        );


    //- Destructor
    virtual ~isotropic()
    {}


    // Member Functions

        //- Read thermophysicalTransport dictionary
        virtual bool read();

        //- Return the heat flux [W/m^2]
        virtual tmp<surfaceScalarField> q() const;

        //- Return null patch heat flux correction [W/m^2]
        virtual tmp<scalarField> qCorr(const label patchi) const;

        //- Return the source term for the energy equation
        virtual tmp<fvScalarMatrix> divq(volScalarField& he) const;

        //- Correct the isotropic viscosity
        virtual void predict();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidThermophysicalTransportModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "isotropic.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

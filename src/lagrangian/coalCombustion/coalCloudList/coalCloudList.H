/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#ifndef coalCloudList_H
#define coalCloudList_H

#include "coalCloud.H"
#include "volFieldsFwd.H"
#include "fvMatricesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class coalCloudList Declaration
\*---------------------------------------------------------------------------*/

class coalCloudList
:
    public PtrList<coalCloud>
{
private:

    //- Reference to the mesh
    const fvMesh& mesh_;


public:

    // Constructor

        coalCloudList
        (
            const volScalarField& rho,
            const volVectorField& U,
            const dimensionedVector& g,
            const SLGThermo& slgThermo
        );


    // Member Functions

        // Evolution

            //- Evolve the cloud collection
            void evolve();


        // Source terms

            //- Return const reference to momentum source
            inline tmp<DimensionedField<vector, volMesh> > UTrans() const;

            //- Return tmp momentum source term
            inline tmp<fvVectorMatrix> SU(volVectorField& U) const;

            //- Sensible enthalpy transfer [J/kg]
            inline tmp<DimensionedField<scalar, volMesh> > hsTrans() const;

            //- Return sensible enthalpy source term [J/kg/m3/s]
            inline tmp<fvScalarMatrix> Sh(volScalarField& hs) const;


            //- Return mass source term for specie i - specie eqn
            inline tmp<fvScalarMatrix> SYi
            (
                const label i,
                volScalarField& Yi
            ) const;

            //- Return total mass transfer [kg/m3]
            inline tmp<DimensionedField<scalar, volMesh> > rhoTrans() const;

            //- Return tmp total mass source for carrier phase
            //  - fully explicit
            inline tmp<DimensionedField<scalar, volMesh> > Srho() const;

            //- Return tmp total mass source for carrier phase specie i
            //  - fully explicit
            inline tmp<DimensionedField<scalar, volMesh> > Srho
            (
                const label i
            ) const;

            //- Return total mass source term [kg/m3/s]
            inline tmp<fvScalarMatrix> Srho(volScalarField& rho) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "coalCloudListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

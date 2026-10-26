/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2026 OpenFOAM Foundation
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
    Foam::chemistryModels::Standard

Description
    Extension to Foam::chemistryModels::standard templated on thermo
    and provides stiff ODE integration functions.

    Integrates a standard OpenFOAM reaction system and evaluation of chemical
    source terms with optional support for mechanism reduction and tabulation.

See also
    Foam::chemistryModels::standard

SourceFiles
    Standard_chemistryModelI.H
    standard_chemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef Standard_chemistryModel_H
#define Standard_chemistryModel_H

#include "standard_chemistryModel.H"
#include "chemistryReductionMethod.H"
#include "chemistryTabulationMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace chemistryModels
{

/*---------------------------------------------------------------------------*\
                           Class Standard Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class Standard
:
    public standard
{
    // Private data

        //- Type of the Jacobian to be calculated
        const jacobianType jacobianType_;

        //- Thermodynamic data of the species
        const PtrList<ThermoType>& specieThermos_;


public:

    //- Runtime type information
    TypeName("standard");


    // Constructors

        //- Construct from thermo
        Standard(const fluidMulticomponentThermo& thermo);

        //- Disallow default bitwise copy construction
        Standard(const Standard&) = delete;


    //- Destructor
    virtual ~Standard();


    // Member Functions

        // ODE solution functions

            //- Calculate the ODE derivatives
            virtual void derivatives
            (
                const scalar t,
                const scalarField& YT,
                const label li,
                scalarField& dYTdt
            ) const;

            //- Calculate the ODE jacobian
            virtual void jacobian
            (
                const scalar t,
                const scalarField& YT,
                const label li,
                scalarField& dYTdt,
                scalarSquareMatrix& J
            ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const Standard&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace chemistryModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Standard_chemistryModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

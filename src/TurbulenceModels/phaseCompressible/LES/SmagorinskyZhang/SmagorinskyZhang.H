/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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
    Foam::LESModels::SmagorinskyZhang

Group
    grpLESTurbulence

Description
    The Smagorinsky SGS model including bubble-generated turbulence

    Reference:
    \verbatim
        "Numerical simulation of the dynamic flow behavior in a bubble column:
        A study of closures for turbulence and interface forces"
        D. Zhang,
        N.G. Deen,
        J.A.M. Kuipers,
        Chemical Engineering Science 61 (2006) pp 7593-7608.
    \endverbatim

    The default model coefficients correspond to the following:
    \verbatim
        SmagorinskyZhangCoeffs
        {
            Ck              0.094;
            Ce              1.048;
            Cmub            0.6;
        }
    \endverbatim

SourceFiles
    SmagorinskyZhang.C

\*---------------------------------------------------------------------------*/

#ifndef SmagorinskyZhang_H
#define SmagorinskyZhang_H

#include "LESModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class SmagorinskyZhang Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class SmagorinskyZhang
:
    public Smagorinsky<BasicTurbulenceModel>
{
    // Private data

        mutable const PhaseCompressibleTurbulenceModel
        <
            typename BasicTurbulenceModel::transportModel
        > *gasTurbulencePtr_;


    // Private Member Functions

        //- Return the turbulence model for the gas phase
        const PhaseCompressibleTurbulenceModel
        <
            typename BasicTurbulenceModel::transportModel
        >&
        gasTurbulence() const;

        // Disallow default bitwise copy construct and assignment
        SmagorinskyZhang(const SmagorinskyZhang&);
        SmagorinskyZhang& operator=(const SmagorinskyZhang&);


protected:

    // Protected data

        // Model coefficients

            dimensionedScalar Cmub_;


    // Protected Member Functions

        virtual void correctNut();


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("SmagorinskyZhang");


    // Constructors

        //- Construct from components
        SmagorinskyZhang
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~SmagorinskyZhang()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "SmagorinskyZhang.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

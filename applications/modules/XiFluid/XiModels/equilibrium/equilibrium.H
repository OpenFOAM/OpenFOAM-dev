/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::XiModels::equilibrium

Description
    Simple equilibrium model for flame wrinkling \c Xi

    Note that this model does not guarantee that the brush will attain a stable
    thickness in unbounded propagation which would require a \c Xi distribution
    through the flame, lower at the front and higher at the rear, see
    Foam::XiEqModels::linearEquilibrium or Foam::XiEqModels::cubicEquilibrium
    for alternative equilibrium models which provide this.

    See Foam::XiModel for more details on flame wrinkling modelling.

    References:
    \verbatim
        Weller, H. G. (1993).
        The development of a new flame area combustion model
        using conditional averaging.
        Thermo-fluids section report TF 9307.

        Weller, H. G., Tabor, G., Gosman, A. D., & Fureby, C. (1998, January).
        Application of a flame-wrinkling LES combustion model
        to a turbulent mixing layer.
        In Symposium (International) on combustion
        (Vol. 27, No. 1, pp. 899-907). Elsevier.
    \endverbatim

SourceFiles
    equilibrium.C

\*---------------------------------------------------------------------------*/

#ifndef equilibrium_H
#define equilibrium_H

#include "XiModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class XiEqModel;

namespace XiModels
{

/*---------------------------------------------------------------------------*\
                          Class equilibrium Declaration
\*---------------------------------------------------------------------------*/

class equilibrium
:
    public XiModel
{
protected:

        autoPtr<XiEqModel> XiEqModel_;

        //- Update coefficients from given dictionary
        virtual bool readCoeffs(const dictionary& dict);


public:

    //- Runtime type information
    TypeName("equilibrium");


    // Constructors

        //- Construct from components
        equilibrium
        (
            const dictionary& dict,
            const psiuMulticomponentThermo& thermo,
            const fluidThermoThermophysicalTransportModel& turbulence,
            const volScalarField& Su
        );

        //- Disallow default bitwise copy construction
        equilibrium(const equilibrium&) = delete;


    //- Destructor
    virtual ~equilibrium();


    // Member Functions

        //- Return the flame diffusivity
        virtual tmp<volScalarField> Db() const;

        //- Correct the flame-wrinkling Xi
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const equilibrium&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace XiModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

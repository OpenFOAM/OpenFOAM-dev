/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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

Global
    MULES

Description
    MULES: Multidimensional universal limiter for explicit solution.

    Solve a convective-only transport equation using an explicit universal
    multi-dimensional limiter.

    Parameters are the variable to solve, the normal convective flux and the
    actual explicit flux of the variable which is also used to return limited
    flux used in the bounded-solution.

SourceFiles
    MULES.C
    MULESTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef MULES_H
#define MULES_H

#include "Switch.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "primitiveFieldsFwd.H"
#include "geometricOneField.H"
#include "zero.H"
#include "zeroField.H"
#include "UPtrList.H"
#include "HashSet.H"
#include "UniformField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace MULES
{

NamespaceName("MULES");

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// MULES control structure
struct control
{
    //- Optional maximum number of limiter iterations
    //  Defaults to 3
    label nIter;

    //- Optional limiter convergence tolerance
    //  1e-2 is usually sufficient
    //  but 1e-3 can be used for very tight convergence
    //  Defaults to 0, i.e. unused and nIter iterations are performed
    scalar tol;

    //- Optional switch to select global bounds only
    //  rather than local and global bounds
    Switch globalBounds;

    //- Optional coefficient to relax the local boundedness constraint
    //  by adding extremaCoeff*local_range to the allowed range
    //  allowing the solution to approach the global bounds in case temporary
    //  smearing has caused the solution to deviate from the entire range.
    //
    //  A value of 0 enforces strict local boundedness,
    //  a value of 1 enforces the global bounds only.
    //  Defaults to 0
    scalar extremaCoeff;

    //- Optional coefficient to relax the local boundedness constraint
    //  in the cells adjacent to boundaries only
    //  Defaults to 0
    scalar boundaryExtremaCoeff;

    //- Optional coefficient to reduce the allowed range of the solution to
    //  smoothingCoeff*local_range.
    //
    //  This option should only be used to improve the smoothness of the
    //  solution at low Courant number (explicit solution with sub-cycling)
    //  otherwise excessive smearing may result.
    scalar smoothingCoeff;

    //- Null constructor
    //  Should be followed by a call to read(dict)
    control()
    {}

    //- Construct from dict and set the controls
    control(const dictionary& dict);

    //- Read dict and set the controls
    void read(const dictionary& dict);
};

template<class RdeltaTType, class RhoType, class SpType, class SuType>
void explicitSolve
(
    const RdeltaTType& rDeltaT,
    const RhoType& rho,
    volScalarField& psi,
    const surfaceScalarField& psiPhi,
    const SpType& Sp,
    const SuType& Su
);

template<class RhoType>
void explicitSolve
(
    const RhoType& rho,
    volScalarField& psi,
    const surfaceScalarField& psiPhi
);

template<class RhoType, class SpType, class SuType>
void explicitSolve
(
    const RhoType& rho,
    volScalarField& psi,
    const surfaceScalarField& psiPhi,
    const SpType& Sp,
    const SuType& Su
);

template<class RhoType, class PsiMaxType, class PsiMinType>
void explicitSolve
(
    const control& controls,
    const RhoType& rho,
    volScalarField& psi,
    const surfaceScalarField& phiBD,
    surfaceScalarField& psiPhi,
    const PsiMaxType& psiMax,
    const PsiMinType& psiMin
);

template
<
    class RhoType,
    class SpType,
    class SuType,
    class PsiMaxType,
    class PsiMinType
>
void explicitSolve
(
    const control& controls,
    const RhoType& rho,
    volScalarField& psi,
    const surfaceScalarField& phiBD,
    surfaceScalarField& psiPhi,
    const SpType& Sp,
    const SuType& Su,
    const PsiMaxType& psiMax,
    const PsiMinType& psiMin
);

template
<
    class RdeltaTType,
    class RhoType,
    class SpType,
    class PsiMaxType,
    class PsiMinType
>
void limiter
(
    const control& controls,
    surfaceScalarField& lambda,
    const RdeltaTType& rDeltaT,
    const RhoType& rho,
    const volScalarField& psi,
    const scalarField& SuCorr,
    const surfaceScalarField& phi,
    const surfaceScalarField& phiCorr,
    const SpType& Sp,
    const PsiMaxType& psiMax,
    const PsiMinType& psiMin
);

template
<
    class RdeltaTType,
    class RhoType,
    class SpType,
    class SuType,
    class PsiMaxType,
    class PsiMinType
>
void limit
(
    const control& controls,
    const RdeltaTType& rDeltaT,
    const RhoType& rho,
    const volScalarField& psi,
    const surfaceScalarField& phi,
    surfaceScalarField& psiPhi,
    const SpType& Sp,
    const SuType& Su,
    const PsiMaxType& psiMax,
    const PsiMinType& psiMin,
    const bool returnCorr
);

template
<
    class RdeltaTType,
    class RhoType,
    class SpType,
    class SuType,
    class PsiMaxType,
    class PsiMinType
>
void limit
(
    const control& controls,
    const RdeltaTType& rDeltaT,
    const RhoType& rho,
    const volScalarField& psi,
    const surfaceScalarField& phi,
    const surfaceScalarField& phiBD,
    surfaceScalarField& psiPhi,
    const SpType& Sp,
    const SuType& Su,
    const PsiMaxType& psiMax,
    const PsiMinType& psiMin,
    const bool returnCorr
);

template
<
    class RhoType,
    class SpType,
    class SuType,
    class PsiMaxType,
    class PsiMinType
>
void limit
(
    const control& controls,
    const RhoType& rho,
    const volScalarField& psi,
    const surfaceScalarField& phi,
    surfaceScalarField& psiPhi,
    const SpType& Sp,
    const SuType& Su,
    const PsiMaxType& psiMax,
    const PsiMinType& psiMin,
    const bool returnCorr
);

template
<
    class RhoType,
    class SpType,
    class SuType,
    class PsiMaxType,
    class PsiMinType
>
void limit
(
    const control& controls,
    const RhoType& rho,
    const volScalarField& psi,
    const surfaceScalarField& phi,
    const surfaceScalarField& phiBD,
    surfaceScalarField& psiPhi,
    const SpType& Sp,
    const SuType& Su,
    const PsiMaxType& psiMax,
    const PsiMinType& psiMin,
    const bool returnCorr
);

void limitSum(UPtrList<scalarField>& psiPhiCorrs);

void limitSum
(
    UPtrList<scalarField>& phiPsiCorrs,
    const UPtrList<scalarField>& phiPsiFixedCorrs
);

template
<
    template<class> class ConstraintList,
    class ConstraintField,
    template<class> class PsiPhiList
>
void limitSum
(
    const ConstraintList<const ConstraintField>& constraints,
    PsiPhiList<surfaceScalarField>& psiPhiCorrs,
    const surfaceScalarField& phi
);

template
<
    template<class> class ConstraintList,
    class ConstraintField,
    template<class> class PsiPhiList
>
void limitSum
(
    const PtrList<surfaceScalarField>& alphaPhiBDs,
    const ConstraintList<const ConstraintField>& constraints,
    PsiPhiList<surfaceScalarField>& psiPhis,
    const surfaceScalarField& phi
);

template
<
    template<class> class AlphaList,
    template<class> class ConstraintList,
    class ConstraintField,
    template<class> class PsiPhiList
>
void limitSum
(
    const AlphaList<const volScalarField>& alphas,
    const ConstraintList<const ConstraintField>& constraints,
    PsiPhiList<surfaceScalarField>& psiPhis,
    const surfaceScalarField& phi
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace MULES
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "MULESlimiter.C"
    #include "MULESTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

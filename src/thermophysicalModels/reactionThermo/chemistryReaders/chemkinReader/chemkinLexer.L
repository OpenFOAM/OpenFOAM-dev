/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

%{

#undef yyFlexLexer

#include "chemkinReader.H"

#include "error.H"
#include "IStringStream.H"
// For EOF only
#include <cstdio>

// flex input buffer size
int Foam::chemkinReader::yyBufSize = YY_BUF_SIZE;

// Dummy yyFlexLexer::yylex() to keep the linker happy. It is not called
//! \cond dummy
int yyFlexLexer::yylex()
{
    FatalErrorInFunction
        << "should not have called this function"
        << abort(Foam::FatalError);

    return 0;
}
//! \endcond


// Dummy yywrap to keep yylex happy at compile time.
// It is called by yylex but is not used as the mechanism to change file.
// See <<EOF>>
//! \cond dummy
#if YY_FLEX_MINOR_VERSION < 6 && YY_FLEX_SUBMINOR_VERSION < 34
extern "C" int yywrap()
#else
int yyFlexLexer::yywrap()
#endif
{
    return 1;
}
//! \endcond


Foam::string foamSpecieString(const char* YYText)
{
    Foam::string specieString(YYText);
    return specieString;
}


Foam::word foamName(const char* YYText)
{
    Foam::string fn(YYText);
    Foam::string::stripInvalid<Foam::word>(fn);
    return fn;
}


Foam::word foamName(const Foam::string& s)
{
    Foam::string fn(s);
    Foam::string::stripInvalid<Foam::word>(fn);
    return fn;
}


 /* ------------------------------------------------------------------------- *\
   ------ cppLexer::yylex()
 \* ------------------------------------------------------------------------- */

#define YY_DECL int Foam::chemkinReader::lex()

%}

one_space             [ \t\f\r]
space                 {one_space}*
some_space            {one_space}+

alpha                 [_A-Za-z]
digit                 [0-9]

exponent_part         [eEdD][-+]?{digit}+
fractional_constant   [-+]?(({digit}*"."{digit}+)|({digit}+"."?))

floatNum              (({fractional_constant}{exponent_part}?)|({digit}+{exponent_part}))


elements              {space}("ELEMENTS"|"ELEM"){space}
species               {space}("SPECIES"|"SPECIE"|"SPEC"){space}
thermoAll             {space}"THERMO"{some_space}"ALL"{space}
thermo                {space}"THERMO"{space}
reactions             {space}("REACTIONS"|"REAC"){space}
end                   {space}"END"{space}

elementName           {space}([A-Z]|([A-Z][A-Z])){space}
startIsotopeMolW      {space}"/"{space}
isotopeMolW           {space}{floatNum}{space}"/"{space}

specieName            {space}[A-Za-z](([A-Za-z0-9)*+-])|("("[^+]))*{space}
Y                {space}{floatNum}{space}

thermoTemp            .{10}

thermoSpecieName      .{18}
thermoDate            .{6}
thermoFormula         .{20}
thermoPhase           S|L|G|C
thermoLowTemp         .{10}
thermoHighTemp        .{10}
thermoCommonTemp      .{8}
thermoFormula2        .{5}
thermoCoeff           .{15}
thermoLineLabel1      " "1{space}
thermoLineLabel2      " "{4}2{space}
thermoLineLabel3      " "{4}3{space}
thermoLineLabel4      " "{4}4{space}

specieDelimiter       {space}"+"{space}
reversibleReactionDelimiter {space}("="|"<=>"){space}
irreversibleReactionDelimiter {space}"=>"{space}
startPDependentSpecie {space}"("{space}"+"{space}
pDependentSpecie      {specieName}")"{space}
reactionCoeffs        {space}{floatNum}{some_space}{floatNum}{some_space}{floatNum}{space}
reactionKeyword       {space}[A-Za-z](([A-Za-z0-9)*-])|("("[^+]))*{space}
reactionKeywordSlash  {reactionKeyword}"/"{space}
thirdBodyEfficiency   {space}{floatNum}{space}"/"{space}
startReactionCoeffs   {space}"/"{space}
endReactionCoeffs     {space}"/"{space}
reactionCoeff         {space}{floatNum}{space}

calPerMol             {space}"CAL/MOLE"{space}
kcalPerMol            {space}"KCAL/MOLE"{space}
joulePerMol           {space}"JOULES/MOLE"{space}
kelvins               {space}"KELVINS"{space}
otherReactionsUnit    {space}("EVOLTS"|"MOLES"|"MOLECULES"){space}

cal                   {space}"CAL"{space}
kcal                  {space}"KCAL"{space}
joule                 {space}"JOUL"{space}
kjoule                {space}"KJOU"{space}
kelvin                {space}("KELV"|"KELVIN"){space}
otherReactionUnit     {space}("MOLE"|"MOLECULE"|"EVOL"|"EVOLTS"){space}


 /* ------------------------------------------------------------------------- *\
                      -----  Exclusive start states -----
 \* ------------------------------------------------------------------------- */

%option stack

%x readElements
%x readIsotopeMolW
%x readSpecies
%x readThermoAll
%x readThermoSpecieName
%x readThermoDate
%x readThermoFormula
%x readThermoPhase
%x readThermoTemps
%x readThermoFormula2
%x readThermoLineLabel1
%x readThermoCoeff1
%x readThermoLineLabel2
%x readThermoCoeff2
%x readThermoLineLabel3
%x readThermoCoeff3
%x readThermoLineLabel4
%x readReactionsUnits
%x readReactionKeyword
%x readSpecieNamePlus
%x readReactionDelimiter
%x readPDependentSpecie
%x readThirdBodyEfficiency
%x readReactionCoeffs
%x readTdepSpecie
%x readReactionOrderSpecie
%x readReactionOrder
%x readReactionUnit
%x CHEMKINError

%%

%{

static const char* stateNames[30] =
{
    "reading CHEMKIN III file",
    "reading elements",
    "reading isotope molecular weight",
    "reading species",
    "reading all thermodynamic data temperatures",
    "reading thermodynamic specie name",
    "reading thermodynamic data date",
    "reading thermodynamic data specie formula",
    "reading thermodynamic data specie phase",
    "reading thermodynamic data temperatures",
    "reading thermodynamic data specie formula (supplement)",
    "reading thermodynamic data line label 1",
    "reading thermodynamic data coefficient set 1",
    "reading thermodynamic data line label 2",
    "reading thermodynamic data coefficient set 2",
    "reading thermodynamic data line label 3",
    "reading thermodynamic data coefficient set 3",
    "reading thermodynamic data line label 4",
    "reading reaction units",
    "reading reaction specie/keyword",
    "reading reaction specie",
    "reading reaction delimiter",
    "reading reaction pressure dependent specie",
    "reading third-body efficiency",
    "reading reaction coeff",
    "reading temperature dependent specie name",
    "reading reaction order specie name",
    "reading reaction order",
    "reading reaction unit",
    "error"
};

static const char* stateExpects[30] =
{
    "'ELEMENTS' or 'ELEM', 'SPECIES' or 'SPEC', 'THERMO' or 'THERMO ALL', 'REACTIONS'",
    "<elementName>, <isotopeName> / or 'END'",
    "<scalar>/",
    "<specieName> or 'END'",
    "<scalar><scalar><scalar> (3F10.0)",
    "<word> (18A1)",
    "<word> (6A1)",
    "<word><label><word><label><word><label><word><label> (4(2A1,I3))",
    "<char> (A1)",
    "<scalar><scalar><scalar> (E10.0E10.0E8.0)",
    "<word><label> (2A1,I3)",
    "'1' (I1)",
    "<scalar><scalar><scalar><scalar><scalar> (5(E15.0))",
    "'2' (I1)",
    "<scalar><scalar><scalar><scalar><scalar> (5(E15.0))",
    "'3' (I1)",
    "<scalar><scalar><scalar><scalar> (4(E15.0))",
    "'4' (I1)",
    "'CAL/MOLE', 'KCAL/MOLE', 'JOULES/MOLE', 'KELVINS', 'EVOLTS', 'MOLES' or 'MOLECULES'",
    "<word> or <label>",
    "'+'",
    "'+', '=', '<=>', '=>', '(+<word>', '<scalar> <scalar> <scalar>'",
    "<word>')'",
    "<scalar>'/'",
    "<scalar>",
    "<word>",
    "<word>",
    "<scalar>",
    "'MOLE', 'MOLECULE', 'CAL', 'KCAL', 'JOUL', 'KJOU', 'KELV', 'KELVIN', 'EVOL' or 'EVOLTS'",
    ""
};

string startError;

static const scalar RRjoule = 8.31451; // J/kg-mol-K
static const scalar RRcal = 1.987316;  // cal/g-mol-K

scalar RRreactions = RRcal;
scalar RRreaction = RRcal;

scalar allCommonT = 1000.0;

word currentElementName;
label currentElementIndex = 0;

word currentSpecieName;
label currentSpecieIndex = 0;
label nSpecieElements = 0;
List<specieElement> currentSpecieComposition(5);

scalar currentLowT = 0;
scalar currentHighT = 0;
scalar currentCommonT = 0;
gasHThermoPhysics::coeffArray highCpCoeffs(scalarList(7));
gasHThermoPhysics::coeffArray lowCpCoeffs(scalarList(7));

gasHReaction::specieCoeffs currentSpecieCoeff;

DynamicList<gasHReaction::specieCoeffs> lhs;
DynamicList<gasHReaction::specieCoeffs> rhs;

scalarList ArrheniusCoeffs(3);
DynamicList<scalar> reactionCoeffs;
scalarList thirdBodyEfficiencies;
label currentThirdBodyIndex = -1;

word reactionCoeffsName = word::null;
HashTable<scalarList> reactionCoeffsTable;

DynamicList<gasHReaction::specieCoeffs> *lrhsPtr = &lhs;

reactionType rType = unknownReactionType;
reactionRateType rrType = Arrhenius;
fallOffFunctionType fofType = unknownFallOffFunctionType;
word pDependentSpecieName = word::null;
label lhsThirdBodyCounter = 0;
label rhsThirdBodyCounter = 0;

bool finishReaction = false;

%}

 /* ------------------------------------------------------------------------- *\
                            ------ Start Lexing ------
 \* ------------------------------------------------------------------------- */

 /* ------------------------------------------------------------------------- *\
    ------ Discard comments being careful to count comments
 \* ------------------------------------------------------------------------- */

<*>{space}"!".* { // Remove one line comments
    }

 /* ------------------------------------------------------------------------- *\
    ------ Read elements
 \* ------------------------------------------------------------------------- */

{elements} {
        BEGIN(readElements);
    }

<readElements>{elementName} {
        currentElementName = foamName(YYText());
        correctElementName(currentElementName);

        if (!elementIndices_.found(currentElementName))
        {
            elementIndices_.insert(currentElementName, currentElementIndex++);
            elementNames_.append(currentElementName);
        }
        else
        {
            WarningInFunction
                << "element " << currentElementName
                << " already in table." << endl;
        }
    }

<readElements>{startIsotopeMolW} {
        BEGIN(readIsotopeMolW);
    }

<readIsotopeMolW>{isotopeMolW} {
        isotopeAtomicWts_.insert(currentElementName, stringToScalar(YYText()));
        BEGIN(readElements);
    }

<readElements>{end} {
        BEGIN(INITIAL);
    }

 /* ------------------------------------------------------------------------- *\
    ------ Read species
 \* ------------------------------------------------------------------------- */

<INITIAL,readElements>{species} {
        BEGIN(readSpecies);
    }

<readSpecies>{specieName} {
        word specieName(foamName(foamSpecieString(YYText())));

        if (specieName == "THERMO")
        {
            specieNames_.shrink();
            speciesTable_ = specieNames_;
            thirdBodyEfficiencies.setSize(specieNames_.size());
            thirdBodyEfficiencies = 1.0;
            BEGIN(readThermoSpecieName);
        }
        else if (specieName == "END")
        {
            specieNames_.shrink();
            speciesTable_ = specieNames_;
            thirdBodyEfficiencies.setSize(specieNames_.size());
            thirdBodyEfficiencies = 1.0;
            BEGIN(INITIAL);
        }
        else
        {
            if (!specieIndices_.found(specieName))
            {
                specieNames_.append(specieName);
                specieIndices_.insert(specieName, currentSpecieIndex++);
            }
            else
            {
                WarningInFunction
                    << "specie " << specieName
                    << " already in table." << endl;
            }
        }
    }

 /* ------------------------------------------------------------------------- *\
    ------ Read thermo
 \* ------------------------------------------------------------------------- */

<INITIAL,readSpecies>{thermo} {
        BEGIN(readThermoSpecieName);
    }

<INITIAL,readSpecies>{thermoAll} {
        BEGIN(readThermoAll);
    }

<readThermoAll>{thermoTemp}{thermoTemp}{thermoTemp} {
        string temperaturesString(YYText());
        // scalar lowestT(stringToScalar(temperaturesString(0, 10)));
        allCommonT = stringToScalar(temperaturesString(10, 10));
        // scalar highestT(stringToScalar(temperaturesString(20, 10)));
        BEGIN(readThermoSpecieName);
    }

<readThermoSpecieName>{thermoSpecieName} {
        string specieString(foamSpecieString(YYText()));
        if (newFormat_)
        {
            specieString.replaceAll(" ", "_");
            size_t strEnd = specieString.find_last_not_of('_');
            currentSpecieName = specieString.substr(0, strEnd + 1);
        }
        else
        {
            size_t spacePos = specieString.find(' ');
            if (spacePos != string::npos)
            {
                currentSpecieName = specieString(0, spacePos);
            }
            else
            {
                currentSpecieName = specieString;
            }
        }

        BEGIN(readThermoDate);
    }

<readThermoDate>{thermoDate} {
        // string thermoDate(YYText());
        // Date is not currently used
        BEGIN(readThermoFormula);
    }

<readThermoFormula>{thermoFormula} {
        string thermoFormula(YYText());

        nSpecieElements = 0;
        currentSpecieComposition.setSize(5);

        for (int i=0; i<4; i++)
        {
            word elementName(foamName(thermoFormula(5*i, 2)));
            label nAtoms = atoi(thermoFormula(5*i + 2, 3).c_str());

            if (elementName.size() && nAtoms)
            {
                correctElementName(elementName);
                currentSpecieComposition[nSpecieElements].name() =
                    elementName;
                currentSpecieComposition[nSpecieElements++].nAtoms() = nAtoms;
            }
        }

        BEGIN(readThermoPhase);
    }

<readThermoPhase>{thermoPhase} {
        char phaseChar = YYText()[0];

        switch (phaseChar)
        {
            case 'S':
                speciePhase_.insert(currentSpecieName, solid);
            break;

            case 'L':
                speciePhase_.insert(currentSpecieName, liquid);
            break;

            case 'G':
                speciePhase_.insert(currentSpecieName, gas);
            break;
        }

        BEGIN(readThermoTemps);
    }

<readThermoTemps>{thermoLowTemp}{thermoHighTemp}{thermoCommonTemp} {
        string temperaturesString(YYText());
        currentLowT = stringToScalar(temperaturesString(0, 10));
        currentHighT = stringToScalar(temperaturesString(10, 10));
        currentCommonT = stringToScalar(temperaturesString(20, 8));

        if (currentCommonT < SMALL)
        {
            currentCommonT = allCommonT;
        }

        BEGIN(readThermoFormula2);
    }

<readThermoFormula2>{thermoFormula2} {
        string thermoFormula(YYText());
        word elementName(foamName(thermoFormula(0, 2)));
        label nAtoms = atoi(thermoFormula(2, 3).c_str());

        if
        (
            elementName.size()
         && elementName.find('0') == string::npos
         && nAtoms
        )
        {
            correctElementName(elementName);
            currentSpecieComposition[nSpecieElements].name() =
                elementName;
            currentSpecieComposition[nSpecieElements++].nAtoms() = nAtoms;
        }

        currentSpecieComposition.setSize(nSpecieElements);

        speciesCompositionTable::iterator specieCompositionIter
        (
            speciesComposition_.find(currentSpecieName)
        );

        if (specieCompositionIter != speciesComposition_.end())
        {
            speciesComposition_.erase(specieCompositionIter);
        }

        speciesComposition_.insert
        (
            currentSpecieName,
            currentSpecieComposition
        );

        BEGIN(readThermoLineLabel1);
    }

<readThermoLineLabel1>{thermoLineLabel1} {
        BEGIN(readThermoCoeff1);
    }

<readThermoCoeff1>{thermoCoeff}{thermoCoeff}{thermoCoeff}{thermoCoeff}{thermoCoeff} {
        string thermoCoeffString(YYText());

        highCpCoeffs[0] = stringToScalar(thermoCoeffString(0, 15));
        highCpCoeffs[1] = stringToScalar(thermoCoeffString(15, 15));
        highCpCoeffs[2] = stringToScalar(thermoCoeffString(30, 15));
        highCpCoeffs[3] = stringToScalar(thermoCoeffString(45, 15));
        highCpCoeffs[4] = stringToScalar(thermoCoeffString(60, 15));

        BEGIN(readThermoLineLabel2);
    }

<readThermoLineLabel2>{thermoLineLabel2} {
        BEGIN(readThermoCoeff2);
    }

<readThermoCoeff2>{thermoCoeff}{thermoCoeff}{thermoCoeff}{thermoCoeff}{thermoCoeff} {
        string thermoCoeffString(YYText());

        highCpCoeffs[5] = stringToScalar(thermoCoeffString(0, 15));
        highCpCoeffs[6] = stringToScalar(thermoCoeffString(15, 15));

        lowCpCoeffs[0] = stringToScalar(thermoCoeffString(30, 15));
        lowCpCoeffs[1] = stringToScalar(thermoCoeffString(45, 15));
        lowCpCoeffs[2] = stringToScalar(thermoCoeffString(60, 15));

        BEGIN(readThermoLineLabel3);
    }

<readThermoLineLabel3>{thermoLineLabel3} {
        BEGIN(readThermoCoeff3);
    }

<readThermoCoeff3>{thermoCoeff}{thermoCoeff}{thermoCoeff}{thermoCoeff}{thermoCoeff} {
        string thermoCoeffString(YYText());

        lowCpCoeffs[3] = stringToScalar(thermoCoeffString(0, 15));
        lowCpCoeffs[4] = stringToScalar(thermoCoeffString(15, 15));
        lowCpCoeffs[5] = stringToScalar(thermoCoeffString(30, 15));
        lowCpCoeffs[6] = stringToScalar(thermoCoeffString(45, 15));

        BEGIN(readThermoLineLabel4);
    }

<readThermoLineLabel4>{thermoLineLabel4} {

        HashPtrTable<gasHThermoPhysics>::iterator specieThermoIter
        (
            speciesThermo_.find(currentSpecieName)
        );

        if (specieThermoIter != speciesThermo_.end())
        {
            speciesThermo_.erase(specieThermoIter);
        }

        speciesThermo_.insert
        (
            currentSpecieName,
            new gasHThermoPhysics
            (
                janafThermo<perfectGas<specie>>
                (
                    specie
                    (
                        currentSpecieName,
                        1.0,
                        molecularWeight(currentSpecieComposition)
                    ),
                    currentLowT,
                    currentHighT,
                    currentCommonT,
                    highCpCoeffs,
                    lowCpCoeffs,
                    true
                ),
                transportDict_.subDict(currentSpecieName)
            )
        );

        BEGIN(readThermoSpecieName);
    }

<readThermoSpecieName>{end} {
        Reaction<gasHThermoPhysics>::TlowDefault = max
        (
            Reaction<gasHThermoPhysics>::TlowDefault,
            currentLowT
        );

        Reaction<gasHThermoPhysics>::ThighDefault = min
        (
            Reaction<gasHThermoPhysics>::ThighDefault,
            currentHighT
        );

        BEGIN(INITIAL);
    }

 /* ------------------------------------------------------------------------- *\
    ------ Read reactions
 \* ------------------------------------------------------------------------- */

<INITIAL,readThermoSpecieName>{reactions} {
        currentSpecieCoeff.stoichCoeff = 1.0;
        currentSpecieCoeff.exponent = currentSpecieCoeff.stoichCoeff;
        BEGIN(readReactionsUnits);
    }

<readReactionsUnits>{calPerMol} {
        RRreactions = RRcal;
    }

<readReactionsUnits>{kcalPerMol} {
        RRreactions = RRcal/1000.0;
    }

<readReactionsUnits>{joulePerMol} {
        RRreactions = RRjoule;
    }

<readReactionsUnits>{kelvins} {
        RRreactions = 1.0;
    }

<readReactionsUnits>{otherReactionsUnit} {
    }

<readReactionsUnits>\n {
        lineNo_++;
        BEGIN(readReactionKeyword);
    }

<readReactionKeyword>{Y} {
        currentSpecieCoeff.stoichCoeff = atof(YYText());
        currentSpecieCoeff.exponent = currentSpecieCoeff.stoichCoeff;
    }

<readReactionKeyword>{reactionKeyword} {

        word keyword(foamName(YYText()));

        HashTable<int>::iterator reactionKeywordIter
        (
            reactionKeywordTable_.find(keyword)
        );

        if (reactionKeywordIter != reactionKeywordTable_.end())
        {
            switch(reactionKeywordIter())
            {
                case duplicateReactionType:
                {
                    BEGIN(readReactionKeyword);
                    break;
                }

                case thirdBodyReactionType:
                {
                    if (rrType != Arrhenius && rrType != thirdBodyArrhenius)
                    {
                        FatalErrorInFunction
                            << "Attempt to set reaction rate type to"
                               " thirdBodyArrhenius when it is already set to "
                            << reactionRateTypeNames[rrType]
                            << " on line " << lineNo_
                            << exit(FatalError);
                    }

                    if (pDependentSpecieName.size())
                    {
                        FatalErrorInFunction
                            << "A non-pressure dependent third-body appears in"
                               " the pressure dependent reaction on line "
                            << lineNo_
                            << exit(FatalError);
                    }

                    rrType = thirdBodyArrhenius;

                    if (lrhsPtr == &lhs)
                    {
                        lhsThirdBodyCounter++;
                    }
                    else
                    {
                        rhsThirdBodyCounter++;
                    }

                    BEGIN(readReactionDelimiter);
                    break;
                }

                case plasmaMomentumTransfer:
                {
                    FatalErrorInFunction
                        << "Plasma momentum-transfer in reaction on line "
                        << lineNo_ << "not yet supported"
                        << exit(FatalError);

                    BEGIN(readReactionKeyword);
                    break;
                }

                case collisionCrossSection:
                {
                    FatalErrorInFunction
                        << "Collision cross-section in reaction on line "
                        << lineNo_ << "not yet supported"
                        << exit(FatalError);

                    BEGIN(readReactionKeyword);
                    break;
                }

                case end:
                {
                    BEGIN(INITIAL);
                    addReaction
                    (
                        lhs,
                        rhs,
                        thirdBodyEfficiencies,
                        rType,
                        rrType,
                        fofType,
                        ArrheniusCoeffs,
                        reactionCoeffsTable,
                        RRreaction
                    );
                    finishReaction = false;
                    rType = unknownReactionType;
                    rrType = Arrhenius;
                    fofType = unknownFallOffFunctionType;
                    thirdBodyEfficiencies = 1.0;
                    pDependentSpecieName = word::null;
                    lrhsPtr = &lhs;
                    break;
                }

                default:
                {
                    FatalErrorInFunction
                        << "keyword " << keyword
                        << " should be followed by parameters"
                        << " on line " << lineNo_
                        << exit(FatalError);
                }
            }
        }
        else
        {
            currentSpecieName = foamName(foamSpecieString(YYText()));

            HashTable<label>::iterator specieIndexIter
            (
                specieIndices_.find(currentSpecieName)
            );

            if (specieIndexIter != specieIndices_.end())
            {
                if (finishReaction)
                {
                    addReaction
                    (
                        lhs,
                        rhs,
                        thirdBodyEfficiencies,
                        rType,
                        rrType,
                        fofType,
                        ArrheniusCoeffs,
                        reactionCoeffsTable,
                        RRreaction
                    );
                    finishReaction = false;
                    rType = unknownReactionType;
                    rrType = Arrhenius;
                    fofType = unknownFallOffFunctionType;
                    thirdBodyEfficiencies = 1.0;
                    pDependentSpecieName = word::null;
                    lrhsPtr = &lhs;
                }

                currentSpecieCoeff.index = specieIndexIter();
                lrhsPtr->append(currentSpecieCoeff);

                BEGIN(readReactionDelimiter);
            }
            else
            {
                BEGIN(readSpecieNamePlus);
            }
        }
    }

<readReactionKeyword>{reactionKeywordSlash} {

        word keyword(foamName(YYText()));

        HashTable<int>::iterator reactionKeywordIter
        (
            reactionKeywordTable_.find(keyword)
        );

        if (reactionKeywordIter != reactionKeywordTable_.end())
        {
            switch(reactionKeywordIter())
            {
                case unimolecularFallOffReactionType:
                {
                    if (!pDependentSpecieName.size())
                    {
                        FatalErrorInFunction
                            << "LOW keyword given for a unimolecular fall-off"
                               " reaction which does not contain a pressure"
                               " dependent specie" << " on line " << lineNo_
                            << exit(FatalError);
                    }

                    if (rrType == Arrhenius)
                    {
                        rrType = unimolecularFallOff;
                    }
                    else
                    {
                        FatalErrorInFunction
                            << "Attempt to set reaction rate type to"
                               " unimolecularFallOff when it is already set to "
                            << reactionRateTypeNames[rrType]
                            << " on line " << lineNo_
                            << exit(FatalError);
                    }

                    if (fofType == unknownFallOffFunctionType)
                    {
                        fofType = Lindemann;
                    }

                    reactionCoeffsName = reactionRateTypeNames[rrType];
                    BEGIN(readReactionCoeffs);
                    break;
                }

                case chemicallyActivatedBimolecularReactionType:
                {
                    if (!pDependentSpecieName.size())
                    {
                        FatalErrorInFunction
                            << "HIGH keyword given for a chemically"
                               " activated bimolecular reaction which does not"
                               " contain a pressure dependent specie"
                            << " on line " << lineNo_
                            << exit(FatalError);
                    }

                    if (rrType == Arrhenius)
                    {
                        rrType = chemicallyActivatedBimolecular;
                    }
                    else
                    {
                        FatalErrorInFunction
                            << "Attempt to set reaction rate type to"
                               " chemicallyActivatedBimolecular when it is"
                               " already set to "
                            << reactionRateTypeNames[rrType]
                            << " on line " << lineNo_
                            << exit(FatalError);
                    }

                    if (fofType == unknownFallOffFunctionType)
                    {
                        fofType = Lindemann;
                    }

                    reactionCoeffsName = reactionRateTypeNames[rrType];
                    BEGIN(readReactionCoeffs);
                    break;
                }

                case TroeReactionType:
                {
                    if (!pDependentSpecieName.size())
                    {
                        FatalErrorInFunction
                            << "TROE keyword given for a"
                               " reaction which does not contain a pressure"
                               " dependent specie" << " on line " << lineNo_
                            << exit(FatalError);
                    }

                    if
                    (
                        fofType == unknownFallOffFunctionType
                     || fofType == Lindemann
                    )
                    {
                        fofType = Troe;
                    }
                    else
                    {
                        FatalErrorInFunction
                            << "Attempt to set fall-off function type to Troe"
                               " when it is already set to "
                            << fallOffFunctionNames[fofType]
                            << " on line " << lineNo_
                            << exit(FatalError);
                    }

                    reactionCoeffsName = fallOffFunctionNames[fofType];
                    BEGIN(readReactionCoeffs);
                    break;
                }

                case SRIReactionType:
                {
                    if (!pDependentSpecieName.size())
                    {
                        FatalErrorInFunction
                            << "SRI keyword given for a"
                               " reaction which does not contain a pressure"
                               " dependent specie" << " on line " << lineNo_
                            << exit(FatalError);
                    }

                    if
                    (
                        fofType == unknownFallOffFunctionType
                     || fofType == Lindemann
                    )
                    {
                        fofType = SRI;
                    }
                    else
                    {
                        FatalErrorInFunction
                            << "Attempt to set fall-off function type to SRI"
                               " when it is already set to "
                            << fallOffFunctionNames[fofType]
                            << " on line " << lineNo_
                            << exit(FatalError);
                    }

                    reactionCoeffsName = fallOffFunctionNames[fofType];
                    BEGIN(readReactionCoeffs);
                    break;
                }

                case LandauTellerReactionType:
                {
                    if (pDependentSpecieName.size())
                    {
                        FatalErrorInFunction
                            << "Landau-Teller reaction rate cannot be used"
                               " for the pressure-dependent reaction on line "
                            << lineNo_
                            << exit(FatalError);
                    }

                    if (rrType == Arrhenius)
                    {
                        rrType = LandauTeller;
                    }
                    else
                    {
                        FatalErrorInFunction
                            << "Attempt to set reaction rate type to"
                               " LandauTeller when it is already set to "
                            << reactionRateTypeNames[rrType]
                            << " on line " << lineNo_
                            << exit(FatalError);
                    }

                    rrType = LandauTeller;
                    reactionCoeffsName = reactionRateTypeNames[rrType];
                    BEGIN(readReactionCoeffs);
                    break;
                }

                case reverseLandauTellerReactionType:
                {
                    if (pDependentSpecieName.size())
                    {
                        FatalErrorInFunction
                            << "Non-equilibrium Landau-Teller reaction rate"
                               " cannot be used"
                               " for the pressure-dependent reaction on line "
                            << lineNo_
                            << exit(FatalError);
                    }

                    if (rType != nonEquilibriumReversible)
                    {
                        FatalErrorInFunction
                            << "Reverse reaction Arrhenius coefficients not"
                               " given for reverse LandauTeller reaction."
                               " Please reorder 'REV' keyword to precede 'RLT'"
                            << " on line " << lineNo_
                            << exit(FatalError);
                    }

                    rrType = LandauTeller;
                    reactionCoeffsName =
                        word(reactionTypeNames[rType])
                      + reactionRateTypeNames[rrType];
                    BEGIN(readReactionCoeffs);
                    break;
                }

                case JanevReactionType:
                {
                    if (rrType == Arrhenius)
                    {
                        rrType = Janev;
                    }
                    else
                    {
                        FatalErrorInFunction
                            << "Attempt to set reaction rate type to"
                               " Janev when it is already set to "
                            << reactionRateTypeNames[rrType]
                            << " on line " << lineNo_
                            << exit(FatalError);
                    }

                    reactionCoeffsName = reactionRateTypeNames[rrType];
                    BEGIN(readReactionCoeffs);
                    break;
                }

                case powerSeriesReactionRateType:
                {
                    if (rrType == Arrhenius)
                    {
                        rrType = powerSeries;
                    }
                    else
                    {
                        FatalErrorInFunction
                            << "Attempt to set reaction rate type to"
                               " powerSeries when it is already set to "
                            << reactionRateTypeNames[rrType]
                            << " on line " << lineNo_
                            << exit(FatalError);
                    }

                    reactionCoeffsName = reactionRateTypeNames[rrType];
                    BEGIN(readReactionCoeffs);
                    break;
                }

                case radiationActivatedReactionType:
                {
                    FatalErrorInFunction
                        << "Radiation activated reaction on line "
                        << lineNo_ << "not yet supported"
                        << exit(FatalError);
                    // reactionCoeffsName = reactionRateTypeNames[rrType];
                    BEGIN(readReactionCoeffs);
                    break;
                }

                case energyLossReactionType:
                {
                    FatalErrorInFunction
                        << "Energy loss in reaction on line "
                        << lineNo_ << "not yet supported"
                        << exit(FatalError);
                    // reactionCoeffsName = reactionRateTypeNames[rrType];
                    BEGIN(readReactionCoeffs);
                    break;
                }

                case nonEquilibriumReversibleReactionType:
                {
                    rType = nonEquilibriumReversible;
                    reactionCoeffsName = reactionTypeNames[rType];
                    BEGIN(readReactionCoeffs);
                    break;
                }

                case speciesOrderForward:
                {
                    lrhsPtr = &lhs;
                    BEGIN(readReactionOrderSpecie);
                    break;
                }

                case speciesOrderReverse:
                {
                    lrhsPtr = &rhs;
                    BEGIN(readReactionOrderSpecie);
                    break;
                }

                case UnitsOfReaction:
                {
                    BEGIN(readReactionUnit);
                    break;
                }

                default:
                {
                    FatalErrorInFunction
                        << "unknown reaction keyword " << keyword
                        << " on line " << lineNo_
                        << exit(FatalError);
                }
            }
        }
        else
        {
            HashTable<label>::iterator specieIndexIter
            (
                specieIndices_.find(keyword)
            );

            if (specieIndexIter != specieIndices_.end())
            {
                currentThirdBodyIndex = specieIndexIter();
            }
            else
            {
                FatalErrorInFunction
                    << "unknown third-body specie " << keyword
                    << " on line " << lineNo_ << nl
                    << "Valid species are : " << nl
                    << specieIndices_.toc() << endl
                    << exit(FatalError);
            }

            BEGIN(readThirdBodyEfficiency);
        }
    }

<readSpecieNamePlus>"+" {
        currentSpecieName += "+";

        HashTable<label>::iterator specieIndexIter
        (
            specieIndices_.find(currentSpecieName)
        );

        if (specieIndexIter != specieIndices_.end())
        {
            if (finishReaction)
            {
                addReaction
                (
                    lhs,
                    rhs,
                    thirdBodyEfficiencies,
                    rType,
                    rrType,
                    fofType,
                    ArrheniusCoeffs,
                    reactionCoeffsTable,
                    RRreaction
                );
                finishReaction = false;
                rType = unknownReactionType;
                rrType = Arrhenius;
                fofType = unknownFallOffFunctionType;
                thirdBodyEfficiencies = 1.0;
                pDependentSpecieName = word::null;
                lrhsPtr = &lhs;
            }

            currentSpecieCoeff.index = specieIndexIter();
            lrhsPtr->append(currentSpecieCoeff);

            BEGIN(readReactionDelimiter);
        }
        else
        {
            FatalErrorInFunction
                << "unknown specie " << currentSpecieName
                << " on line " << lineNo_ << nl
                << "Valid species are : " << nl
                << specieIndices_.toc() << endl
                << exit(FatalError);
        }
    }

<readReactionDelimiter>{specieDelimiter} {
        currentSpecieCoeff.stoichCoeff = 1.0;
        currentSpecieCoeff.exponent = currentSpecieCoeff.stoichCoeff;
        BEGIN(readReactionKeyword);
    }

<readReactionDelimiter>{irreversibleReactionDelimiter} {
        currentSpecieCoeff.stoichCoeff = 1.0;
        currentSpecieCoeff.exponent = currentSpecieCoeff.stoichCoeff;
        rType = irreversible;
        lrhsPtr = &rhs;
        BEGIN(readReactionKeyword);
    }

<readReactionDelimiter>{reversibleReactionDelimiter} {
        currentSpecieCoeff.stoichCoeff = 1.0;
        currentSpecieCoeff.exponent = currentSpecieCoeff.stoichCoeff;
        rType = reversible;
        lrhsPtr = &rhs;
        BEGIN(readReactionKeyword);
    }

<readReactionDelimiter>& {
    }

<readReactionDelimiter>{reactionCoeffs} {
        string reactionCoeffsString(YYText());
        reactionCoeffsString.replaceAll("d", "e");
        reactionCoeffsString.replaceAll("D", "e");
        IStringStream reactionCoeffsStream(reactionCoeffsString);
        reactionCoeffsStream.lineNumber() = lineNo_;

        reactionCoeffsStream
            >> ArrheniusCoeffs[0]
            >> ArrheniusCoeffs[1]
            >> ArrheniusCoeffs[2];

        finishReaction = true;
        currentSpecieCoeff.stoichCoeff = 1.0;
        currentSpecieCoeff.exponent = currentSpecieCoeff.stoichCoeff;
        RRreaction = RRreactions;

        if (lhsThirdBodyCounter || rhsThirdBodyCounter)
        {
            if (!lhsThirdBodyCounter || !rhsThirdBodyCounter)
            {
                FatalErrorInFunction
                    << "Third body not present on both sides of reaction"
                       " on line " << lineNo_
                << exit(FatalError);
            }

            if (lhsThirdBodyCounter != 1)
            {
                FatalErrorInFunction
                    << "More than 1 third body present on l.h.s. side"
                       " of reaction on line " << lineNo_
                    << exit(FatalError);
            }

            if (rhsThirdBodyCounter != 1)
            {
                FatalErrorInFunction
                    << "More than 1 third body present on r.h.s. side"
                       " of reaction on line " << lineNo_
                    << exit(FatalError);
            }

            lhsThirdBodyCounter = 0;
            rhsThirdBodyCounter = 0;
        }

        BEGIN(readReactionKeyword);
    }


<readThirdBodyEfficiency>{thirdBodyEfficiency} {
        thirdBodyEfficiencies[currentThirdBodyIndex] = stringToScalar(YYText());
        BEGIN(readReactionKeyword);
    }

<readReactionDelimiter>{startPDependentSpecie} {
        BEGIN(readPDependentSpecie);
    }

<readPDependentSpecie>{pDependentSpecie} {

        word rhsPDependentSpecieName = pDependentSpecieName;
        pDependentSpecieName =
            foamName(foamSpecieString(YYText()));
        pDependentSpecieName =
            pDependentSpecieName(0, pDependentSpecieName.size() - 1);

        if (rrType == thirdBodyArrhenius)
        {
            FatalErrorInFunction
                << "The pressure-dependent third-body '"
                << pDependentSpecieName
                << "' is given in non-pressure-dependent third-body reaction"
                << " on line " << lineNo_
                << exit(FatalError);
        }

        if (lrhsPtr == &lhs)
        {
            lhsThirdBodyCounter++;
        }
        else
        {
            if (pDependentSpecieName != rhsPDependentSpecieName)
            {
                FatalErrorInFunction
                    << "The third-body reactant '"
                    << pDependentSpecieName
                    << "' is not the same as the third-body product '"
                    << rhsPDependentSpecieName
                    << "' in pressure-dependent reaction on line " << lineNo_
                    << exit(FatalError);
            }

            rhsThirdBodyCounter++;
        }

        if (pDependentSpecieName != "M")
        {
            HashTable<label>::iterator specieIndexIter
            (
                specieIndices_.find(pDependentSpecieName)
            );

            if (specieIndexIter != specieIndices_.end())
            {
                thirdBodyEfficiencies = 0.0;
                thirdBodyEfficiencies[specieIndexIter()] = 1.0;
            }
            else
            {
                FatalErrorInFunction
                    << "unknown third-body specie " << pDependentSpecieName
                    << " on line " << lineNo_ << nl
                    << "Valid species are : " << nl
                    << specieIndices_.toc() << endl
                    << exit(FatalError);
            }
        }

        BEGIN(readReactionDelimiter);
    }

<readReactionCoeffs>{reactionCoeff} {
        reactionCoeffs.append(stringToScalar(YYText()));
    }

<readReactionCoeffs>{endReactionCoeffs} {
        reactionCoeffsTable.insert(reactionCoeffsName, reactionCoeffs.shrink());
        reactionCoeffs.clear();
        BEGIN(readReactionKeyword);
    }

<readTdepSpecie>{specieName} {
        word specieName(foamName(foamSpecieString(YYText())));
        FatalErrorInFunction
            << "Temperature-dependent reaction on line "
            << lineNo_ << "not yet supported"
            << exit(FatalError);
        BEGIN(readReactionKeyword);
    }

<readReactionOrderSpecie>{specieName} {
        currentSpecieName =
            word(foamName(foamSpecieString(YYText())));

        HashTable<label>::iterator specieIndexIter
        (
            specieIndices_.find(currentSpecieName)
        );

        if (specieIndexIter != specieIndices_.end())
        {
            currentSpecieIndex = specieIndexIter();
        }
        else
        {
            FatalErrorInFunction
                << "unknown specie " << currentSpecieName
                << " given in reaction-order specification"
                << " on line " << lineNo_ << nl
                << "Valid species are : " << nl
                << specieIndices_.toc() << endl
                << exit(FatalError);
        }

        BEGIN(readReactionOrder);
    }

<readReactionOrder>{reactionCoeff}{endReactionCoeffs} {

        DynamicList<gasHReaction::specieCoeffs>& lrhs = *lrhsPtr;

        bool found = false;

        forAll(lrhs, i)
        {
            if (lrhs[i].index == currentSpecieIndex)
            {
                lrhs[i].exponent = stringToScalar(YYText());
                found = true;
                break;
            }
        }

        if (!found)
        {
            word side("l.h.s.");

            if (lrhsPtr == &rhs)
            {
                side = "r.h.s.";
            }

            FatalErrorInFunction
                << "Specie " << currentSpecieName
                << " on line " << lineNo_
                << " not present in " << side << " of reaction " << nl << lrhs
                << exit(FatalError);
        }

        BEGIN(readReactionKeyword);
    }

<readReactionUnit>{cal} {
        RRreaction = RRcal;
    }

<readReactionUnit>{kcal} {
        RRreaction = RRcal/1000.0;
    }

<readReactionUnit>{joule} {
        RRreaction = RRjoule;
    }

<readReactionUnit>{kjoule} {
        RRreaction = RRjoule/1000.0;
    }

<readReactionUnit>{kelvin} {
        RRreaction = 1.0;
    }

<readReactionUnit>{otherReactionsUnit} {
    }

<readReactionUnit>{endReactionCoeffs} {
        BEGIN(readReactionKeyword);
    }


 /* ------------------ Ignore remaining space and \n s. --------------------- */

<*>\n      { lineNo_++; }

 /* ------ Ignore remaining space and \n s.  Any other characters are errors. */

<*>. {
        startError = YYText();
        yy_push_state(CHEMKINError);
    }

<CHEMKINError>.* {
        yy_pop_state();
        FatalErrorInFunction
            << "while " << stateNames[YY_START] << " on line " << lineNo_ << nl
            << "    expected " << stateExpects[YY_START]
            << " but found '" << startError << YYText() << "'"
            << exit(FatalError);
    }


 /*  ------------------------ On EOF terminate. ----------------------------  */

<<EOF>> {
        if (finishReaction)
        {
            addReaction
            (
                lhs,
                rhs,
                thirdBodyEfficiencies,
                rType,
                rrType,
                fofType,
                ArrheniusCoeffs,
                reactionCoeffsTable,
                RRreaction
            );
        }
        yyterminate();
    }
%%

 /* ------------------------------------------------------------------------- *\
    ------ End of chemkinLexer.L
 \* ------------------------------------------------------------------------- */

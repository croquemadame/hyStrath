/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2021 hyStrath
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of hyStrath, a derivative work of OpenFOAM.

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

#include "specieReactionRates.H"
#include "volFields.H"
#include "fvcVolumeIntegrate.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ChemistryModelType>
void Foam::functionObjects::specieReactionRates<ChemistryModelType>::
writeFileHeader
(
    Ostream& os
) const
{
    writeHeader(os, "Specie reaction rates");
    volRegion::writeFileHeader(*this, os);
    writeHeaderValue(os, "nSpecie", chemistryModel_.nSpecie());
    writeHeaderValue(os, "nReaction", chemistryModel_.nReaction());

    writeCommented(os, "Time");
    writeTabbed(os, "Reaction");

    const wordList& speciesNames =
        chemistryModel_.thermo().composition().species();

    forAll (speciesNames, si)
    {
        writeTabbed(os, speciesNames[si]);
    }

    os  << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ChemistryModelType>
Foam::functionObjects::specieReactionRates<ChemistryModelType>::
specieReactionRates
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    volRegion(fvMeshFunctionObject::mesh_, dict),
    writeFile(obr_, name, typeName, dict),
    chemistryModel_
    (
        fvMeshFunctionObject::mesh_.lookupObject<ChemistryModelType>
        (
            "chemistryProperties"
        )
    )
{
    writeFileHeader(file());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ChemistryModelType>
Foam::functionObjects::specieReactionRates<ChemistryModelType>::
~specieReactionRates()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ChemistryModelType>
bool Foam::functionObjects::specieReactionRates<ChemistryModelType>::read
(
    const dictionary& dict
)
{
    regionFunctionObject::read(dict);

    return true;
}


template<class ChemistryModelType>
bool Foam::functionObjects::specieReactionRates<ChemistryModelType>::execute()
{
    return true;
}


template<class ChemistryModelType>
bool Foam::functionObjects::specieReactionRates<ChemistryModelType>::write()
{
    const label nSpecie = chemistryModel_.nSpecie();
    const label nReaction = chemistryModel_.nReaction();

    volRegion::update();        // Ensure cached values are valid

    // Region volume
    //const scalar V = this->V();

    const scalar volTotal = this->volRegion::V();
    const bool useAll = this->volRegion::useAllCells();

    for (label ri=0; ri<nReaction; ri++)
    {
        writeCurrentTime(file());
        file() << token::TAB << ri;

        for (label si=0; si<nSpecie; si++)
        {
            volScalarField::Internal RR
            (
                chemistryModel_.calculateRR(ri, si)
            );

            scalar sumVRRi = 0;

            if (useAll)
            {
                sumVRRi = fvc::domainIntegrate(RR).value();
            }
            else
            {
                sumVRRi = gSum
                (
                    scalarField(fvMeshFunctionObject::mesh_.V()*RR, cellIDs())
                );
            }

            file() << token::TAB << sumVRRi/volTotal;
        }

        file() << nl;
    }

    file() << nl << endl;

    return true;
}


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

#include "addToRunTimeSelectionTable.H"
#include "rho2ChemistryModel.H"

namespace Foam
{
    typedef functionObjects::specieReactionRates<rho2ChemistryModel>
        rhoSpecieReactionRates;

    defineTemplateTypeNameAndDebugWithName
    (
        rhoSpecieReactionRates,
        "rhoSpecieReactionRates",
        0
    );

    addToRunTimeSelectionTable
    (
        functionObject,
        rhoSpecieReactionRates,
        dictionary
    );
}


// ************************************************************************* //

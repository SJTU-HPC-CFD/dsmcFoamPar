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

Class
    dsmcDynamicLoadBalancing

Description

\*----------------------------------------------------------------------------*/

#include "dsmcDynamicLoadBalancing.H"
#include "processorPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "wallPolyPatch.H"
#include "dsmcCloud.H"

#include <chrono>

namespace Foam
{

namespace
{
scalar elapsedSeconds(const std::chrono::steady_clock::time_point& start)
{
    return std::chrono::duration_cast<std::chrono::duration<scalar>>
    (
        std::chrono::steady_clock::now() - start
    ).count();
}

scalar lookupBalanceUntilTime(const dictionary& dict)
{
    if (dict.found("balanceUntilTime"))
    {
        return dict.lookupOrDefault<scalar>("balanceUntilTime", VGREAT);
    }

    return dict.lookupOrDefault<scalar>("loadBalancingUntilTime", VGREAT);
}
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

scalar dsmcDynamicLoadBalancing::totalBalanceWallTime_ = 0.0;
scalar dsmcDynamicLoadBalancing::totalPrepareWallTime_ = 0.0;
scalar dsmcDynamicLoadBalancing::totalReconstructMeshWallTime_ = 0.0;
scalar dsmcDynamicLoadBalancing::totalReconstructFieldsWallTime_ = 0.0;
scalar dsmcDynamicLoadBalancing::totalDecomposeWallTime_ = 0.0;
scalar dsmcDynamicLoadBalancing::totalProcessorMeshSyncWallTime_ = 0.0;
scalar dsmcDynamicLoadBalancing::totalBackupWallTime_ = 0.0;
label dsmcDynamicLoadBalancing::totalBalanceCount_ = 0;
bool dsmcDynamicLoadBalancing::reportTimingEnabled_ = false;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Constructor
dsmcDynamicLoadBalancing::dsmcDynamicLoadBalancing
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud
)
:
    time_(t),
    mesh_(refCast<const fvMesh>(mesh)),
    cloud_(cloud),
    dsmcLoadBalanceDict_
    (
        IOobject
        (
            "loadBalanceDict",
            time_.system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    performBalance_(false),
    enableBalancing_(Switch(Foam::hyCompat::lookup(dsmcLoadBalanceDict_, "enableBalancing"))),
    reportTiming_
    (
        dsmcLoadBalanceDict_.lookupOrDefault<Switch>
        (
            "reportBalancingTiming",
            enableBalancing_
        )
    ),
    balanceUntilTime_(lookupBalanceUntilTime(dsmcLoadBalanceDict_)),
    originalEndTime_(time_.time().endTime().value()),
    maxImbalance_(Foam::hyCompat::toScalar
    (
        dsmcLoadBalanceDict_, "maximumAllowableImbalance"
    )),
    balanceCheckInterval_
    (
        max
        (
            dsmcLoadBalanceDict_.lookupOrDefault<label>
            (
                "balanceCheckInterval",
                1
            ),
            label(1)
        )
    ),
    lastBalanceCheckTimeIndex_(-1),
    limitTimeDirBackups_
    (
        dsmcLoadBalanceDict_.lookupOrDefault<label>
        (
            "limitTimeDirBackups",
            -1
        )
    )
{
    reportTimingEnabled_ = reportTimingEnabled_ || reportTiming_;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcDynamicLoadBalancing::~dsmcDynamicLoadBalancing()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcDynamicLoadBalancing::update()
{
    updateProperties();

    if (!enableBalancing_)
    {
        performBalance_ = false;
        return;
    }

    if (!Pstream::parRun())
    {
        return;
    }

    const label currentTimeIndex = time_.timeIndex();

    if
    (
        currentTimeIndex == lastBalanceCheckTimeIndex_
     || (currentTimeIndex % balanceCheckInterval_) != 0
    )
    {
        return;
    }

    lastBalanceCheckTimeIndex_ = currentTimeIndex;

    const scalar& allowableImbalance = maxImbalance_;

    // First determine current level of imbalance - do this for all
    // parallel runs, even if balancing is disabled.
    scalar nGlobalParticles = cloud_.size();
    Foam::reduce(nGlobalParticles, sumOp<scalar>());

    scalar idealNParticles =
        scalar(nGlobalParticles)/scalar(Pstream::nProcs());

    scalar nParticles = cloud_.size();
    scalar localImbalance = mag(nParticles - idealNParticles);
    Foam::reduce(localImbalance, maxOp<scalar>());
    scalar maxImbalance = localImbalance/idealNParticles;

    Info<< "    DLB imbalance check at time " << time_.timeName()
        << " (timeIndex " << currentTimeIndex
        << ", interval " << balanceCheckInterval_ << ")" << nl
        << "    Maximum imbalance = " << 100*maxImbalance << "%" << nl
        << endl;

    // explanation of modes:
    // 1. enableBalancing = true:
    //   1. if time <= balanceUntilTime -> load balance
    //   2. if time > balanceUntilTime -> do not load balance
    // 2. enableBalancing = false -> do not load balance
    if
    (
           enableBalancing_
        && time_.time().value() <= balanceUntilTime_
        && maxImbalance > allowableImbalance
    )
    {
        Info<< "    DLB trigger: forcing write of current time before "
            << "mesh repartition" << nl << endl;

        writeParticleWeightField();
        time_.writeNow();

        performBalance_ = true;

        originalEndTime_ = time_.time().endTime().value();

        scalar currentTime = time_.time().value();

        time_.setEndTime(currentTime);
    }
}


void dsmcDynamicLoadBalancing::copyPolyMeshToLatestTimeFolder() const
{
    for (label i=0; i<Pstream::nProcs(); i++)
    {
        const word processorName = "processor" + name(i) + "/";
        const word copyPolyMesh =
            word("latest=$(find ")
          + processorName
          + word(" -mindepth 1 -maxdepth 1 -type d | sed 's#.*/##' | ")
          + word("grep -E '^[0-9.+-eE]+$' | sort -g | tail -n 1);")
          + word("starting=$(find ")
          + processorName
          + word(" -mindepth 1 -maxdepth 1 -type d | sed 's#.*/##' | ")
          + word("grep -E '^[0-9.+-eE]+$' | sort -g | head -n 1);")
          + word("if [ -n \"$latest\" ] && [ -n \"$starting\" ] && ")
          + word("[ ! -d ")
          + processorName
          + word("$latest/polyMesh ]; then meshSource=")
          + processorName
          + word("$starting/polyMesh; ")
          + word("if [ ! -d \"$meshSource\" ]; then meshSource=")
          + processorName
          + word("constant/polyMesh; fi; mkdir -p ")
          + processorName
          + word("$latest; cp -r \"$meshSource\" ")
          + processorName
          + word("$latest/; fi");

        Foam::system(copyPolyMesh);
    }
}


void dsmcDynamicLoadBalancing::writeParticleWeightField() const
{
    volScalarField particleWeights
    (
        IOobject
        (
            "dsmcParticleCount",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    );

    forAllConstIter(dsmcCloud, cloud_, iter)
    {
        const label celli = iter().cell();

        if (celli >= 0 && celli < particleWeights.size())
        {
            particleWeights[celli] += 1.0;
        }
    }

    particleWeights.correctBoundaryConditions();
    particleWeights.write();
}


void dsmcDynamicLoadBalancing::perform(const label noRefinement)
{
    if (enableBalancing_ && performBalance_)
    {
        const auto balanceStart = std::chrono::steady_clock::now();
        scalar prepareWallTime = 0.0;
        scalar reconstructMeshWallTime = 0.0;
        scalar reconstructFieldsWallTime = 0.0;
        scalar decomposeWallTime = 0.0;
        scalar processorMeshSyncWallTime = 0.0;
        scalar backupWallTime = 0.0;

        if (Pstream::master())
        {
            if (noRefinement == 0)
            {
                const auto sectionStart = std::chrono::steady_clock::now();
                copyPolyMeshToLatestTimeFolder();
                prepareWallTime += elapsedSeconds(sectionStart);
            }

            {
                const auto sectionStart = std::chrono::steady_clock::now();
                Foam::system("reconstructParMesh -latestTime");
                reconstructMeshWallTime = elapsedSeconds(sectionStart);
            }

            {
                const auto sectionStart = std::chrono::steady_clock::now();
                Foam::system("reconstructPar -latestTime");
                reconstructFieldsWallTime = elapsedSeconds(sectionStart);
            }

            {
                const auto sectionStart = std::chrono::steady_clock::now();
                Foam::system
                (
                    "latest=$(foamListTimes -latestTime);"
                    "if [ ! -d $latest/polyMesh ]; then cp -r constant/polyMesh $latest/; fi"
                );
                Foam::system("rm -r processor*");
                prepareWallTime += elapsedSeconds(sectionStart);
            }

            const word decomposePar =
                word("timeDirs=$(foamListTimes -noZero);")
                // check if there are any time dirs, if not this indicates a
                // fatal error
                + word("if [ -z \"$timeDirs\" ];")
                + word("then echo \"error: no time directories to decompose\";")
                // decompose the latest time dir
                + word("else decomposeTool=\"${FOAM_USER_APPBIN}/decomposeDSMCLoadBalancePar\";")
                + word("if [ ! -x \"$decomposeTool\" ]; then ")
                + word("decomposeTool=decomposeDSMCLoadBalancePar; fi;")
                + word("\"$decomposeTool\" -force -latestTime -copyUniform;")
                + word("fi");
            {
                const auto sectionStart = std::chrono::steady_clock::now();
                Foam::system(decomposePar);
                decomposeWallTime = elapsedSeconds(sectionStart);
            }

            {
                const auto sectionStart = std::chrono::steady_clock::now();
                Foam::system
                (
                    "latest=$(foamListTimes -processor -latestTime); "
                    "for procDir in processor*; do "
                    "if [ -d \"$procDir/$latest/polyMesh\" ]; then "
                    "mkdir -p \"$procDir/constant\"; "
                    "rm -rf \"$procDir/constant/polyMesh\"; "
                    "cp -r \"$procDir/$latest/polyMesh\" \"$procDir/constant/\"; "
                    "fi; "
                    "done"
                );
                processorMeshSyncWallTime = elapsedSeconds(sectionStart);
            }

            // backup time dirs must be stored in resultFolders to prevent them
            // from being cleared. They can be moved back when the simulation
            // has finished.
            mkDir("resultFolders");

            // respect the specified limit of max. concurrent time dir backups
            if (limitTimeDirBackups_ >= 0)
            {
                // impose the limit by moving back the time dirs that are
                // currently already backuped up and then only keeping the
                // most recent time dirs <= the imposed limit.
                const word backupTimeDirsWithLimit =
                    // move the time dirs currently backed up to the case dir
                    // so the foamListTimes utility can be used
                    // make sure directory is not empty to prevent mv from
                    // printing a warning
                    word("if [ \"$(ls -A resultFolders)\" ];")
                    + word("then mv resultFolders/* .; fi;")
                    // total number of time directories
                    + word("timeDirs=$(foamListTimes);")
                    + word("nTimeDirs=$(echo $timeDirs | tr -cd ' ' | wc -c);")
                    + word("nTimeDirs=$((nTimeDirs+1));")
                    // convert OpenFOAM label to shell variable for limit
                    + word("limitNTimeDirs=")
                    + name(limitTimeDirBackups_)
                    + word(";")
                    // if the total number of time directories is larger than
                    // limit remove the oldest time directories first
                    + word("diffNTimeDirs=$((nTimeDirs-limitNTimeDirs));")
                    + word("if [ $diffNTimeDirs -gt 0 ];")
                    + word("then timeDirs=$(echo $timeDirs | ")
                    + word("cut -d ' ' -f$((diffNTimeDirs+1))-$nTimeDirs);")
                    + word("fi;")
                    // move the time dirs that were chosen to be eligible for
                    // backup to the backup dir
                    + word("mv $timeDirs resultFolders/;")
                    // clear all other time dirs
                    + word("foamListTimes -rm");
                const auto sectionStart = std::chrono::steady_clock::now();
                Foam::system(backupTimeDirsWithLimit);
                backupWallTime = elapsedSeconds(sectionStart);
            }
            else
            {
                // keep all time dirs in the backup dir
                const auto sectionStart = std::chrono::steady_clock::now();
                Foam::system
                (
                    "timeDirs=$(foamListTimes); mv $timeDirs resultFolders/"
                );
                backupWallTime = elapsedSeconds(sectionStart);
            }

            performBalance_ = false;
        }

        reduce(prepareWallTime, maxOp<scalar>());
        reduce(reconstructMeshWallTime, maxOp<scalar>());
        reduce(reconstructFieldsWallTime, maxOp<scalar>());
        reduce(decomposeWallTime, maxOp<scalar>());
        reduce(processorMeshSyncWallTime, maxOp<scalar>());
        reduce(backupWallTime, maxOp<scalar>());

        const scalar balanceWallTime = elapsedSeconds(balanceStart);

        totalBalanceWallTime_ += balanceWallTime;
        totalPrepareWallTime_ += prepareWallTime;
        totalReconstructMeshWallTime_ += reconstructMeshWallTime;
        totalReconstructFieldsWallTime_ += reconstructFieldsWallTime;
        totalDecomposeWallTime_ += decomposeWallTime;
        totalProcessorMeshSyncWallTime_ += processorMeshSyncWallTime;
        totalBackupWallTime_ += backupWallTime;
        totalBalanceCount_++;

        if (reportTiming_ && Pstream::master())
        {
            Info<< "DLB wall-time summary [s]:" << nl
                << "    total              = " << balanceWallTime << nl
                << "    prepare            = " << prepareWallTime << nl
                << "    reconstruct mesh   = " << reconstructMeshWallTime << nl
                << "    reconstruct fields = " << reconstructFieldsWallTime << nl
                << "    decompose          = " << decomposeWallTime << nl
                << "    processor mesh sync= " << processorMeshSyncWallTime << nl
                << "    backup             = " << backupWallTime << nl
                << endl;
        }

        time_.setEndTime(originalEndTime_);
    }
}


void dsmcDynamicLoadBalancing::resetTiming()
{
    totalBalanceWallTime_ = 0.0;
    totalPrepareWallTime_ = 0.0;
    totalReconstructMeshWallTime_ = 0.0;
    totalReconstructFieldsWallTime_ = 0.0;
    totalDecomposeWallTime_ = 0.0;
    totalProcessorMeshSyncWallTime_ = 0.0;
    totalBackupWallTime_ = 0.0;
    totalBalanceCount_ = 0;
    reportTimingEnabled_ = false;
}

void dsmcDynamicLoadBalancing::updateProperties()
{
    enableBalancing_ = Switch(Foam::hyCompat::lookup(Foam::hyCompat::lookup(dsmcLoadBalanceDict_, "enableBalancing")));
    reportTiming_ =
        dsmcLoadBalanceDict_.lookupOrDefault<Switch>
        (
            "reportBalancingTiming",
            enableBalancing_
        );
    reportTimingEnabled_ = reportTimingEnabled_ || reportTiming_;

    // if balancing is active this additional option allows to specify a time
    // after which balancing is deactivated (this is useful in conjunction with
    // resetAtOutput / resetAtOutputUntilTime and averaging across solver
    // restarts)
    balanceUntilTime_ = lookupBalanceUntilTime(dsmcLoadBalanceDict_);
    maxImbalance_ = Foam::hyCompat::toScalar
    (
        dsmcLoadBalanceDict_, "maximumAllowableImbalance"
    );
    balanceCheckInterval_ = max
    (
        dsmcLoadBalanceDict_.lookupOrDefault<label>
        (
            "balanceCheckInterval",
            1
        ),
        label(1)
    );
    limitTimeDirBackups_ = dsmcLoadBalanceDict_.lookupOrDefault<label>
    (
        "limitTimeDirBackups",
        -1
    );
}


}  // End namespace Foam

// ************************************************************************* //

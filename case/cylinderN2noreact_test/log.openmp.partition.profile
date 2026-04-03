/*---------------------------------------------------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2506                                  |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
Build  : 95a5dfacc7-20250714 OPENFOAM=2506 version=v2506
Arch   : "LSB;label=32;scalar=64"
Exec   : /home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/platforms/linux64IcxDPInt32Opt/bin/dsmcFoam+
Date   : Apr 02 2026
Time   : 21:03:41
Host   : SuperXCX-ROG
PID    : 2372554
I/O    : uncollated
Case   : /home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinderN2noreact_test
nProcs : 1
trapFpe: Floating point exception trapping enabled (FOAM_SIGFPE).
fileModificationChecking : Monitoring run-time modified files using timeStampMaster (fileModificationSkew 5, maxFileModificationPolls 20)
allowSystemOperations : Allowing user-supplied system call operations

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
Create time

--> FOAM Warning : 
    From static Foam::IOstreamOption::compressionType Foam::IOstreamOption::compressionEnum(const word &, const compressionType)
    in file db/options/IOstreamOption.C at line 189
    Unknown compression specifier 'uncompressed' using compression off
Create mesh for time = 0


Constructing dsmcCloud
Selecting the coordinate system model:	dsmcCartesian

Selecting the time-step model:	dsmcConstantTimeStepModel

Constant time-step model:
- time-step [sec]	6.632426865e-08

Selecting the porous measurement model:	dsmcNoPorousMediumMeasurements


Creating dsmcControllers

Creating fields: 

Selecting field: dsmcVolFields

TimeData Statistics: 
Resetting fields at output until time: 4.841671611e-06
 measurement option: write
 nSamples: 1, time interval: 6.632426865e-08
 nAverages: 83, time interval: 5.504914298e-06
 total no. of sampling steps: 0
 total no. of averaging Steps: 0


Creating the boundary models: 

Selecting dsmcPatchBoundaryModel dsmcDeletionPatch
Selecting dsmcPatchBoundaryModel dsmcDeletionPatch
Selecting dsmcPatchBoundaryModel dsmcDiffuseWallPatch
Selecting dsmcGeneralBoundaryModel dsmcFreeStreamInflowPatch
Selecting BinaryCollisionModel LarsenBorgnakkeVariableHardSphere
OpenMP enabled for dsmcCloud with 8 thread-local RNG streams using collision strategy 'partition'

Constructing constant properties for
    N2

Creating dsmcReactions

There are no chemical reactions defined.
Selecting collisionPartnerSelectionModel noTimeCounter

Initialising the measurement fields

Initialising dsmcVolFields field

Starting time loop

ExecutionTime = 9.8 s  ClockTime = 7 s

ExecutionTime = 13.24 s  ClockTime = 9 s

ExecutionTime = 16.3 s  ClockTime = 10 s

ExecutionTime = 19.68 s  ClockTime = 12 s

ExecutionTime = 22.57 s  ClockTime = 13 s

ExecutionTime = 25.79 s  ClockTime = 14 s

ExecutionTime = 28.96 s  ClockTime = 16 s

ExecutionTime = 32.17 s  ClockTime = 17 s

ExecutionTime = 35.43 s  ClockTime = 19 s

Time = 6.632426865e-07

    Collisions                      = 27886

Cloud name: dsmc
    Number of dsmc particles        = 1329919
    Number of molecules             = 7.979514e+15
    Mass in system                  = 3.71047401e-10
    Average linear kinetic energy   = 1.83630123e-19
    Average rotational energy       = 1.339602569e-21
    Average vibrational energy      = 2.327089075e-20
    Average electronic energy       = 0
    Average total energy            = 2.082406163e-19
ExecutionTime = 39.04 s  ClockTime = 21 s

ExecutionTime = 42.09 s  ClockTime = 22 s

ExecutionTime = 45.35 s  ClockTime = 24 s

ExecutionTime = 48.68 s  ClockTime = 25 s

ExecutionTime = 51.95 s  ClockTime = 27 s

ExecutionTime = 55.22 s  ClockTime = 28 s

ExecutionTime = 58.86 s  ClockTime = 30 s

ExecutionTime = 62.08 s  ClockTime = 32 s

ExecutionTime = 65.49 s  ClockTime = 33 s

ExecutionTime = 68.9 s  ClockTime = 35 s

Time = 1.326485373e-06

    Collisions                      = 37270

Cloud name: dsmc
    Number of dsmc particles        = 1353501
    Number of molecules             = 8.121006e+15
    Mass in system                  = 3.77626779e-10
    Average linear kinetic energy   = 1.810121667e-19
    Average rotational energy       = 1.576231594e-21
    Average vibrational energy      = 2.327089014e-20
    Average electronic energy       = 0
    Average total energy            = 2.058592884e-19
ExecutionTime = 72.63 s  ClockTime = 37 s

ExecutionTime = 75.88 s  ClockTime = 38 s

ExecutionTime = 79.57 s  ClockTime = 37 s

ExecutionTime = 83.2 s  ClockTime = 39 s

ExecutionTime = 87.13 s  ClockTime = 41 s

ExecutionTime = 90.6 s  ClockTime = 42 s

ExecutionTime = 94.24 s  ClockTime = 44 s

ExecutionTime = 98 s  ClockTime = 46 s

ExecutionTime = 101.66 s  ClockTime = 48 s

ExecutionTime = 105.37 s  ClockTime = 49 s

Time = 1.98972806e-06

    Collisions                      = 47843

Cloud name: dsmc
    Number of dsmc particles        = 1377013
    Number of molecules             = 8.262078e+15
    Mass in system                  = 3.84186627e-10
    Average linear kinetic energy   = 1.786571242e-19
    Average rotational energy       = 1.896478576e-21
    Average vibrational energy      = 2.327092335e-20
    Average electronic energy       = 0
    Average total energy            = 2.038245261e-19
ExecutionTime = 109.54 s  ClockTime = 51 s

ExecutionTime = 112.78 s  ClockTime = 53 s

ExecutionTime = 116.58 s  ClockTime = 55 s

ExecutionTime = 120.4 s  ClockTime = 57 s

ExecutionTime = 124.16 s  ClockTime = 58 s

ExecutionTime = 128.05 s  ClockTime = 60 s

ExecutionTime = 131.85 s  ClockTime = 62 s

ExecutionTime = 135.54 s  ClockTime = 64 s

ExecutionTime = 139.45 s  ClockTime = 66 s

ExecutionTime = 143.33 s  ClockTime = 68 s

Time = 2.652970746e-06

    Collisions                      = 57766

Cloud name: dsmc
    Number of dsmc particles        = 1400526
    Number of molecules             = 8.403156e+15
    Mass in system                  = 3.90746754e-10
    Average linear kinetic energy   = 1.765101772e-19
    Average rotational energy       = 2.293073136e-21
    Average vibrational energy      = 2.327098868e-20
    Average electronic energy       = 0
    Average total energy            = 2.02074239e-19
ExecutionTime = 147.6 s  ClockTime = 67 s

ExecutionTime = 151.33 s  ClockTime = 69 s

ExecutionTime = 155.49 s  ClockTime = 71 s

ExecutionTime = 159.56 s  ClockTime = 73 s

ExecutionTime = 163.65 s  ClockTime = 75 s

ExecutionTime = 167.73 s  ClockTime = 77 s

ExecutionTime = 171.79 s  ClockTime = 79 s

ExecutionTime = 175.74 s  ClockTime = 81 s

ExecutionTime = 180.1 s  ClockTime = 83 s

ExecutionTime = 183.98 s  ClockTime = 85 s

Time = 3.316213433e-06

    Collisions                      = 67665

Cloud name: dsmc
    Number of dsmc particles        = 1423986
    Number of molecules             = 8.543916e+15
    Mass in system                  = 3.97292094e-10
    Average linear kinetic energy   = 1.745903029e-19
    Average rotational energy       = 2.726957511e-21
    Average vibrational energy      = 2.327098649e-20
    Average electronic energy       = 0
    Average total energy            = 2.005882469e-19
ExecutionTime = 188.33 s  ClockTime = 87 s

ExecutionTime = 192.23 s  ClockTime = 89 s

ExecutionTime = 196.62 s  ClockTime = 91 s

ExecutionTime = 200.8 s  ClockTime = 93 s

ExecutionTime = 205.01 s  ClockTime = 95 s

ExecutionTime = 209.18 s  ClockTime = 97 s

ExecutionTime = 213.48 s  ClockTime = 99 s

ExecutionTime = 217.97 s  ClockTime = 99 s

ExecutionTime = 222.1 s  ClockTime = 101 s

ExecutionTime = 226.41 s  ClockTime = 103 s

Time = 3.979456119e-06

    Collisions                      = 77405

Cloud name: dsmc
    Number of dsmc particles        = 1447207
    Number of molecules             = 8.683242e+15
    Mass in system                  = 4.03770753e-10
    Average linear kinetic energy   = 1.727927302e-19
    Average rotational energy       = 3.190770424e-21
    Average vibrational energy      = 2.327101655e-20
    Average electronic energy       = 0
    Average total energy            = 1.992545172e-19
ExecutionTime = 230.88 s  ClockTime = 105 s

ExecutionTime = 235.01 s  ClockTime = 107 s

ExecutionTime = 239.37 s  ClockTime = 109 s

ExecutionTime = 243.48 s  ClockTime = 111 s

ExecutionTime = 247.95 s  ClockTime = 114 s

ExecutionTime = 252.37 s  ClockTime = 116 s

ExecutionTime = 257.08 s  ClockTime = 118 s

ExecutionTime = 262.4 s  ClockTime = 120 s

ExecutionTime = 267.11 s  ClockTime = 122 s

ExecutionTime = 271.49 s  ClockTime = 125 s

Time = 4.642698805e-06

    Collisions                      = 87100

Cloud name: dsmc
    Number of dsmc particles        = 1470678
    Number of molecules             = 8.824068e+15
    Mass in system                  = 4.10319162e-10
    Average linear kinetic energy   = 1.711331237e-19
    Average rotational energy       = 3.66084576e-21
    Average vibrational energy      = 2.327101398e-20
    Average electronic energy       = 0
    Average total energy            = 1.980649834e-19
ExecutionTime = 275.92 s  ClockTime = 127 s

ExecutionTime = 280.66 s  ClockTime = 129 s


Start sample at time 4.775347343e-06 and step 72

ExecutionTime = 285.07 s  ClockTime = 129 s

ExecutionTime = 289.57 s  ClockTime = 131 s

ExecutionTime = 294.09 s  ClockTime = 133 s

ExecutionTime = 298.88 s  ClockTime = 136 s

ExecutionTime = 303.7 s  ClockTime = 138 s

ExecutionTime = 308.5 s  ClockTime = 140 s

ExecutionTime = 313.17 s  ClockTime = 143 s

ExecutionTime = 317.81 s  ClockTime = 145 s

Time = 5.305941492e-06

    Collisions                      = 96664

Cloud name: dsmc
    Number of dsmc particles        = 1494051
    Number of molecules             = 8.964306e+15
    Mass in system                  = 4.16840229e-10
    Average linear kinetic energy   = 1.695364584e-19
    Average rotational energy       = 4.169580903e-21
    Average vibrational energy      = 2.327129187e-20
    Average electronic energy       = 0
    Average total energy            = 1.969773312e-19
ExecutionTime = 322.85 s  ClockTime = 148 s

ExecutionTime = 327.43 s  ClockTime = 150 s

ExecutionTime = 332.06 s  ClockTime = 152 s

Collision profiling summary:
    collision calls               = 83
    precompute candidates [s]     = 0.172884688
    rebuild partition [s]         = 0.011068006
    selection/collide [s]         = 11.10244779
    total profiled [s]            = 11.28640049


Sample average steps = 11

Total wall force: 
 wallForce_x 0.02384753499
 wallForce_y -0.006947281972
 wallForce_z -4.353206319e-05
Evolve profiling summary:
    evolve calls                  = 83
    move/build/coord [s]          = 134.1551294
    collision phase [s]           = 11.29030194
    reaction/output [s]           = 2.2372e-05
    post fields/output [s]        = 13.56880046
    total profiled [s]            = 159.0142542

ExecutionTime = 351.75 s  ClockTime = 168 s

End


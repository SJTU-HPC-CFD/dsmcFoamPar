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
Time   : 21:00:17
Host   : SuperXCX-ROG
PID    : 2370764
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
OpenMP enabled for dsmcCloud with 8 thread-local RNG streams using collision strategy 'dynamic'

Constructing constant properties for
    N2

Creating dsmcReactions

There are no chemical reactions defined.
Selecting collisionPartnerSelectionModel noTimeCounter

Initialising the measurement fields

Initialising dsmcVolFields field

Starting time loop

ExecutionTime = 9.83 s  ClockTime = 7 s

ExecutionTime = 13.05 s  ClockTime = 9 s

ExecutionTime = 16.24 s  ClockTime = 10 s

ExecutionTime = 19.17 s  ClockTime = 12 s

ExecutionTime = 22.46 s  ClockTime = 13 s

ExecutionTime = 25.43 s  ClockTime = 15 s

ExecutionTime = 28.51 s  ClockTime = 16 s

ExecutionTime = 31.7 s  ClockTime = 18 s

ExecutionTime = 34.97 s  ClockTime = 19 s

Time = 6.632426865e-07

    Collisions                      = 27937

Cloud name: dsmc
    Number of dsmc particles        = 1329945
    Number of molecules             = 7.97967e+15
    Mass in system                  = 3.71054655e-10
    Average linear kinetic energy   = 1.836335884e-19
    Average rotational energy       = 1.339817587e-21
    Average vibrational energy      = 2.327085575e-20
    Average electronic energy       = 0
    Average total energy            = 2.082442617e-19
ExecutionTime = 38.47 s  ClockTime = 21 s

ExecutionTime = 41.4 s  ClockTime = 22 s

ExecutionTime = 44.55 s  ClockTime = 24 s

ExecutionTime = 47.79 s  ClockTime = 25 s

ExecutionTime = 51.17 s  ClockTime = 27 s

ExecutionTime = 54.6 s  ClockTime = 29 s

ExecutionTime = 57.87 s  ClockTime = 30 s

ExecutionTime = 61.23 s  ClockTime = 32 s

ExecutionTime = 64.67 s  ClockTime = 33 s

ExecutionTime = 68.09 s  ClockTime = 32 s

Time = 1.326485373e-06

    Collisions                      = 37193

Cloud name: dsmc
    Number of dsmc particles        = 1353419
    Number of molecules             = 8.120514e+15
    Mass in system                  = 3.77603901e-10
    Average linear kinetic energy   = 1.809960603e-19
    Average rotational energy       = 1.577342231e-21
    Average vibrational energy      = 2.327085575e-20
    Average electronic energy       = 0
    Average total energy            = 2.058442582e-19
ExecutionTime = 71.66 s  ClockTime = 34 s

ExecutionTime = 74.93 s  ClockTime = 36 s

ExecutionTime = 78.34 s  ClockTime = 37 s

ExecutionTime = 81.72 s  ClockTime = 39 s

ExecutionTime = 85.25 s  ClockTime = 41 s

ExecutionTime = 88.79 s  ClockTime = 42 s

ExecutionTime = 92.18 s  ClockTime = 44 s

ExecutionTime = 96.23 s  ClockTime = 46 s

ExecutionTime = 99.61 s  ClockTime = 48 s

ExecutionTime = 103.3 s  ClockTime = 49 s

Time = 1.98972806e-06

    Collisions                      = 47473

Cloud name: dsmc
    Number of dsmc particles        = 1376992
    Number of molecules             = 8.261952e+15
    Mass in system                  = 3.84180768e-10
    Average linear kinetic energy   = 1.7865642e-19
    Average rotational energy       = 1.893826246e-21
    Average vibrational energy      = 2.327085575e-20
    Average electronic energy       = 0
    Average total energy            = 2.03821102e-19
ExecutionTime = 107.26 s  ClockTime = 51 s

ExecutionTime = 110.43 s  ClockTime = 53 s

ExecutionTime = 114.23 s  ClockTime = 55 s

ExecutionTime = 117.64 s  ClockTime = 57 s

ExecutionTime = 121.25 s  ClockTime = 58 s

ExecutionTime = 125.31 s  ClockTime = 60 s

ExecutionTime = 128.99 s  ClockTime = 62 s

ExecutionTime = 132.64 s  ClockTime = 64 s

ExecutionTime = 136.31 s  ClockTime = 63 s

ExecutionTime = 140.09 s  ClockTime = 65 s

Time = 2.652970746e-06

    Collisions                      = 58056

Cloud name: dsmc
    Number of dsmc particles        = 1400462
    Number of molecules             = 8.402772e+15
    Mass in system                  = 3.90728898e-10
    Average linear kinetic energy   = 1.765285643e-19
    Average rotational energy       = 2.290600206e-21
    Average vibrational energy      = 2.327105515e-20
    Average electronic energy       = 0
    Average total energy            = 2.020902197e-19
ExecutionTime = 144.39 s  ClockTime = 67 s

ExecutionTime = 147.85 s  ClockTime = 69 s

ExecutionTime = 152.01 s  ClockTime = 71 s

ExecutionTime = 155.54 s  ClockTime = 73 s

ExecutionTime = 159.31 s  ClockTime = 75 s

ExecutionTime = 163.16 s  ClockTime = 77 s

ExecutionTime = 167.24 s  ClockTime = 79 s

ExecutionTime = 170.96 s  ClockTime = 81 s

ExecutionTime = 174.98 s  ClockTime = 83 s

ExecutionTime = 178.79 s  ClockTime = 85 s

Time = 3.316213433e-06

    Collisions                      = 67566

Cloud name: dsmc
    Number of dsmc particles        = 1424012
    Number of molecules             = 8.544072e+15
    Mass in system                  = 3.97299348e-10
    Average linear kinetic energy   = 1.745954299e-19
    Average rotational energy       = 2.732101026e-21
    Average vibrational energy      = 2.327092112e-20
    Average electronic energy       = 0
    Average total energy            = 2.005984521e-19
ExecutionTime = 183.15 s  ClockTime = 87 s

ExecutionTime = 186.79 s  ClockTime = 89 s

ExecutionTime = 190.94 s  ClockTime = 91 s

ExecutionTime = 194.66 s  ClockTime = 93 s

ExecutionTime = 198.66 s  ClockTime = 92 s

ExecutionTime = 202.6 s  ClockTime = 94 s

ExecutionTime = 206.66 s  ClockTime = 96 s

ExecutionTime = 210.47 s  ClockTime = 98 s

ExecutionTime = 214.72 s  ClockTime = 100 s

ExecutionTime = 218.75 s  ClockTime = 102 s

Time = 3.979456119e-06

    Collisions                      = 77514

Cloud name: dsmc
    Number of dsmc particles        = 1447325
    Number of molecules             = 8.68395e+15
    Mass in system                  = 4.03803675e-10
    Average linear kinetic energy   = 1.727941233e-19
    Average rotational energy       = 3.199653939e-21
    Average vibrational energy      = 2.327098438e-20
    Average electronic energy       = 0
    Average total energy            = 1.992647616e-19
ExecutionTime = 223.09 s  ClockTime = 105 s

ExecutionTime = 227.42 s  ClockTime = 107 s

ExecutionTime = 231.73 s  ClockTime = 109 s

ExecutionTime = 236.6 s  ClockTime = 112 s

ExecutionTime = 241.51 s  ClockTime = 115 s

ExecutionTime = 246.25 s  ClockTime = 117 s

ExecutionTime = 250.62 s  ClockTime = 119 s

ExecutionTime = 254.7 s  ClockTime = 122 s

ExecutionTime = 259.15 s  ClockTime = 124 s

ExecutionTime = 263.22 s  ClockTime = 123 s

Time = 4.642698805e-06

    Collisions                      = 87439

Cloud name: dsmc
    Number of dsmc particles        = 1470744
    Number of molecules             = 8.824464e+15
    Mass in system                  = 4.10337576e-10
    Average linear kinetic energy   = 1.711266326e-19
    Average rotational energy       = 3.6751094e-21
    Average vibrational energy      = 2.327126714e-20
    Average electronic energy       = 0
    Average total energy            = 1.980730091e-19
ExecutionTime = 267.41 s  ClockTime = 125 s

ExecutionTime = 271.71 s  ClockTime = 128 s


Start sample at time 4.775347343e-06 and step 72

ExecutionTime = 275.85 s  ClockTime = 130 s

ExecutionTime = 280.05 s  ClockTime = 132 s

ExecutionTime = 284.38 s  ClockTime = 134 s

ExecutionTime = 288.61 s  ClockTime = 136 s

ExecutionTime = 293.06 s  ClockTime = 139 s

ExecutionTime = 297.22 s  ClockTime = 141 s

ExecutionTime = 301.52 s  ClockTime = 143 s

ExecutionTime = 305.75 s  ClockTime = 145 s

Time = 5.305941492e-06

    Collisions                      = 96744

Cloud name: dsmc
    Number of dsmc particles        = 1494088
    Number of molecules             = 8.964528e+15
    Mass in system                  = 4.16850552e-10
    Average linear kinetic energy   = 1.695632475e-19
    Average rotational energy       = 4.172903337e-21
    Average vibrational energy      = 2.327119841e-20
    Average electronic energy       = 0
    Average total energy            = 1.970073492e-19
ExecutionTime = 310.16 s  ClockTime = 148 s

ExecutionTime = 314.62 s  ClockTime = 150 s

ExecutionTime = 318.96 s  ClockTime = 152 s

Collision profiling summary:
    collision calls               = 83
    precompute candidates [s]     = 0.216618341
    rebuild partition [s]         = 2.5743e-05
    selection/collide [s]         = 9.346432881
    total profiled [s]            = 9.563076965


Sample average steps = 11

Total wall force: 
 wallForce_x 0.02373118164
 wallForce_y -0.006896768357
 wallForce_z 5.521424658e-05
Evolve profiling summary:
    evolve calls                  = 83
    move/build/coord [s]          = 135.6382817
    collision phase [s]           = 9.564011592
    reaction/output [s]           = 2.107e-05
    post fields/output [s]        = 13.82691264
    total profiled [s]            = 159.029227

ExecutionTime = 338.12 s  ClockTime = 168 s

End


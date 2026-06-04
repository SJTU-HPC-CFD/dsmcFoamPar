1. 激活环境只需要 source /home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/docs/env.sh
   编译只需要 source /home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/docs/compile.sh
2. replicated mesh模式，不需要decomposePar，运行不需要-parallel
3. 保留完整dsmcFoam+工作日志log输出
4. 每阶段工作的进展、测试数据结果、分析，都记录在/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/docs/detail中
5. 所有测试都要保证正确性（粒子数、碰撞数、能量一致）
0. 这轮工作的参照是我之前做的（在OF-v2506版本）dsmcFoam+的性能优化工作，主要是针对mpi+omp混合并行、omp-负载平衡dlb、mpi-负载平衡dlb（replicated mesh）的工作，需要你在OFv1706版本上复现之前的工作，并且在此基础上进行进一步的优化和测试，最终目标是实现dsmcFoam+在mpi+omp混合并行、omp-负载平衡dlb、mpi-负载平衡dlb（replicated mesh）三种模式下的性能提升，并且保证结果的正确性。
1. 激活环境只需要 source /home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb/doc/scripts/env.sh
   编译只需要 source /home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb/doc/scripts/build-dsmcFoam.sh
2. 保留完整dsmcFoam+工作日志log输出到对应case文件夹，如果有批量测试，则建立子文件夹存放log输出和对应controlDict
3. 每阶段工作的进展、测试数据结果、性能分析，都记录在/home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb/doc/worklog/v2506/detail中
4. 所有测试都要保证正确性（粒子数、碰撞数、能量一致）
5. 测试算例在/home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb/run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh，包含了原始的decomposePar划分的mpi8核算例，和你后需要开展的omp8核、mpi优化（replicated mesh）8核算例，不需要初始化，直接运行即可，如果你不慎修改了算例文件，可以从/home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb/run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmeshbkp中找到原始的算例文件进行替换
6. 性能基准在/home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb/run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmeshbkp对应的log日志中，你可以进行对比分析。
7. 代码参考路径是/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx，不要修改它
8. 工作顺序是: a.omp负载平衡dlb（dlb-omp）优化和测试，b.mpi负载平衡dlb（dlb-mpi）优化和测试，c.mpi+omp混合并行优化和测试，d.性能分析和总结。
9. 阅读我之前的工作日志：/home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb/doc/worklog/v2506，了解完整工作内容，但是在你后续的工作中，以参考代码/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx的实际内容为准
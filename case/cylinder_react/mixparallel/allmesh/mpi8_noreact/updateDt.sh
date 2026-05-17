#!/bin/sh

# 定义文件路径
newdtfile="log.dsmcInitialise+"
controlDict="./system/controlDict"
fieldPropertiesDict="./system/fieldPropertiesDict"

# 从newdtfile中读取step的值
deltaT=$(grep "deltaT" "$newdtfile" | awk '{print $2}')
steadyStateTime=$(grep "steadyStateTime" "$newdtfile" | awk '{print $2}')
totalTime=$(grep "totalTime" "$newdtfile" | awk '{print $2}')

# 检查是否成功读取到值
if [ -z "$deltaT" ]; then
    echo "错误：无法从log1中读取varA的值"
    exit 1
fi

# 使用sed替换controlDict中以"deltaT"开头的行中的值
sed -i "/^deltaT/s/deltaT.*/deltaT          $deltaT;/" "$controlDict"
# 使用sed替换controlDict中以"endTime"开头的行中的值
sed -i "/^endTime/s/endTime.*/endTime         $totalTime;/" "$controlDict"
# 使用sed替换controlDict中以"writeInterval"开头的行中的值
sed -i "/^writeInterval/s/writeInterval.*/writeInterval   $totalTime;/" "$controlDict"
# 使用sed替换fieldPropertiesDict中以"resetAtOutputUntilTime"开头的行中的值
sed -i "s/[[:space:]]*resetAtOutputUntilTime[[:space:]]\+[0-9.e-]\+/ 	          resetAtOutputUntilTime $steadyStateTime/g" "$fieldPropertiesDict"






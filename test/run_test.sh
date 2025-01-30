#!/bin/bash
parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:../external/mt-kahypar/build/lib
if [ "$OMPI_COMM_WORLD_RANK" == 0 ]; then
	exec ./test $*
else
	exec ./test $* >/dev/null 2>&1
fi
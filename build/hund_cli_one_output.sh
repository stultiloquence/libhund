#!/bin/bash
if [ "$OMPI_COMM_WORLD_RANK" == 0 ]; then
	exec $(dirname "${BASH_SOURCE[0]}")/hundcli $*
else
	exec $(dirname "${BASH_SOURCE[0]}")/hundcli $* >/dev/null 2>&1
fi
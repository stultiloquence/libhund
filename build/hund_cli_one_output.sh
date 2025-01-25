#!/bin/bash
if [ "$OMPI_COMM_WORLD_RANK" == 0 ]; then
	exec ./hundcli $*
else
	exec ./hundcli $* >/dev/null 2>&1
fi
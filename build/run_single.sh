#!/bin/bash
if [ "$OMPI_COMM_WORLD_RANK" == 0 ]; then
	exec ./runhund $*
else
	exec ./runhund $* >/dev/null 2>&1
fi
#!/bin/bash

#AUTHORS: Hardip Patel

# This utility runs a provided command, logs its execution details, 
# and can clean up specified files post-execution.

#TODO: Script has issues when the command to be run contains double quotes
#TODO: Script has issues when output has to be redirected to a file using '>' operator
#Both issues can be resolved by checking the command string and it does have these characteristics, then create a temp bash file and run the command from there.

# Initialize variables to hold argument values
command=""
taskname=""
donefile=""
force=false

# Use getopts for option parsing
while getopts ":c:t:d:f:" opt; do
  case $opt in
    c) # Command to execute
      command="$OPTARG"
      ;;
    t) # Task name
      taskname="$OPTARG"
      ;;
    d) # File to log execution details
      donefile="$OPTARG"
      ;;
    f) # Boolean force flag
      force="$OPTARG"
      if [[ "$force" != "true" && "$force" != "false" ]]; then
        echo "Error: -f option accepts only 'true' or 'false' as arguments."
        exit 1
      fi
      ;;
    \?) # Invalid option
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :) # Option without its required argument
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

# Ensure all required parameters are set
if [[ -z "$command" || -z "$taskname" || -z "$donefile" ]]; then
  echo "Usage: $0 -c <command> -t <task_name> -d <donefile_path> -f <true|false>"
  exit 1
fi

##convert command string to array
command=(${command})

EXECHOST=$(hostname)
JOBID=${PBS_JOBID}

# Determine if the task is already running
running=$(dirname ${donefile})/$(basename ${donefile} .done).running

if [[ -e "$running" ]]; then
	exit 0
fi

# Run the task if conditions require
if [[ "$force" == "true" ||  ! -e "$donefile"  || ! -s "$donefile" || "`tail -n1 $donefile`" != "EXIT_STATUS: 0" ]]; then
	touch $running
  echo COMMAND: "${command[@]}" >> $donefile
	echo EXECHOST: "$EXECHOST" >> $donefile
	echo JOBID: "$JOBID" >> $donefile
	echo TASKNAME: "$taskname" >>$donefile
	echo STARTTIME: $(date +%s) >>$donefile
    /usr/bin/time --format='ELAPSED: %e\nCPU: %S\nUSER: %U\nCPUPERCENT: %P\nMAXRM: %M Kb\nAVGRM: %t Kb\nAVGTOTRM: %K Kb\nPAGEFAULTS: %F\nRPAGEFAULTS: %R\nSWAP: %W\nWAIT: %w\nFSI: %I\nFSO: %O\nSMI: %r\nSMO: %s\nEXITCODE: %x' -o $donefile -a -- "${command[@]}"
    ret=$?
    echo ENDTIME: $(date +%s) >>$donefile
    echo EXIT_STATUS: $ret >>$donefile
	rm -f $running
    if [ "$ret" -ne 0 ]; then
        echo ERROR_command: "${command[@]}"
        echo ERROR_exitcode: $taskname failed with $ret exit code.
        exit $ret
    fi
fi

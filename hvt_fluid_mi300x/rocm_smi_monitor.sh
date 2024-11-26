#!/bin/bash

# Check if program and arguments are provided
if [ "$#" -lt 2 ]; then
  echo "Usage: $0 <program> <output_log_file> [program_arguments...]"
  exit 1
fi

# Extract program, log file, and program arguments
PROGRAM=$1
LOGFILE=$2
shift 2 # Remove program and log file arguments, leave the rest as program args

# Define a temporary file for monitoring script PID
MONITOR_PID_FILE=$(mktemp)

# Start the GPU monitoring with watch and save output to the log file
start=`date +%s`

(while true; do rocm-smi; timestamp=`date +%s`; echo "Time: "$((timestamp-start))"s"; sleep 1; done) > "$LOGFILE" &
MONITOR_PID=$!
echo "$MONITOR_PID" > "$MONITOR_PID_FILE"

# Run the target program
"$PROGRAM" "$@"

# After the program finishes, stop the monitoring
kill "$(cat "$MONITOR_PID_FILE")" 2>/dev/null

# Clean up
rm -f "$MONITOR_PID_FILE"

echo "GPU monitoring has stopped. Results are saved in $LOGFILE."

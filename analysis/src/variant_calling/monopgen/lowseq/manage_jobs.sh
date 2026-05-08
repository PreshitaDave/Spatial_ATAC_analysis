#!/bin/bash

###############################################################################
# Job Management Script for Monopogen Cluster Jobs
# Usage: ./manage_jobs.sh [command]
###############################################################################

echo "=========================================="
echo "Monopogen Job Manager"
echo "=========================================="
echo ""

case "${1:-status}" in
    
    status|list)
        echo "[1] CURRENT JOB STATUS"
        echo "---------"
        qstat -u $USER 2>/dev/null | grep "mono_"
        echo ""
        echo "[2] JOB COUNT"
        echo "---------"
        JOB_COUNT=$(qstat -u $USER 2>/dev/null | grep "mono_" | wc -l)
        echo "Active mono_* jobs: $JOB_COUNT"
        echo ""
        ;;
        
    cancel-all)
        echo "⚠️  WARNING: This will cancel ALL running jobs"
        read -p "Continue? (yes/no): " confirm
        if [ "$confirm" = "yes" ]; then
            echo "Canceling all jobs..."
            qstat -u $USER 2>/dev/null | grep "mono_" | awk '{print $1}' | xargs -I {} qdel {}
            echo "✓ All jobs canceled"
        else
            echo "Canceled (no jobs deleted)"
        fi
        echo ""
        ;;
        
    cancel-chr)
        CHR=$2
        if [ -z "$CHR" ]; then
            echo "ERROR: Please specify chromosome"
            echo "Usage: ./manage_jobs.sh cancel-chr <chromosome>"
            echo "Example: ./manage_jobs.sh cancel-chr chr8"
            exit 1
        fi
        echo "Finding jobs for $CHR..."
        qstat -u $USER 2>/dev/null | grep "mono_${CHR}" | awk '{print $1}' | while read job_id; do
            echo "  - Canceling job $job_id..."
            qdel "$job_id"
        done
        echo "✓ Jobs for $CHR canceled"
        echo ""
        ;;
        
    submit-one)
        CHR=$2
        if [ -z "$CHR" ]; then
            echo "ERROR: Please specify chromosome"
            echo "Usage: ./manage_jobs.sh submit-one <chromosome>"
            echo "Example: ./manage_jobs.sh submit-one chr8"
            exit 1
        fi
        SCRIPT="qsub_${CHR}.sh"
        if [ ! -f "$SCRIPT" ]; then
            echo "ERROR: $SCRIPT not found"
            exit 1
        fi
        echo "Submitting: $SCRIPT"
        job_id=$(qsub "$SCRIPT")
        echo "✓ Job submitted: $job_id"
        echo ""
        ;;
        
    submit-all)
        echo "Submitting all qsub_*.sh scripts..."
        SUBMITTED=0
        for script in qsub_chr*.sh; do
            if [ -f "$script" ]; then
                echo "  - Submitting $script"
                job_id=$(qsub "$script")
                echo "    Job ID: $job_id"
                ((SUBMITTED++))
            fi
        done
        echo ""
        echo "✓ Submitted $SUBMITTED jobs"
        echo ""
        ;;
        
    *)
        cat << 'HELP'

COMMANDS:
  status (or no argument)     - Show current running jobs
  list                        - Show current running jobs
  
  cancel-all                  - Cancel ALL mono_* jobs (interactive)
  cancel-chr <chr>           - Cancel jobs for specific chromosome
  
  submit-one <chr>           - Submit single chromosome job
  submit-all                 - Submit ALL qsub_chr*.sh scripts
  

EXAMPLES:
  ./manage_jobs.sh                      # Show job status
  ./manage_jobs.sh cancel-chr chr8      # Cancel chr8 jobs
  ./manage_jobs.sh submit-one chr8      # Submit only chr8
  ./manage_jobs.sh submit-all           # Submit all chromosomes

QUICK TIPS:
  • Check job queue:     qstat -u $USER
  • Check specific job:  qstat -j <job_id>
  • Manually delete:     qdel <job_id>
  • Delete all:          qdel '*'  (be careful!)

HELP
        ;;
        
esac

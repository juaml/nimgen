TEMPLATE_DAG = """
JOB  A  step_1.submit
JOB  B  step_2.submit
JOB  C  step_3.submit
JOB  D  step_4.submit

PARENT A CHILD B
PARENT B CHILD C
PARENT C CHILD D
"""

TEMPLATE_JOB = """
# The environment
transfer_executable = False
initial_dir={}
executable = $(initial_dir)/run_in_venv.sh
universe = vanilla
getenv = True

# Resources
request_cpus    = {}
request_memory  = {}
request_disk    = {}

"""

TEMPLATE_QUEUED_JOB = """
arguments={}

log = {logs}/{job_id}$(Cluster).$(Process).log
output = {logs}/{job_id}$(Cluster).$(Process).out
error = {logs}/{job_id}$(Cluster).$(Process).err

queue
"""

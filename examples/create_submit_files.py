from itertools import product
import os
from pathlib import Path
import random

TEMPLATE_DAG = """
JOB  A  {}
JOB  B  {}
JOB  C  {}
JOB  D  {}

PARENT A CHILD B 
PARENT B CHILD C
PARENT C CHILD D
"""

TEMPLATE_JOB = """
# The environment
executable = /usr/bin/bash
transfer_executable = False
initial_dir={}
universe = vanilla
getenv = True

# Resources
request_cpus    = {}
request_memory  = {}
request_disk    = {}

"""

TEMPLATE_QUEUED_JOB = """
arguments={} 

log = {logs}/$(Cluster).$(Process).log
output = {logs}/$(Cluster).$(Process).out
error = {logs}/$(Cluster).$(Process).err

queue
"""


def submission_job1(
    parcel_file,    
    marker_file,    
    n_surrogate_maps,
    submission_file='job1.submit',
    initial_dir="../../compute_gene_lists/",
    request_cpus=8,
    request_memory=64,
    request_disk=100,
):

    job1 = TEMPLATE_JOB.format(
        initial_dir, request_cpus, request_memory, request_disk
    )

    with open(submission_file, "w") as f:
        f.write(job1)

    arguments = f"./run_in_venv.sh ./job1.py {parcel_file} {marker_file} {n_surrogate_maps}"
    submit_dir = os.path.dirname(submission_file)
    with open(submission_file, "a") as f:
        f.write(TEMPLATE_QUEUED_JOB.format(arguments, logs=f"{submit_dir}/logs"))

    
def submission_job2(
    parcel_file,
    marker_file,
    n_surrogate_maps=5,
    submission_file='job2.submit',
    initial_dir="../../compute_gene_lists/",
    request_cpus=1,
    request_memory=64,
    request_disk=10,
):

    job2 = TEMPLATE_JOB.format(
        initial_dir, request_cpus, request_memory, request_disk
    )

    with open(submission_file, "w") as f:
        f.write(job2)

    arguments = [
        f"./run_in_venv.sh ./job2.py {parcel_file} {marker_file} {smap}" for smap in range(n_surrogate_maps)
    ]

    submit_dir = os.path.dirname(submission_file)
    for arg in arguments:
        with open(submission_file, "a") as f:
            f.write(TEMPLATE_QUEUED_JOB.format(arg, logs=f"{submit_dir}/logs" ))


def submission_job3(
    parcel_file,
    marker_file,
    submission_file='job3.submit',
    initial_dir="../../compute_gene_lists/",
    request_cpus=1,
    request_memory=64,
    request_disk=10,
):
    job3 = TEMPLATE_JOB.format(
        initial_dir, request_cpus, request_memory, request_disk
    )

    with open(submission_file, "w") as f:
        f.write(job3)

    arguments =  f"./run_in_venv.sh ./job3.py {parcel_file} {marker_file}"

    submit_dir = os.path.dirname(submission_file)
    with open(submission_file, "a") as f:
        f.write(TEMPLATE_QUEUED_JOB.format(arguments, logs=f"{submit_dir}/logs" ))


def submission_job4(
    parcel_file,
    marker_file,
    partial=True, pca=True, n_comps=[1, 3, 5, 10], custom_covariates=False,
    submission_file='job4.submit',
    initial_dir="../../compute_gene_lists/",
    request_cpus=1,
    request_memory=64,
    request_disk=10,
):
    job3 = TEMPLATE_JOB.format(
        initial_dir, request_cpus, request_memory, request_disk
    )

    with open(submission_file, "w") as f:
        f.write(job3)

    arguments = [ 
        f"./run_in_venv.sh ./job4.py {parcel_file} {marker_file} partial_correlation=True perform_pca=True n_comps={comp} " for comp in n_comps
     ]

    submit_dir = os.path.dirname(submission_file)
    for arg in arguments:
        with open(submission_file, "a") as f:
            f.write(TEMPLATE_QUEUED_JOB.format(arg, logs=f"{submit_dir}/logs" ))


def dag_script(dag_name, job1, job2, job3, job4):
    dag_script = TEMPLATE_DAG.format(
        job1, job2, job3, job4
    )

    with open(dag_name, 'w') as f:
        f.write(dag_script)


def run_dag_script(project_dir, parcels_and_markers, n_surrogate_maps, 
partial=True, pca=True, n_comps=[1, 3, 5], custom_covariates=False,
submit_files_dir='code/submit_files'):

    for pm in parcels_and_markers:
        parcel_file = os.path.join(project_dir, pm[0])
        marker_file = os.path.join(project_dir, pm[1])
        parcel_filename = Path(pm[0]).stem.split('.')[0]
        marker_filename = Path(pm[1]).stem.split('.')[0]
        submit_dir = os.path.join( project_dir, submit_files_dir, f'{marker_filename}_{parcel_filename}_{n_surrogate_maps}')
        submit_log_dir = os.path.join(submit_dir, 'logs')
        if not os.path.isdir(submit_dir) : os.makedirs(submit_dir)
        if not os.path.isdir(submit_log_dir) : os.makedirs(submit_log_dir)

        submission_job1(
            parcel_file, marker_file, n_surrogate_maps,
            initial_dir=os.path.join(project_dir,'code','compute_gene_lists'),
            submission_file= os.path.join(submit_dir, 'job1.submit')
        )

        submission_job2(
            parcel_file,
            marker_file,
            n_surrogate_maps,
            initial_dir=os.path.join(project_dir,'code','compute_gene_lists'),
            submission_file= os.path.join(submit_dir, 'job2.submit'),
            request_cpus=1,
            request_memory=1,
            request_disk=1,
        )

        submission_job3(
            parcel_file,
            marker_file,
            initial_dir=os.path.join(project_dir,'code','compute_gene_lists'),
            submission_file= os.path.join(submit_dir, 'job3.submit'),
            request_cpus=1,
            request_memory=32,
            request_disk=10,
        )


        submission_job4(
            parcel_file,
            marker_file,
            partial=True, pca=True, n_comps=[1, 3, 5], custom_covariates=False,
            initial_dir=os.path.join(project_dir,'code','compute_gene_lists'),
            submission_file= os.path.join(submit_dir, 'job4.submit'),
            request_cpus=1,
            request_memory=32,
            request_disk=10,
        )


        dag_script(
            f"{submit_dir}/nimgen_dag.dag",
            "job1.submit",
            "job2.submit",
            "job3.submit",
            "job4.submit"
        )


if __name__ == "__main__":

    # tuple (parcel,marker)
    parcels_and_markers = [
        ("input/parcellations/Resliced_FinalAtlas_436_FixedAffine/Resliced_FinalAtlas_436_FixedAffine.nii", "input/markers/HCP/FC/436Parcels_HCP1_acc_v18.nii"),

        ("input/parcellations/DK_RESAMPLED_83_vAbagen/DK_RESAMPLED_83_vAbagen.nii.gz", "input/markers/HCP/VBM/S1200_AverageT1wDividedByT2w.nii.gz"),

        ("input/parcellations/Schaefer2018_100Parcels_17Networks_order_FSLMNI152_2mm/Schaefer2018_100Parcels_17Networks_order_FSLMNI152_2mm.nii.gz", "input/markers/HCP/VBM/S1200_AverageT1wDividedByT2w.nii.gz"),
        ]

    project_dir = '/data/project/nimgen/nimgen_validation_study'
    run_dag_script(
        project_dir, parcels_and_markers, 1000, partial=True, pca=True, n_comps=[1, 3, 5, 10], custom_covariates=False, submit_files_dir = "submit_files/"
    )



